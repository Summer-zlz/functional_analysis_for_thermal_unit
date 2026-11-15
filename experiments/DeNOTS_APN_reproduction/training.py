"""Three rolling validation folds, validation-only selection and final 80/20 fit."""
import json
from itertools import product
from pathlib import Path
import random
import time
import bootstrap
import numpy as np
import pandas as pd
import torch
from data_loader import Standardizer, batches, collate, load_dataset
from evaluate import metrics
from official_models import make_model

CONFIG = {"seed": 42, "hidden_dimensions": [32, 64], "learning_rate": .001,
          "optimizer": "Adam", "loss": "MSE of standardized response at valid observations",
          "max_epochs": 100, "patience": 10, "batch_size": 16, "gradient_clip": 5.,
          "denots_depth": 10, "denots_tolerance": .001, "denots_vf": "adaptive", "denots_nf": True,
          "apn_patches": 8, "apn_time_dimension": 16, "apn_dropout": .1,
          "final_epochs": "rounded median of selected candidate's three best validation epochs"}


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(4)


def device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def candidate_parameters(model_name):
    grid = {"hidden_dim": CONFIG["hidden_dimensions"]}
    if model_name == "DeNOTS" and "denots_depth_grid" in CONFIG:
        grid["denots_depth"] = CONFIG["denots_depth_grid"]
    if model_name == "APN" and "apn_patches_grid" in CONFIG:
        grid["apn_patches"] = CONFIG["apn_patches_grid"]
        grid["apn_time_dimension"] = CONFIG["apn_time_dimensions_grid"]
    return [dict(zip(grid, values)) for values in product(*grid.values())]


def candidate_key(parameters):
    key = f"width{parameters['hidden_dim']}"
    for name, short in (("denots_depth", "depth"), ("apn_patches", "patch"),
                        ("apn_time_dimension", "te")):
        if name in parameters:
            key += f"_{short}{parameters[name]}"
    return key


def predict(model, trajectories, scaler, chosen_device, model_name, fold):
    model.eval()
    rows = []
    with torch.no_grad():
        for samples in batches(trajectories, CONFIG["batch_size"]):
            batch = collate(samples, scaler, chosen_device)
            output = model(batch["x"], batch["time"], batch["mask"])
            assert output.shape == batch["mask"].shape
            output = output.cpu().numpy().astype(np.float64)*scaler.y_std+scaler.y_mean
            for i, sample in enumerate(samples):
                for k, t in enumerate(sample.time):
                    rows.append({"model": model_name, "dataset": sample.dataset, "fold": fold,
                                 "id": sample.id, "time": t, "actual": sample.y[k],
                                 "predicted": output[i, k], "error": output[i, k]-sample.y[k]})
    answer = pd.DataFrame(rows)
    assert len(answer) == sum(len(s.time) for s in trajectories)
    assert np.isfinite(answer[["actual", "predicted", "error"]].to_numpy()).all()
    return answer


def train_fit(model_name, parameters, training, validation, scaler, folder, epochs=None):
    set_seed(CONFIG["seed"])
    chosen_device = device()
    model_parameters = {**parameters, "input_dim": len(scaler.x_mean)}
    model = make_model(model_name, **model_parameters).to(chosen_device)
    optimizer = torch.optim.Adam(model.parameters(), lr=CONFIG["learning_rate"])
    rng = np.random.default_rng(CONFIG["seed"])
    folder.mkdir(parents=True, exist_ok=True)
    scaler.save(folder / "standardizer.json")
    best_rmse, best_epoch, stale, best_state = np.inf, 0, 0, None
    history, duration = [], time.perf_counter()
    limit = CONFIG["max_epochs"] if epochs is None else epochs
    for epoch in range(1, limit+1):
        model.train()
        total_loss, total_points = 0., 0
        for samples in batches(training, CONFIG["batch_size"], shuffle=True, rng=rng):
            batch = collate(samples, scaler, chosen_device)
            optimizer.zero_grad(set_to_none=True)
            predicted = model(batch["x"], batch["time"], batch["mask"])
            loss = (predicted[batch["mask"]]-batch["y"][batch["mask"]]).square().mean()
            if not torch.isfinite(loss):
                raise RuntimeError("Nonfinite training loss")
            loss.backward()
            norm = torch.nn.utils.clip_grad_norm_(model.parameters(), CONFIG["gradient_clip"], error_if_nonfinite=True)
            optimizer.step()
            points = int(batch["mask"].sum().item())
            total_loss += loss.item()*points
            total_points += points
        record = {"epoch": epoch, "standardized_train_MSE": total_loss/total_points,
                  "gradient_norm_before_clip": norm.item()}
        if validation is not None:
            predictions = predict(model, validation, scaler, chosen_device, model_name, "validation")
            score = metrics(predictions.actual, predictions.predicted)["RMSE"]
            record["validation_RMSE"] = score
            if score < best_rmse-1e-8:
                best_rmse, best_epoch, stale = score, epoch, 0
                best_state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
            else:
                stale += 1
        else:
            best_epoch = epoch
        history.append(record)
        pd.DataFrame(history).to_csv(folder / "history.csv", index=False)
        if epoch == 1 or epoch % 10 == 0:
            print(f"{model_name} {training[0].dataset} {folder.parent.name}/{folder.name} "
                  f"epoch={epoch} loss={record['standardized_train_MSE']:.6f} "
                  f"val_RMSE={record.get('validation_RMSE', float('nan')):.6f}", flush=True)
        if validation is not None and stale >= CONFIG["patience"]:
            break
    if validation is not None:
        assert best_state is not None
        model.load_state_dict(best_state)
    elapsed = time.perf_counter()-duration
    torch.save({"model": model_name, "hidden_dim": parameters["hidden_dim"],
                "model_parameters": model_parameters, "state_dict": model.state_dict(),
                "standardizer": scaler.to_dict(), "configuration": CONFIG,
                "best_epoch": best_epoch, "completed_epochs": epoch}, folder / "checkpoint.pt")
    result = {**parameters, "input_dim": len(scaler.x_mean),
              "best_epoch": best_epoch, "completed_epochs": epoch,
              "training_seconds": elapsed, "device": str(chosen_device), "status": "complete"}
    (folder / "training_status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return model, result


def run_dataset(model_name, label, results):
    samples, splits, names = load_dataset(label)
    folder = Path(results) / model_name / label
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "feature_names.json").write_text(json.dumps(names, indent=2), encoding="utf-8")
    split_rows = [{"dataset": label, "split": role, "id": identity}
                  for role, identities in splits.items() for identity in identities]
    pd.DataFrame(split_rows).to_csv(folder / "split_ids.csv", index=False)
    split_summary = [{"dataset": label, "split": role, "trajectories": len(ids),
        "observations": sum(len(samples[i].time) for i in ids), "id_min": min(ids), "id_max": max(ids)}
        for role, ids in splits.items()]
    pd.DataFrame(split_summary).to_csv(folder / "split_summary.csv", index=False)
    candidates = []
    parameter_grid = candidate_parameters(model_name)
    parameter_columns = list(parameter_grid[0])
    for parameters in parameter_grid:
        key = candidate_key(parameters)
        for fold in (1, 2, 3):
            train = [samples[i] for i in splits[f"fold{fold}_train"]]
            val = [samples[i] for i in splits[f"fold{fold}_validation"]]
            subfolder = folder / key / f"fold{fold}"
            print(f"Starting {model_name} {label} {key} Fold {fold}", flush=True)
            try:
                scaler = Standardizer(train)
                model, status = train_fit(model_name, parameters, train, val, scaler, subfolder)
                prediction = predict(model, val, scaler, device(), model_name, f"fold{fold}_validation")
                prediction.to_csv(subfolder / "validation_predictions.csv", index=False, float_format="%.17g")
                candidates.append({"model": model_name, "dataset": label, "fold": fold,
                    "candidate": key,
                    "trajectories": len(val), "observations": len(prediction),
                    **status, **metrics(prediction.actual, prediction.predicted)})
                del model
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            except Exception as error:
                candidates.append({"model": model_name, "dataset": label, "fold": fold,
                                   "candidate": key, **parameters,
                                   "status": "error", "message": str(error)})
                print(f"Excluded fold {key}/{fold}: {error}", flush=True)
            pd.DataFrame(candidates).to_csv(folder / "cv_results.csv", index=False)
    table = pd.DataFrame(candidates)
    summaries = []
    for parameters in parameter_grid:
        key = candidate_key(parameters)
        rows = table[table.candidate == key].sort_values("fold")
        complete = len(rows) == 3 and (rows.status == "complete").all()
        eligible = complete and np.isfinite(rows[["MAE", "RMSE", "MAPE"]].to_numpy()).all()
        row = {"model": model_name, "dataset": label, "candidate": key, **parameters,
               "eligible": bool(eligible), "status": "complete" if eligible else "excluded_error"}
        for k in (1, 2, 3):
            fold_rows = rows[rows.fold == k]
            row[f"fold{k}_RMSE"] = float(fold_rows.RMSE.iloc[0]) if "RMSE" in fold_rows else np.nan
        row["mean_validation_RMSE"] = float(rows.RMSE.mean()) if eligible else np.nan
        summaries.append(row)
    scores = pd.DataFrame(summaries)
    scores.to_csv(folder / "cv_candidate_summary.csv", index=False)
    scores = scores[scores.eligible].sort_values(["mean_validation_RMSE"]+parameter_columns)
    if scores.empty:
        raise RuntimeError(f"No complete finite three-fold candidate: {model_name} {label}")
    winner = scores.iloc[0]
    parameters = {name: int(winner[name]) for name in parameter_columns}
    selected = table[table.candidate == winner.candidate].sort_values("fold")
    assert len(selected) == 3 and np.isfinite(selected[["MAE", "RMSE", "MAPE"]].to_numpy()).all()
    selected.to_csv(folder / "selected_fold_metrics.csv", index=False)
    epochs = int(np.floor(np.median(selected.best_epoch)+.5))
    best = {"model": model_name, "dataset": label, "candidate": winner.candidate, **parameters,
            "mean_validation_RMSE": float(selected.RMSE.mean()), "final_epochs": epochs,
            "epoch_rule": CONFIG["final_epochs"], "tie_rule": "ascending "+", ".join(parameter_columns)}
    (folder / "best_parameters.json").write_text(json.dumps(best, indent=2), encoding="utf-8")
    train = [samples[i] for i in splits["final_train"]]
    test = [samples[i] for i in splits["final_test"]]
    scaler = Standardizer(train)
    model, status = train_fit(model_name, parameters, train, None, scaler, folder / "final", epochs=epochs)
    for period, trajectories in (("train", train), ("test", test)):
        predictions = predict(model, trajectories, scaler, device(), model_name, f"final_{period}")
        predictions.to_csv(folder / f"{period}_predictions.csv", index=False, float_format="%.17g")
    (folder / "run_status.json").write_text(json.dumps({**status, **best}, indent=2), encoding="utf-8")
    print(f"Completed {model_name} {label}; selected {winner.candidate}, final epochs={epochs}", flush=True)

