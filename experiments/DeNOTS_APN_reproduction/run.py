"""Standalone rolling validation and direct group fitting for DeNOTS/APN."""
import argparse
import json
from pathlib import Path
import bootstrap
import numpy as np
import pandas as pd
import torch
from data_loader import Trajectory, Standardizer, load_dataset
from evaluate import metrics
from official_models import make_model
from training import CONFIG, device, predict, run_dataset, set_seed, train_fit

LABELS = ("UP", "DOWN", "STABLE")
GROUPS = ("UP1", "UP2", "DOWN1", "DOWN2", "STABLE")


def configure():
    CONFIG.update({"hidden_dimensions": [8, 16], "denots_depth_grid": [1, 5],
        "apn_patches_grid": [1, 2], "apn_time_dimensions_grid": [4, 8],
        "seed": 42, "learning_rate": .001})


def group_data(group):
    if group == "STABLE":
        samples, _, _ = load_dataset(group)
        return group, [samples[i] for i in sorted(samples)]
    source, cluster = group[:-1], int(group[-1])
    data = pd.read_csv(bootstrap.PROJECT / f"coaldata_{source.lower()}_unit1.csv")
    assert (data.groupby("id").cluster.nunique() == 1).all()
    data = data[data.cluster == cluster].sort_values(["id", "time"], kind="stable")
    samples = [Trajectory(group, int(identity), rows.time.to_numpy(dtype=np.float64),
        rows.iloc[:, 2:-2].to_numpy(dtype=np.float64), rows.coalseries.to_numpy(dtype=np.float64))
        for identity, rows in data.groupby("id", sort=True)]
    return source, samples


def fit_direct(root):
    for name in ("DeNOTS", "APN"):
        for group in GROUPS:
            source, samples = group_data(group)
            best = json.loads((root / "rolling" / name / source / "best_parameters.json").read_text())
            params = {key: best[key] for key in ("hidden_dim", "denots_depth", "apn_patches",
                                                "apn_time_dimension") if key in best}
            folder = root / "direct_fit" / name / group
            folder.mkdir(parents=True, exist_ok=True)
            model, status = train_fit(name, params, samples, None, Standardizer(samples),
                                      folder / "fit", epochs=int(best["final_epochs"]))
            output = predict(model, samples, Standardizer(samples), device(), name, "all_trajectory_fit")
            output.to_csv(folder / "fitted_predictions.csv", index=False, float_format="%.17g")
            (folder / "run_status.json").write_text(json.dumps(status, indent=2))


def summarize(root):
    rows, selected = [], []
    for name in ("DeNOTS", "APN"):
        for label in LABELS:
            folder = root / "rolling" / name / label
            selected.append(json.loads((folder / "best_parameters.json").read_text()))
            for period in ("train", "test"):
                table = pd.read_csv(folder / f"{period}_predictions.csv")
                rows.append({"model": name, "dataset": label, "period": period,
                    "trajectories": table.id.nunique(), "observations": len(table),
                    **metrics(table.actual, table.predicted)})
        for group in GROUPS:
            table = pd.read_csv(root / "direct_fit" / name / group / "fitted_predictions.csv")
            rows.append({"model": name, "dataset": group, "period": "all_sample_fit",
                "trajectories": table.id.nunique(), "observations": len(table),
                **metrics(table.actual, table.predicted)})
    pd.DataFrame(rows).to_csv(root / "metrics.csv", index=False)
    pd.DataFrame(selected).to_csv(root / "best_parameters.csv", index=False)


def verify(root):
    set_seed()
    count = 0
    for name in ("DeNOTS", "APN"):
        for kind, labels in (("rolling", LABELS), ("direct_fit", GROUPS)):
            for label in labels:
                folder = root / kind / name / label
                if kind == "rolling":
                    data, splits, _ = load_dataset(label)
                    samples_by_period = {p: [data[i] for i in splits[f"final_{p}"]]
                                         for p in ("train", "test")}
                    checkpoint_path = folder / "final" / "checkpoint.pt"
                    saved_splits = pd.read_csv(folder / "split_ids.csv")
                    for role, ids in splits.items():
                        assert sorted(saved_splits.loc[saved_splits.split == role, "id"]) == ids
                    scores = pd.read_csv(folder / "cv_candidate_summary.csv")
                    best = json.loads((folder / "best_parameters.json").read_text())
                    keys = [k for k in ("hidden_dim", "denots_depth", "apn_patches", "apn_time_dimension") if k in scores]
                    assert scores[scores.eligible].sort_values(["mean_validation_RMSE"]+keys).iloc[0].candidate == best["candidate"]
                else:
                    _, samples = group_data(label)
                    samples_by_period = {"fitted": samples}
                    checkpoint_path = folder / "fit" / "checkpoint.pt"
                checkpoint = torch.load(checkpoint_path, map_location=device(), weights_only=False)
                model = make_model(name, **checkpoint["model_parameters"]).to(device())
                model.load_state_dict(checkpoint["state_dict"])
                scaler = Standardizer.from_dict(checkpoint["standardizer"])
                training = samples_by_period.get("train", samples_by_period.get("fitted"))
                assert scaler.to_dict() == Standardizer(training).to_dict()
                for period, samples in samples_by_period.items():
                    filename = "fitted_predictions.csv" if period == "fitted" else f"{period}_predictions.csv"
                    saved = pd.read_csv(folder / filename)
                    fresh = predict(model, samples, scaler, device(), name, "verify")
                    assert np.array_equal(saved.id, fresh.id)
                    for col in ("time", "actual"):
                        assert np.allclose(saved[col], fresh[col], rtol=0, atol=1e-12)
                    assert np.allclose(saved.predicted, fresh.predicted, rtol=0, atol=2e-5)
                    count += 1
    (root / "verification.json").write_text(json.dumps({"status": "passed",
        "replayed_prediction_files": count, "tolerance": 2e-5}, indent=2))
    print(f"Verified {count} prediction files from saved checkpoints.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("all", "rolling", "direct", "verify"), default="all")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    configure()
    root = args.output or bootstrap.HERE / "rerun_results"
    root.mkdir(parents=True, exist_ok=True)
    if args.stage in ("all", "rolling"):
        (root / "configuration.json").write_text(json.dumps(CONFIG, indent=2))
        for name in ("DeNOTS", "APN"):
            for label in LABELS:
                run_dataset(name, label, root / "rolling")
    if args.stage in ("all", "direct"):
        fit_direct(root)
    if args.stage in ("all", "verify"):
        verify(root)
    if args.stage in ("all", "direct"):
        summarize(root)
