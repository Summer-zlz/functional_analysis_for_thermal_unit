"""Whole-trajectory datasets with independently computed integer boundaries."""
from dataclasses import dataclass
import json
from pathlib import Path
import bootstrap
import numpy as np
import pandas as pd
import torch

PROJECT = bootstrap.PROJECT
@dataclass
class Trajectory:
    dataset: str
    id: int
    time: np.ndarray
    x: np.ndarray
    y: np.ndarray


def load_dataset(label):
    stable = label == "STABLE"
    data = pd.read_csv(PROJECT / f"coaldata_{label.lower()}_unit1.csv",
                       usecols=None if stable else range(20))
    assert data.shape[1] == (22 if stable else 20) and data.columns[-1] == "coalseries"
    assert list(data.columns[:3]) == ["id", "time", "mw"]
    assert np.isfinite(data.to_numpy()).all() and not data.duplicated(["id", "time"]).any()
    data = data.sort_values(["id", "time"], kind="stable")
    trajectories = {}
    for identity, group in data.groupby("id", sort=True):
        assert identity == int(identity)
        t = group.time.to_numpy(dtype=np.float64)
        assert len(t) >= 2 and np.all(np.diff(t) > 0)
        trajectories[int(identity)] = Trajectory(label, int(identity), t,
            group.iloc[:, 2:-1].to_numpy(dtype=np.float64),
            group.iloc[:, -1].to_numpy(dtype=np.float64))
    ids = np.array(sorted(trajectories))
    boundaries = {p: int(np.floor(len(ids)*p)) for p in (0, .1, .2, .5, .6, .7, .8, 1)}
    periods = {"fold1_train": (0, .5), "fold1_validation": (.5, .6),
               "fold2_train": (.1, .6), "fold2_validation": (.6, .7),
               "fold3_train": (.2, .7), "fold3_validation": (.7, .8),
               "final_train": (0, .8), "final_test": (.8, 1)}
    splits = {role: ids[boundaries[a]:boundaries[b]].tolist()
              for role, (a, b) in periods.items()}
    assert not set(splits["final_train"]) & set(splits["final_test"])
    for k in (1, 2, 3):
        assert not set(splits[f"fold{k}_train"]) & set(splits[f"fold{k}_validation"])
        assert set(splits[f"fold{k}_train"]+splits[f"fold{k}_validation"]) <= set(splits["final_train"])
    return trajectories, splits, list(data.columns[2:-1])


class Standardizer:
    def __init__(self, trajectories):
        x = np.concatenate([s.x for s in trajectories])
        y = np.concatenate([s.y for s in trajectories])
        self.x_mean = x.mean(axis=0)
        self.x_std = x.std(axis=0)
        self.x_std[self.x_std < 1e-12] = 1
        self.y_mean, self.y_std = float(y.mean()), float(y.std())
        assert self.y_std > 0
        self.time_scale = max(float(s.time[-1]) for s in trajectories)
        assert self.time_scale > 0
        self.training_ids = [s.id for s in trajectories]

    def to_dict(self):
        return {"x_mean": self.x_mean.tolist(), "x_std": self.x_std.tolist(),
                "y_mean": self.y_mean, "y_std": self.y_std,
                "time_scale": self.time_scale, "training_ids": self.training_ids}

    def save(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, values):
        result = cls.__new__(cls)
        for name, value in values.items():
            setattr(result, name, np.asarray(value) if name in ("x_mean", "x_std") else value)
        return result


def collate(trajectories, standardizer, device):
    count, length = len(trajectories), max(len(s.time) for s in trajectories)
    x = np.zeros((count, length, len(standardizer.x_mean)), dtype=np.float32)
    t = np.zeros((count, length), dtype=np.float32)
    y = np.zeros((count, length), dtype=np.float32)
    mask = np.zeros((count, length), dtype=bool)
    for i, sample in enumerate(trajectories):
        m = len(sample.time)
        x[i, :m] = (sample.x-standardizer.x_mean)/standardizer.x_std
        # Repeating last x/time fills storage only; mask excludes all padding.
        x[i, m:] = x[i, m-1]
        t[i, :m] = sample.time/standardizer.time_scale
        t[i, m:] = t[i, m-1]
        y[i, :m] = (sample.y-standardizer.y_mean)/standardizer.y_std
        mask[i, :m] = True
    return {name: torch.as_tensor(value, device=device) for name, value in
            {"x": x, "time": t, "y": y, "mask": mask}.items()}


def batches(trajectories, size, shuffle=False, rng=None):
    indices = np.arange(len(trajectories))
    if shuffle:
        rng.shuffle(indices)
    for start in range(0, len(indices), size):
        yield [trajectories[i] for i in indices[start:start+size]]
