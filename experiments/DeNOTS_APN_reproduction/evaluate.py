"""Point-weighted metrics in original response units."""
import numpy as np


def metrics(actual, predicted):
    actual, predicted = np.asarray(actual, dtype=np.float64), np.asarray(predicted, dtype=np.float64)
    assert actual.shape == predicted.shape and np.isfinite(actual).all() and np.isfinite(predicted).all()
    error = predicted - actual
    return {"MAE": float(np.abs(error).mean()), "RMSE": float(np.sqrt(np.square(error).mean())),
            "MAPE": float(100 * np.abs(error / actual).mean()) if np.all(actual != 0) else float("nan")}
