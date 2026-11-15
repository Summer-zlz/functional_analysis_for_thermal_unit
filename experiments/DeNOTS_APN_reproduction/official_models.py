"""Load pinned author implementations and add trajectory-to-response interfaces."""
import __future__
import difflib
from pathlib import Path
import sys
import types
from types import SimpleNamespace
import bootstrap
import torch
from torch import nn

HERE = bootstrap.HERE
VENDOR = bootstrap.VENDOR


def package(name, path):
    if name not in sys.modules:
        module = types.ModuleType(name)
        module.__path__ = [str(path)]
        sys.modules[name] = module


def load_source(name, path, replacements=()):
    if name in sys.modules:
        return sys.modules[name]
    source = Path(path).read_text(encoding="utf-8")
    updated = source
    for original, replacement in replacements:
        assert original in updated, f"Unexpected pinned source in {path}"
        updated = updated.replace(original, replacement)
    # Future annotations and mechanical match->if conversions support the installed
    # Python 3.9; source files remain untouched and numerical operations unchanged.
    module = types.ModuleType(name)
    module.__file__ = str(path)
    module.__package__ = name if Path(path).name == "__init__.py" else name.rpartition(".")[0]
    if Path(path).name == "__init__.py":
        module.__path__ = [str(Path(path).parent)]
    sys.modules[name] = module
    code = compile(updated, str(path), "exec", flags=__future__.annotations.compiler_flag)
    exec(code, module.__dict__)
    if updated != source:
        diff = "".join(difflib.unified_diff(source.splitlines(True), updated.splitlines(True),
                                          fromfile=str(path.relative_to(VENDOR)), tofile="in-memory compatibility"))
        out = HERE / "results" / "compatibility"
        out.mkdir(parents=True, exist_ok=True)
        (out / (name.replace(".", "_")+".patch")).write_text(diff, encoding="utf-8")
    return module


def denots_classes():
    root = VENDOR / "denots" / "src"
    for name, directory in [("src", root), ("src.nn", root/"nn"),
                            ("src.nn.forecasting", root/"nn"/"forecasting"),
                            ("src.interp", root/"interp")]:
        package(name, directory)
    load_source("src.mask_utils", root/"mask_utils.py")
    load_source("src.interp.base", root/"interp"/"base.py")
    load_source("src.interp.natural_cubic", root/"interp"/"natural_cubic.py", [
        ('''        match t.ndim:
            case 3:
                fit_vmapped = torch.vmap(fit_vmapped, -1, -1)
                eval_vmapped = torch.vmap(eval_vmapped, (-1, -1, None), -1)
            case 2:
                fit_vmapped = torch.vmap(fit_vmapped, (None, -1), -1)
            case _:
                raise ValueError(t)''',
         '''        if t.ndim == 3:
            fit_vmapped = torch.vmap(fit_vmapped, -1, -1)
            eval_vmapped = torch.vmap(eval_vmapped, (-1, -1, None), -1)
        elif t.ndim == 2:
            fit_vmapped = torch.vmap(fit_vmapped, (None, -1), -1)
        else:
            raise ValueError(t)''')])
    load_source("src.interp.forecasting", root/"interp"/"forecasting.py", [
        ('''        match self.forecast_method:
            case "zero":
                return torch.where(out_of_bounds, 0, x)
            case "last":
                return x
            case _:
                raise ValueError()''',
         '''        if self.forecast_method == "zero":
            return torch.where(out_of_bounds, 0, x)
        elif self.forecast_method == "last":
            return x
        else:
            raise ValueError()''')])
    load_source("src.nn.vfs.interprnn", root/"nn"/"vfs"/"interprnn.py")
    load_source("src.nn.vfs", root/"nn"/"vfs"/"__init__.py")
    load_source("src.nn.forecasting.base", root/"nn"/"forecasting"/"base.py")
    model = load_source("src.nn.forecasting.denots", root/"nn"/"forecasting"/"denots.py")
    return model.DeNOTS, sys.modules["src.interp.natural_cubic"].NaturalCubicSpline


class DeNOTSRegression(nn.Module):
    def __init__(self, hidden_dim, input_dim=17, depth=10):
        super().__init__()
        backbone, interpolator = denots_classes()
        self.backbone = backbone(input_dim=input_dim, hidden_dim=hidden_dim,
            interp=interpolator(), vf_type="adaptive", forecast_method="last",
            nf=True, depth=depth, tol=1e-3)
        self.readout = nn.Linear(hidden_dim, 1)

    def forward(self, x, time, mask):
        # The author spline uses +inf for absent knots. Query times repeat the
        # endpoint in storage only; every padded prediction is excluded by mask.
        knots = time.masked_fill(~mask, torch.inf)
        import torchode as to
        knots = knots*self.backbone.depth
        self.backbone.interp.fit(knots, x)
        self.backbone.interp.tmax = knots.masked_fill(~mask, -torch.inf).amax(dim=1)
        problem = to.InitialValueProblem(y0=self.backbone.start_proj(x[:, 0]),
            t_start=knots[:, 0], t_eval=time*self.backbone.depth)
        solution = self.backbone.solver.solve(problem, term=self.backbone.term)
        if (solution.status != 0).any():
            raise RuntimeError(f"DeNOTS ODE failed: {solution.status.tolist()}")
        return self.readout(solution.ys).squeeze(-1)


class APNRegression(nn.Module):
    def __init__(self, hidden_dim, input_dim=17, patches=8, time_dimension=16):
        super().__init__()
        module = load_source("official_apn", VENDOR/"apn"/"models"/"APN.py", [
            ("from utils.globals import logger\n", ""),
            ("from utils.ExpConfigs import ExpConfigs\n", "")])
        configs = SimpleNamespace(d_model=hidden_dim, apn_te_dim=time_dimension, enc_in=input_dim,
            apn_npatch=patches, apn_nlayer=1, apn_attn_heads=1, dropout=.1)
        self.backbone = module.IMTS_SubModel(configs)
        # Retain TAPA, query aggregation and author per-channel/time decoder.
        # Combine covariate-channel outputs into one response; no coal-consumption channel.
        self.readout = nn.Linear(input_dim, 1)

    def forward(self, x, time, mask):
        marks = time.unsqueeze(-1)
        channel_mask = mask.unsqueeze(-1).expand_as(x).float()
        features = self.backbone(x, marks, channel_mask, marks)
        return self.readout(features).squeeze(-1)


def make_model(name, hidden_dim, input_dim=17, denots_depth=10,
               apn_patches=8, apn_time_dimension=16):
    if name == "DeNOTS":
        return DeNOTSRegression(hidden_dim, input_dim, denots_depth)
    if name == "APN":
        return APNRegression(hidden_dim, input_dim, apn_patches, apn_time_dimension)
    raise ValueError(name)
