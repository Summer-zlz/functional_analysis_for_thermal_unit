# FHuber 复现

## 文件

- `run_fhuber.R`：入口，完成三折时序交叉验证、独立测试及全样本拟合。
- `fhuber_time_cv.R`：FHuber 核心训练与预测函数。
- `coaldata_*_unit1.csv`：三个数据集副本。
- `results/`：五个核心 CSV，分别为交叉验证结果、最优参数、训练/测试指标、全样本指标及选中变量。

## 运行

需要 R，以及 `fda`、`Matrix` 包，无需 Python。在本文件夹打开 PowerShell，执行：

```powershell
& 'C:\Program Files\R\R-4.1.3\bin\Rscript.exe' run_fhuber.R
```

参数：`psi=10^seq(-2,10,length.out=15)`、`alpha=0.0001`、`lambda=c(0.01,0.1,0.5,1)`、`max_iter=500000`，保留原提前停止条件。

按完整 ID 进行三折滚动验证，按平均验证 RMSE 选参，在前 80% 重训、后 20% 测试。全样本拟合沿用所选参数，UP/DOWN 各训练两个原有簇，STABLE 使用全部样本。

精简版直接使用原 R 求解器，完整计算耗时较长，浮点运算可能产生微小差异。运行会覆盖汇总表，并新增逐折记录、预测和模型文件；如需保留当前结果，请先备份 `results/`。MAPE 单位为 %，达到迭代上限不代表严格收敛。
