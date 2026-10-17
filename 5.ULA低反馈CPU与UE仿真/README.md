# ULA 低反馈 CPU/UE 仿真

本文件夹包含 ULA 64 阵元下 CPU 端与 UE 端共相位的正式仿真代码、统一输入和 400 位置结果，按本仓库既有编号目录方式新增。

详见 [ULA64_仿真结果说明.md](ULA64_仿真结果说明.md)。正式结果在 `results/full400_ULA_CPU_20261017/` 和 `results/full400_ULA_UE_20261017/`，跨端核验与 B=8、30 dB 的对照表在 `results/` 下。

MATLAB R2024a 复现：

```matlab
verify_ula;
run_ula_cpu('Mode','full');
run_ula_ue('Mode','full');
audit_ula_results;
```

仿真依赖 `shared/` 下固定的几何、位置、散射和噪声输入。正式重跑前请指定新的 `OutputDir`，以保留本次结果。
