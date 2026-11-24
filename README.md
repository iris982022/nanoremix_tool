# NanoRemix

**从 Oxford Nanopore RNA002 / RNA004 直接 RNA 测序数据预测 RNA 修饰。**

[![软件测试](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml)

源码 **v0.3.1** · 预训练模型 **v0.3.0** · 建议使用 **Python 3.10 或 3.11**

[完整中文 README](docs/README.zh-CN.md) · [下载模型](https://github.com/iris982022/nanoremix_tool/releases/tag/v0.3.0) · [输入格式](docs/INPUTS.md) · [使用限制](docs/KNOWN_LIMITATIONS.zh-CN.md)

## 能做什么

按实验使用的 ONT 试剂选择 RNA002 或 RNA004，再按任务选择 Binary 或 Multi。特征表每行对应一条 read 上的一个候选位置。

| 配置 | 支持的修饰 | 每条候选的三代特征输入 |
| --- | --- | --- |
| [RNA002 Binary](docs/README.zh-CN.md#rna002-binary) | m6A、m5C、假尿苷（`pseudoU`），各与对应未修饰碱基比较 | 9 nt 上下文 + 20 个信号特征 |
| [RNA002 Multi](docs/README.zh-CN.md#rna002-multi) | 14 种修饰，在同一参考碱基家族内比较 | 9 nt 上下文 + 20 个信号特征 |
| [RNA004 Binary](docs/README.zh-CN.md#rna004-binary) | m6A、m5C、假尿苷（`pseudoU`），各与对应未修饰碱基比较 | 17 nt 上下文 + 36 个信号特征 |
| [RNA004 Multi](docs/README.zh-CN.md#rna004-multi) | 12 种修饰，在同一参考碱基家族内比较 | 17 nt 上下文 + 36 个信号特征 |

- **纯三代**使用序列和电流特征，或读取符合格式的预计算三代概率。特征推理只需下载对应配置的模型包。
- **二代＋三代融合**额外使用同一候选中心的 81 nt 序列及共享 SecondGen 模型，也可直接提供完整的 15 类二代概率。两代指 NanoRemix 的模型分支，融合无需额外进行 Illumina 测序。
- **策略与输出**可指定目标修饰和 S1/S2 策略，默认用 S2，也可分别输出两种策略的结果。每条路线输出逐 read 候选的预测标签、类别概率和筛选记录，运行清单记录输入、参数和模型来源。

> 手头只有 FAST5 / POD5 时，请先完成 basecalling、比对和特征提取。内置 RNA004 旧 Step6 转换仅接受显式正链表，RNA004 正负链均可在正确准备特征后用于预测。[上游流程](docs/README.zh-CN.md#preprocessing) · [RNA004 使用限制](docs/KNOWN_LIMITATIONS.zh-CN.md)

## 快速上手

用两条人工概率记录检查安装。这个示例无需下载模型权重。

```bash
git clone https://github.com/iris982022/nanoremix_tool.git
cd nanoremix_tool
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .

nanoremix predict --input examples/minimal_binary_probabilities.csv \
  --chemistry RNA004 --task binary --mods m6A \
  --input-kind probabilities --mode thirdgen --strategy s2 \
  --output results/quickstart
```

在 `results/quickstart/thirdgen/s2/predictions.tsv` 中可看到两条预测，依次为 `m6A` 和 `A`。每次运行请选择新的输出目录。

## 运行自己的数据

1. 按上表选择配置并准备输入表。[输入字段与示例](docs/README.zh-CN.md#inputs)
2. 按需安装推理依赖，在仓库根目录解压对应模型包，得到 `models/` 目录。[安装与模型预检](docs/README.zh-CN.md#start-install)
3. 参考选定配置的命令运行 `nanoremix predict`，查看预测表、筛选记录和运行清单。[四套配置的运行命令](docs/README.zh-CN.md#start-run) · [输出说明](docs/README.zh-CN.md#outputs)

用 `nanoremix catalog` 查询各配置的支持类别，用 `nanoremix doctor` 检查模型文件和依赖。[全部参数](docs/README.zh-CN.md#options) · [融合计算](docs/README.zh-CN.md#calculation)
