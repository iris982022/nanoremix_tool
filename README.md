# NanoRemix

**从 Oxford Nanopore RNA002 / RNA004 直接 RNA 测序数据预测 RNA 修饰。**

[![软件测试](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml)

源码 **v0.3.1** · 预训练模型 **v0.3.0** · 建议使用 **Python 3.10 或 3.11**

[完整中文 README](docs/README.zh-CN.md) · [下载模型](https://github.com/iris982022/nanoremix_tool/releases/tag/v0.3.0) · [输入格式](docs/INPUTS.md) · [使用限制](docs/KNOWN_LIMITATIONS.zh-CN.md)

## 1. NanoRemix 能做什么

- 用 `prepare` 把 RNA002 的 Step4/extracted 特征表或 RNA004 的正链 Step6 表转换成预测输入。
- 用 `catalog` 查询支持类别，用 `doctor` 检查模型与依赖，用 `validate` 试跑并检查输入。
- 用 `predict` 运行纯三代、二代＋三代融合或纯二代预测，输出逐 read 候选的标签、概率、筛选记录和运行清单。

## 2. 支持哪些试剂和修饰

按实验使用的 ONT 试剂选择 RNA002 或 RNA004，再按任务选择 Binary 或 Multi。特征表每行对应一条 read 上的一个候选位置。

| 配置 | 支持的修饰 | 三代特征输入 |
| --- | --- | --- |
| [RNA002 Binary](docs/README.zh-CN.md#rna002-binary) | m6A、m5C、假尿苷（`pseudoU`），各与对应未修饰碱基比较 | 9 nt 上下文 + 20 个信号特征 |
| [RNA002 Multi](docs/README.zh-CN.md#rna002-multi) | 14 种。m6A、m1A、Inosine、Am、m5C、hm5C、ac4C、FormylC、m7G、o8G、pseudoU、m1psi、m5U、Um | 9 nt 上下文 + 20 个信号特征 |
| [RNA004 Binary](docs/README.zh-CN.md#rna004-binary) | m6A、m5C、假尿苷（`pseudoU`），各与对应未修饰碱基比较 | 17 nt 上下文 + 36 个信号特征 |
| [RNA004 Multi](docs/README.zh-CN.md#rna004-multi) | 12 种。m6A、m1A、Am、m5C、hm5C、ac4C、o8G、Gm、pseudoU、m1psi、m5U、Um | 17 nt 上下文 + 36 个信号特征 |

Multi 在同一参考碱基家族内比较。三代使用对应试剂的序列、电流特征或预计算概率，支持 S1/S2，默认 S2。二代使用同中心的 81 nt 序列与共享 SecondGen 模型，已有完整 15 类二代概率可直接提供。融合使用两个分支的证据，无需额外进行 Illumina 测序。[各路线的输入与支持类别](docs/README.zh-CN.md#choose)

> 手头只有 FAST5 / POD5 时，请先完成 basecalling、比对和特征提取。内置 RNA004 旧 Step6 转换仅接受显式正链表，RNA004 正负链均可在正确准备特征后用于预测。[上游流程](docs/README.zh-CN.md#preprocessing) · [RNA004 使用限制](docs/KNOWN_LIMITATIONS.zh-CN.md)

## 3. 快速上手

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

## 4. 运行自己的数据

1. 将原始数据转换为特征表，或从已有特征表和概率表开始。[原始数据预处理](docs/README.zh-CN.md#preprocessing)
2. 核对列名、候选身份与融合所需的 81 nt 上下文或二代概率。[文件格式与示例](docs/README.zh-CN.md#inputs)
3. 安装对应模型，先检查，再按所选配置预测。[安装与模型预检](docs/README.zh-CN.md#start-install) · [四套配置的运行命令](docs/README.zh-CN.md#start-run)
4. 核查审计状态、预测标签和证据来源，再进行比较或汇总。[结果分析](docs/README.zh-CN.md#outputs)

[所有命令的参数、默认值与允许值](docs/README.zh-CN.md#options) · [融合计算](docs/README.zh-CN.md#calculation)
