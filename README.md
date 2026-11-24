# NanoRemix

**从 Oxford Nanopore RNA002 / RNA004 直接 RNA 测序数据预测 RNA 修饰。**

[![软件测试](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/iris982022/nanoremix_tool/actions/workflows/tests.yml)

源码 **v0.3.1** · 预训练模型 **v0.3.0** · 建议使用 **Python 3.10 或 3.11**

[预处理](#3-如何准备输入) · [常用参数](#4-常用参数) · [快速上手](#5-快速上手) · [完整手册](docs/README.zh-CN.md) · [下载模型](https://github.com/iris982022/nanoremix_tool/releases/tag/v0.3.0)

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

Multi 在同一参考碱基家族内比较。三代模型使用序列和电流特征，也可直接读入预计算三代概率。二代模型使用同中心的 81 nt 序列，已有完整 15 类二代概率可直接提供。融合使用两个分支的证据，无需额外进行 Illumina 测序。[各路线的输入与支持类别](docs/README.zh-CN.md#choose)

## 3. 如何准备输入

从 FAST5 或 POD5 开始时，先用对应的上游工具完成 basecalling、比对和信号特征提取。已有规范特征表或概率表时，可直接进入预测。

| 原始数据 | 上游步骤 | NanoRemix 转换 |
| --- | --- | --- |
| RNA002 FAST5 | Guppy（`--fast5_out`）→ 拆成单 read FAST5 → Tombo/GraphMap → 提取 Step4 TSV 和坐标侧表 | `prepare --chemistry RNA002 --format step4`，生成 9 nt 上下文和 20 个信号值 |
| RNA004 POD5 | Dorado（`--emit-moves`）→ 比对并排序 BAM → Remora 提取信号特征 → 合并为 Step6 CSV | `prepare --chemistry RNA004 --format step6`，从正链 Step6 生成 17 nt 上下文和 36 个信号值 |

特征表使用带表头的 CSV/TSV，可 gzip 压缩。`9mer_final` 保存短序列上下文，四组信号列为 `mean_output1…K`、`median_output1…K`、`std_output1…K` 和 `length_output1…K`，RNA002 的 K=5，RNA004 的 K=9。按序列做融合时另需同中心、同方向的 81 nt 序列。[完整预处理命令](docs/README.zh-CN.md#preprocessing) · [输入字段与示例](docs/README.zh-CN.md#inputs)

> RNA004 的旧 Step6 转换要求显式 `strand=+`。负链数据需要先用经过验证的流程构造符合规范的 prepared 表，再运行预测。[RNA004 使用限制](docs/KNOWN_LIMITATIONS.zh-CN.md)

## 4. 常用参数

以下选项用于 `predict`，也适用于 `validate`。用 `--input` 指定输入表，预测时用 `--output` 指定新目录。

| 参数 | 怎么选 |
| --- | --- |
| `--chemistry RNA002` / `RNA004` | 选实验使用的试剂版本。 |
| `--task binary` / `multi` | 选修饰与未修饰配对，或同一参考碱基家族内的多类比较。 |
| `--mods m6A,m5C` / `all` | 选修饰列表或当前路线支持的全部修饰，默认加入同家族未修饰对照。Multi 可用 `--targets` 直接列出参与比较的类别。 |
| `--mode thirdgen` / `--data-type biological` | 分别运行纯三代或融合。`--mode` 与 `--data-type` 二选一。 |
| `--strategy s1` / `s2` / `both` | 默认 `s2`，`both` 分别输出两种策略的结果。 |
| `--input-kind features` / `probabilities` | 默认读特征表。已有三代概率时选 `probabilities`。 |
| `--sequence-column sequence81` | 指定用于二代推理的 81 nt 序列列。已有独立二代概率表时可用 `--second-probabilities` 和 `--join-on` 按观测键匹配。 |
| `--model-root models` | 指定模型根目录。加载已核实的权重时需同时使用 `--trust-model-files`，可用 `--model-checksums models/SHA256.json` 检查文件。 |
| `--min-confidence 0.8` / `--output-filter modified` | 分别设置获胜类别的条件概率下限、仅保存已接受的修饰预测。这里的 0.8 是示例，默认下限为 0，默认保存所有候选。 |

[全部参数、默认值和使用条件](docs/README.zh-CN.md#options)

## 5. 快速上手

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

## 6. 运行自己的数据

1. 按试剂和任务安装对应模型，用 `doctor` 检查文件和依赖。[安装与模型预检](docs/README.zh-CN.md#start-install)
2. 用 `validate` 检查输入和推理，再用 `predict` 保存结果。[四套配置的运行命令](docs/README.zh-CN.md#start-run)
3. 查看 `predictions.tsv`、`selection_audit.tsv` 和 `run_manifest.json`，核对候选状态、预测标签和运行参数后再汇总。[结果分析](docs/README.zh-CN.md#outputs)
