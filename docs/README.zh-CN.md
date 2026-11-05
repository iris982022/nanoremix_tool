<div align="center">

# 🧬 NanoRemix

**支持 RNA002 / RNA004 的 RNA 修饰预测工具**

Binary 与 Multi · 指定修饰筛选 · 纯三代与二代三代融合

[English](../README.md) · [中文功能总览](OVERVIEW.zh-CN.md) · [输入格式](INPUTS.md) · [模型文件](MODELS.md) · [完整选项](OPTIONS.md)

**v0.3.0** · 配套预训练模型 · Python 3.10/3.11

</div>

> **下载选择：** 源码仓库不保存预训练权重。从 GitHub Release 按需下载
> chemistry/task 对应的三代模型 ZIP；二代序列融合另加共享二代模型 ZIP。
> 每个模型 ZIP 更小，仍包含原始推理需要的权重与校验索引，无须重新训练。
> 安装方式见 [MODELS.md](MODELS.md)。

## 🧩 功能与输出规则

输入预处理后的序列信号特征，或已有模型的概率表，指定修饰类型后输出预测。
未指定的修饰不参与最终竞争；对应未修饰碱基保留为对照。

| 数据类型 | 默认结果 | 可选额外结果 |
| --- | --- | --- |
| 生物数据 `biological` | 二代与三代融合 `fused` | `--also-thirdgen` 纯三代；`--also-secondgen` 纯二代 |
| 合成数据 `synthetic` | 纯三代 `thirdgen` | 不需要、不读取二代输入 |

这里“二代/三代”是 **NanoRemix 模型版本**，不是测序技术代际。
二代提供序列证据，三代提供与化学版本匹配的序列信号证据。
生物数据缺少二代输入时会报错，不会自动冒充融合结果。
融合结果中，二代不支持的类别仍保留三代证据，并在结果中标注来源。

| 常用选项 | 含义 |
| --- | --- |
| `--chemistry RNA002 / RNA004` | 指定化学版本 |
| `--task binary / multi` | 修饰与未修饰配对比较 / 参考碱基家族内多类别比较 |
| `--mods m6A,m5C,pseudoU` | 只考虑指定修饰及其未修饰对照；`all` 为全部支持修饰 |
| `--strategy s1 / s2 / both` | 独立位移视图 / 联合多位置模型 / 分别输出两者；默认 `s2` |
| `--output-filter all / retained / modified` | 所有观测 / 接受的观测 / 接受的修饰观测 |

S1 内部多个位置模型的组合不等于二代与三代融合。
`both` 分别输出 S1/S2，不生成二者的共识结果。

## 🚀 安装与无权重示例

在仓库根目录运行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .

nanoremix predict \
  --input examples/minimal_binary_probabilities.csv \
  --task binary --chemistry RNA004 --mods m6A \
  --input-kind probabilities --data-type synthetic --strategy s2 \
  --output results/quickstart
```

Windows PowerShell 激活命令为 `.venv\Scripts\Activate.ps1`。
示例是两行人工概率，输出分别为 m6A 和 A；不能作为准确率证据。
每次运行使用新的输出目录。

## 📥 输入

支持 UTF-8 CSV/TSV，可 gzip 压缩。一行对应一次候选观测，例如某条 read 覆盖一个位点。

| 输入路线 | 必需内容 |
| --- | --- |
| 特征输入，默认 | `9mer_final`、信号描述量、对应三代模型 |
| Binary 概率输入 | `center_ref_base`；各请求策略/修饰的 `s1_prob_<mod>` 或 `s2_prob_<mod>` |
| Multi 概率输入 | `center_ref_base`；对应策略的完整原生类别概率，不能只输入筛选后的子集 |
| 生物数据的二代来源 | 完整 15 类二代概率；或中心对齐的 81 nt 序列和二代模型 |

特征列 `9mer_final` 对 RNA002 是 **9 nt**，对 RNA004 是 **17 nt**。
信号为 `mean_output1..K`、`median_output1..K`、`std_output1..K`、
`length_output1..K`，K 分别为 5、9。必须采用原模型一致的预处理和归一化。
本包不直接处理原始 POD5/FAST5，不负责 basecalling 或事件比对。

`record_key` 标识单次观测；`sample_id/read_id/site_id` 和参考坐标可一并保留。
独立二代概率表按明确身份列连接，不能按行顺序拼接。
81 nt 上下文不能从 9/17 nt 短上下文凭空生成。
详细列名、FASTA 和坐标规则见 [INPUTS.md](INPUTS.md)。

## ⚙️ 使用真实模型

建议在 Python 3.10/3.11 环境中安装模型依赖，把已确认模型放入 `models/` 后检查：

```bash
python -m pip install -e '.[models,second]' -c constraints-models.txt
nanoremix doctor --model-root models --model-checksums models/SHA256.json --strict
```

从 GitHub Release 下载所需模型压缩包，在仓库根目录解压，逐步生成/补充
`models/`。纯三代需一个 chemistry/task 模型包；用81 nt现场运行二代的融合
还需共享二代模型包。已有完整二代概率时不需下载二代权重。
检查会列出缺失文件、无效 manifest、空文件或 Git LFS 指针等问题。
它不反序列化模型；通过静态检查后，仍须用真实模型与已知输出验证预测。

生物数据，输入特征表内含 `sequence81` 列：

```bash
nanoremix predict \
  --input prepared004.csv --task multi --chemistry RNA004 \
  --mods m6A,m1A,m5C,pseudoU --data-type biological --strategy s2 \
  --sequence-column sequence81 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files \
  --output results/biological
```

加 `--also-thirdgen --also-secondgen` 可保存两种纯代际结果。
若已有完整二代概率，用
`--second-probabilities second.tsv --join-on record_key`
替代 `--sequence-column`，不必再次加载二代网络。

合成数据，只用三代：

```bash
nanoremix predict \
  --input prepared002.csv --task binary --chemistry RNA002 \
  --mods m6A,pseudoU --data-type synthetic --strategy s2 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files \
  --output results/synthetic
```

上述预测命令需要自行提供预处理数据。可用
`python examples/make_feature_demo.py --output demo-features`
生成 `RNA002_features.csv` / `RNA004_features.csv`，检查模型执行。
这些人工特征只用于运行检查，不是生物学准确率证据。只加载可信来源模型。

## 📤 输出

| 路径 | 内容 |
| --- | --- |
| `fused/s2/predictions.tsv` | 默认生物数据融合结果 |
| `thirdgen/s2/predictions.tsv` | 合成数据纯三代结果，或生物数据额外三代结果 |
| `secondgen/none/predictions.tsv` | 可选纯二代结果，只保存一次 |
| 各结果旁的 `selection_audit.tsv` | 全部观测的保留/拒绝状态及原因 |
| `run_manifest.json` | 输入、参数、类别范围、模型哈希和处理计数 |

S1 使用相应的 `s1/` 目录。预测表包含观测身份、`predicted_label`、
选定类别的 `prob_<label>`、`confidence`、`status`、`evidence_source`。
拒绝的行不输出已接受标签。

概率以选定类别和参考碱基家族为条件；`confidence` 是获胜类别概率，
不总是修饰概率。未纳入的真实修饰可能被分到已选类别，不能据此判断其不存在。
结果不会自动汇总为位点修饰率。纯二代输出不能预测二代不支持的修饰。

## 📚 进一步使用

`nanoremix catalog` 查询支持类别，`nanoremix predict --help` 查看全部参数。
模型、校准参数和目录格式见 [MODELS.md](MODELS.md)，高级选项见 [OPTIONS.md](OPTIONS.md)。
该仓库面向推理使用，不等于自动复现论文所有图表。
公开分发前须由权利人确认 [许可证](../LICENSE)。
