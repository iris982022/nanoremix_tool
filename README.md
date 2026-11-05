# NanoRemix

**RNA002 / RNA004 RNA修饰预测工具：4套核心配置，各支持纯三代和二代＋三代融合。**

工具版本 **v0.3.1** · 建议Python **3.10/3.11** · 使用已发布的 **v0.3.0模型权重**

**已有正确prepared特征表并安装匹配模型，即可预测。** 只有原始FAST5/POD5时，先按下文完成basecalling、比对、信号特征和上下文准备。RNA004旧Step6转换入口仅接受显式正链表；正确准备的正/负链prepared表均可交给predict，但当前没有经过验证的任意未知样本通用负链转换器。

[四种方式](#choose) · [安装和模型下载](#start-install) · [① RNA002 Binary](#rna002-binary) · [② RNA002 Multi](#rna002-multi) · [③ RNA004 Binary](#rna004-binary) · [④ RNA004 Multi](#rna004-multi) · [参数与筛选](#options) · [计算方法](#calculation) · [上游详细流程](#preprocessing)

<a id="choose"></a>
<a id="modifications"></a>

## 四种方式一览：应该选哪一种

| 方式 | 数据对应chemistry | 做什么 | 正式支持的修饰数 | 纯三代输入 |
|---|---|---|---:|---|
| [① RNA002 Binary](#rna002-binary) | SQK-RNA002 | 正常碱基与一个专用修饰二分类 | 3：m6A、m5C、pseudoU | 9 nt上下文＋20个信号值 |
| [② RNA002 Multi](#rna002-multi) | SQK-RNA002 | 在同一正常碱基家族内比较多种修饰 | 14，另有4个正常类 | 9 nt上下文＋20个信号值 |
| [③ RNA004 Binary](#rna004-binary) | SQK-RNA004 | 正常碱基与一个专用修饰二分类 | 3：m6A、m5C、pseudoU | 17 nt上下文＋36个信号值 |
| [④ RNA004 Multi](#rna004-multi) | SQK-RNA004 | 在同一正常碱基家族内比较多种修饰 | 12，另有4个正常类 | 17 nt上下文＋36个信号值 |

按实验使用的ONT试剂/化学版本选择RNA002或RNA004，不能凭文件名决定或互换权重。只研究m6A/m5C/pseudoU的专用pair时选Binary；需要m1A、Am等多类比较时选对应Multi。

**四种方式都可以选：**

- **纯三代**：`--mode thirdgen`，读取对应chemistry的序列＋信号特征，只需对应三代模型包。
- **二代＋三代融合**：`--data-type biological`，在三代输入之外增加同中心的连续81 nt序列和SecondGen模型；也可用完整二代概率代替序列推理。
- 每条路线可选S1/S2，默认 `--strategy s2`；`--strategy both` 分别输出，不计算共识。4套配置×2条主要路线＝8条路线，进一步×2个策略＝16种route/strategy组合。

这里“二代/三代”是**NanoRemix内部模型分支**：二代用连续81 nt序列，三代用序列＋Nanopore电流特征；**无需额外Illumina测序**。Binary/Multi、RNA002/RNA004、S1/S2和纯三代/融合分别是不同选择。

<a id="start-install"></a>

## 安装与下载所需模型

获取源码，在根目录安装；建议新建Python3.10/3.11环境：

```bash
git clone https://github.com/iris982022/nanoremix_tool.git
cd nanoremix_tool
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[models,second]' -c constraints-models.txt
nanoremix --version
```

也可 **Code → Download ZIP** 后解压进入源码根目录，从创建环境开始。Windows PowerShell激活命令为 `.venv\Scripts\Activate.ps1`。

从 [v0.3.0模型Release](https://github.com/iris982022/nanoremix_tool/releases/tag/v0.3.0) 选择附件，**下载后解压到项目根目录，形成 `models/`**：

| 配置 | 对应三代模型下载 |
|---|---|
| RNA002 Binary | [RNA002_Binary ZIP](https://github.com/iris982022/nanoremix_tool/releases/download/v0.3.0/NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip) |
| RNA002 Multi | [RNA002_Multi ZIP](https://github.com/iris982022/nanoremix_tool/releases/download/v0.3.0/NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip) |
| RNA004 Binary | [RNA004_Binary ZIP](https://github.com/iris982022/nanoremix_tool/releases/download/v0.3.0/NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip) |
| RNA004 Multi | [RNA004_Multi ZIP](https://github.com/iris982022/nanoremix_tool/releases/download/v0.3.0/NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip) |
| 现场计算二代的融合 | 在对应三代包之外，再下载 [SecondGen ZIP](https://github.com/iris982022/nanoremix_tool/releases/download/v0.3.0/NanoRemix_v0.3.0_models_SecondGen_20261104.zip) |

纯三代只取对应一包；已有完整二代概率时无需SecondGen。比如RNA004 Multi融合，选中两包后：

```bash
unzip -o NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip -d .
unzip -o NanoRemix_v0.3.0_models_SecondGen_20261104.zip -d .

nanoremix doctor --model-root models --model-checksums models/SHA256.json \
  --chemistry RNA004 --task multi --mods all --strategy s2 \
  --data-type biological --strict
```

doctor只检查资产和依赖。纯三代资产预检用 `--data-type synthetic` 选择其检查范围；实际生物样本纯三代预测仍用 `--mode thirdgen`。

<a id="start-run"></a>
<a id="start-inputs"></a>

## 四种方式逐一使用

以下均在源码根目录、已激活的NanoRemix环境中运行，先下载并解压对应模型。`prepared002.csv`、`prepared004.csv`是你自己的真实特征表；`*_with81.csv`是在同一表增加有效 `sequence81` 后的文件。**每次选择尚不存在的输出目录。** 上游版本来自现有日志、配置及环境元数据，不能视为每个历史作业的完整锁定记录；basecalling、Tombo/Remora和现代推理分别使用独立环境。

主表为UTF-8 CSV/TSV，可gzip压缩；**一行＝一个read在一个候选中心的观测**。建议提供非空、全表唯一 `record_key`，并保留 `sample_id/read_id/site_id`、参考、坐标和方向。以下人工CSV只演示字段，不是可以评价准确率的实验数据。

<a id="rna002-binary"></a>

## 方式一：RNA002 Binary

### 用途与修饰

用于RNA002数据的专用二分类：**A/m6A、C/m5C、U/pseudoU**。可一次指定三种，实际按中心参考家族运行相应pair；不支持G家族修饰，也不把m1A等加入Binary。`--mods m6A`只比较A/m6A；`--mods all`选这三种及各自正常对照。

### 输入是什么样

| 内容 | RNA002 Binary要求 |
|---|---|
| `9mer_final` | 9 nt，例如 `A_C_G_U_A_C_G_U_A`，候选中心是第5位A |
| 数值列 | `mean_output1..5`、`median_output1..5`、`std_output1..5`、`length_output1..5`，共20个值 |
| 融合额外输入 | 同方向、同中心的81 nt列 `sequence81`，或完整15类二代概率 |
| 格式文件 | [prepared＋81 nt人工CSV](docs/guide_20261104/examples/RNA002_prepared_with81.csv)、[上游Step4人工TSV](docs/guide_20261104/examples/RNA002_step4.tsv) |

### 输入怎样得到：上游步骤与版本

1. 原始multi-read FAST5 → **Guppy 6.3.8+d9e0f64**，RNA002配置 `rna_r9.4.1_70bps_hac.cfg`（FLO-MIN106 / SQK-RNA002），保留basecalled FAST5（`--fast5_out`）。
2. **ont-fast5-api 4.1.1** 的 `multi_to_single_fast5` → single-read FAST5。
3. **Tombo 1.0＋GraphMap 0.5.2**，对匹配的参考cDNA执行resquiggle；保留corrected Events/Alignment，并检查失败记录。
4. 包内 [extract_tombo_step4.py](upstream_reference/RNA002/python/extract_tombo_step4.py)，原特征环境Python3.6.9、NumPy1.19.2、h5py2.10.0，使用 `--rna --clip 10` → Step4 `.feature.tsv`＋同步 `.coordinates.tsv`。Step4含 `id/base/sequence_input/mean/std/median/length`：序列是5 nt，四组统计各是5个竖线分隔值。
5. NanoRemix将5行连续、正确重叠的5-mer构造为9 nt，**保留中间行20个信号值**：

```bash
nanoremix prepare --input sample002.step4.feature.tsv \
  --output prepared002.csv --chemistry RNA002 --format step4
```

纯三代到此即可。融合需再按sidecar身份和正确参考坐标取得连续81 nt，写入 `sequence81`；event index不能当作参考坐标，不能将短k-mer拼凑成81 nt。见[RNA002上游完整命令与连接示例](docs/guide_20261104/PIPELINE_命令与输入衔接.md)。

### 后续用NanoRemix运行

解压 **RNA002_Binary** 模型包；融合再解压 **SecondGen** 包。以下分别运行纯三代，以及融合并同时保存两代分支：

```bash
# 纯三代
nanoremix predict --input prepared002.csv --output results/002_binary_third \
  --chemistry RNA002 --task binary --mods m6A,m5C,pseudoU \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files

# 二代＋三代融合；使用额外具有正确sequence81的表
nanoremix predict --input prepared002_with81.csv --output results/002_binary_fused \
  --chemistry RNA002 --task binary --mods m6A,m5C,pseudoU \
  --data-type biological --sequence-column sequence81 --strategy s2 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

<a id="rna002-multi"></a>

## 方式二：RNA002 Multi

### 用途与修饰

用于RNA002的多类预测：模型原生**18类＝14种修饰＋A/C/G/U**；最终在候选中心对应的正常碱基家族内，比较你指定的类别。

| 中心参考家族 | 支持的修饰（另含对应正常对照） |
|---|---|
| A | m6A、m1A、Inosine、Am |
| C | m5C、hm5C、ac4C、FormylC |
| G | m7G、o8G |
| U | pseudoU、m1psi、m5U、Um |

`--mods all`使用全部上述修饰；例如 `--mods m6A,m1A`比较A/m6A/m1A。此配置的三代模型不支持Gm或Cm。融合默认extended仍保留14种修饰，但**hm5C、FormylC、o8G、m1psi只有三代证据**，不会因此获得二代支持。

### 输入是什么样

| 内容 | RNA002 Multi要求 |
|---|---|
| `9mer_final` | 9 nt，例如 `A_C_G_U_A_C_G_U_A`，中心第5位；使用正常参考序列 |
| 数值列 | mean/median/std/length各5列，列名同RNA002 Binary，共20个值 |
| 融合额外输入 | 同候选的连续81 nt `sequence81`，或完整15类二代概率 |
| 格式文件 | [prepared＋81 nt人工CSV](docs/guide_20261104/examples/RNA002_prepared_with81.csv)；与RNA002 Binary共享列格式，使用不同模型 |

### 输入怎样得到：上游步骤与版本

流程同RNA002 Binary，**不需要先运行Binary预测**：

1. FAST5 → **Guppy 6.3.8+d9e0f64**（`rna_r9.4.1_70bps_hac.cfg`，保留FAST5）。
2. **ont-fast5-api 4.1.1** 拆single-read FAST5 → **Tombo 1.0＋GraphMap 0.5.2** 做resquiggle。
3. 包内 `extract_tombo_step4.py --rna --clip 10`（Python3.6.9、NumPy1.19.2、h5py2.10.0）→ 5-mer Step4特征与坐标sidecar。
4. 以下转换得到9 nt＋20值prepared。不要提前按待测修饰删邻接事件，也不要对多个read先求均值。

```bash
nanoremix prepare --input sample002.step4.feature.tsv \
  --output prepared002.csv --chemistry RNA002 --format step4
```

融合时按sidecar连接真实参考坐标，再提取同中心81 nt；具体工具命令和坐标连接见[RNA002 pipeline](docs/guide_20261104/PIPELINE_命令与输入衔接.md)。

### 后续用NanoRemix运行

解压 **RNA002_Multi** 模型包；融合再解压 **SecondGen** 包。

```bash
# 纯三代：14种修饰；mods也可改成受支持子集
nanoremix predict --input prepared002.csv --output results/002_multi_third \
  --chemistry RNA002 --task multi --mods all \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files

# 二代＋三代融合，并额外保存纯三代/二代支持类别的纯二代结果
nanoremix predict --input prepared002_with81.csv --output results/002_multi_fused \
  --chemistry RNA002 --task multi --mods all \
  --data-type biological --sequence-column sequence81 --strategy s2 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

<a id="rna004-binary"></a>

## 方式三：RNA004 Binary

### 用途与修饰

用于RNA004数据的专用二分类：**A/m6A、C/m5C、U/pseudoU**，与RNA002 Binary支持同一组三种修饰，但特征长度、信号定义和权重都按RNA004，不能混用。可用 `--mods m6A`、`--mods m6A,m5C,pseudoU` 或 `--mods all`。

### 输入是什么样

| 内容 | RNA004 Binary要求 |
|---|---|
| `9mer_final` | **17 nt**，例如 `A_C_G_U_A_C_G_U_A_C_G_U_A_C_G_U_A`，中心第9位A；列名沿用历史名称 |
| 数值列 | `mean_output1..9`、`median_output1..9`、`std_output1..9`、`length_output1..9`，共36个值 |
| 融合额外输入 | 同方向、同中心的连续81 nt `sequence81`，或完整15类二代概率 |
| 格式文件 | [prepared＋81 nt人工CSV](docs/guide_20261104/examples/RNA004_prepared_with81.csv)、[正链Step6人工CSV](docs/guide_20261104/examples/RNA004_step6_positive.csv) |

### 输入怎样得到：上游步骤与版本

1. 原始POD5 → **Dorado 1.3.0**，模型 **`rna004_130bps_hac@v5.2.0`**，basecall必须保留moves（`--emit-moves`）。POD5读取环境记录为 **pod5 0.3.15**、Python3.9.20。
2. 同套 **Dorado aligner＋samtools** 比对、过滤、排序、索引。原脚本用 `samtools view -F 2308 -q 10` 保留主比对且MAPQ≥10；历史samtools版本**未锁定**，新运行须保存 `samtools --version`。
3. **ont-remora 3.3.0** 的 `SigMapRefiner`＋包内 `9mer_levels_v1.txt`，配合 **pysam 0.22.1**，运行 [feature_extraction_rna004.py](upstream_reference/RNA004/preprocess_pipeline/feature_extraction_rna004.py)（`--kmer 9 --step 1 --resquiggle True`）→ NPZ＋read/reference/start/end/strand sidecar。
4. [merge_step6.py](upstream_reference/RNA004/preprocess_human_pipeline/merge_step6.py) → Step6 CSV，含 `read_id/ref_name/ref_start/ref_end/strand/sequence_input` 和36个统计值。**此时sequence_input仅9 nt，Step6还不是prepared。**
5. 将连续9行重叠9-mer构造17 nt，保留中间行36个值。下面的旧入口**只适用于全部显式strand=+、符合相邻重叠条件的表**：

```bash
nanoremix prepare --input sample004.positive.step6.csv \
  --output prepared004.csv --chemistry RNA004 --format step6
```

**负链/剪接处理边界：**含任何负链、混合链、缺失/无效strand或同read多个参考时，旧入口中止且不写输出。不要通过改链标签或丢掉负链冒充完整流程。负链数据需经过验证的适配器联合构造RNA方向17 nt及数值槽位；已经正确prepared的正/负链可正常predict。现存catalog适配脚本依赖已知truth，不能原样当作未知样本通用转换器。

融合还需同中心81 nt；跨内含子用成熟转录本或BAM CIGAR剪接路径，不能取普通基因组相邻81 nt代替。全部命令、regions格式和适配限制见[RNA004 pipeline](docs/guide_20261104/PIPELINE_命令与输入衔接.md)及[预处理详细说明](#preprocessing)。

### 后续用NanoRemix运行

解压 **RNA004_Binary** 模型包；融合再解压 **SecondGen** 包。

```bash
# 纯三代
nanoremix predict --input prepared004.csv --output results/004_binary_third \
  --chemistry RNA004 --task binary --mods m6A,m5C,pseudoU \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files

# 二代＋三代融合
nanoremix predict --input prepared004_with81.csv --output results/004_binary_fused \
  --chemistry RNA004 --task binary --mods m6A,m5C,pseudoU \
  --data-type biological --sequence-column sequence81 --strategy s2 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

<a id="rna004-multi"></a>

## 方式四：RNA004 Multi

### 用途与修饰

用于RNA004的多类预测：模型原生**16类＝12种修饰＋A/C/G/U**，在中心参考碱基对应家族内比较所选类别。

| 中心参考家族 | 支持的修饰（另含对应正常对照） |
|---|---|
| A | m6A、m1A、Am |
| C | m5C、hm5C、ac4C |
| G | o8G、Gm |
| U | pseudoU、m1psi、m5U、Um |

`--mods all`选全部上述修饰，也可用 `--mods m6A,m1A`等子集。此配置的三代模型没有Inosine、FormylC、m7G或Cm。默认extended融合保留12种修饰，但**hm5C、o8G、m1psi只有三代证据**。

### 输入是什么样

| 内容 | RNA004 Multi要求 |
|---|---|
| `9mer_final` | 17 nt，例如 `A_C_G_U_A_C_G_U_A_C_G_U_A_C_G_U_A`，中心第9位 |
| 数值列 | mean/median/std/length各9列，列名同RNA004 Binary，共36个值 |
| 融合额外输入 | 同候选、同RNA方向的连续81 nt `sequence81`，或完整15类二代概率 |
| 格式文件 | [prepared＋81 nt人工CSV](docs/guide_20261104/examples/RNA004_prepared_with81.csv)；与RNA004 Binary共享列格式，使用不同模型 |

### 输入怎样得到：上游步骤与版本

流程同RNA004 Binary，**不需要先运行Binary**：

1. POD5（**pod5 0.3.15**）→ **Dorado 1.3.0 / `rna004_130bps_hac@v5.2.0`**，用 `--emit-moves` 输出basecalled BAM。
2. **Dorado aligner＋samtools** → 主比对过滤（`-F 2308 -q 10`）、sort/index；samtools历史版本未锁定，实际运行另存版本号。
3. **ont-remora 3.3.0＋pysam 0.22.1＋`9mer_levels_v1.txt`**，运行包内 `feature_extraction_rna004.py --kmer 9 --step 1 --resquiggle True` → NPZ及坐标sidecar。
4. `merge_step6.py` → 9 nt Step6＋36值 → 正确RNA方向的17 nt prepared；保留中间候选的36值，不把RNA002的MAD归一化再次套到RNA004。

全部正链且符合旧入口条件时：

```bash
nanoremix prepare --input sample004.positive.step6.csv \
  --output prepared004.csv --chemistry RNA004 --format step6
```

**同样受RNA004旧Step6限制：**负链/混合链或strand缺失会中止；正确prepared的负链可predict，但本包尚无任意未知样本的通用负链适配器。融合增加同中心连续81 nt，跨剪接时走转录本或正确CIGAR路径。详见[RNA004上游命令与输入衔接](docs/guide_20261104/PIPELINE_命令与输入衔接.md)。

### 后续用NanoRemix运行

解压 **RNA004_Multi** 模型包；融合再解压 **SecondGen** 包。

```bash
# 纯三代：12种修饰
nanoremix predict --input prepared004.csv --output results/004_multi_third \
  --chemistry RNA004 --task multi --mods all \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files

# 二代＋三代融合，并保存两代分支
nanoremix predict --input prepared004_with81.csv --output results/004_multi_fused \
  --chemistry RNA004 --task multi --mods all \
  --data-type biological --sequence-column sequence81 --strategy s2 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

## 四种方式共用：选修饰、换二代输入、看结果

### 修饰标签怎样写

命令使用上面表中的拼写，逗号列表中不加空格，例如 `--mods m6A,m5C,pseudoU`；正常碱基不要写进mods。默认 `--canonical matching` 自动加入对应正常对照。Multi可改用 `--targets A,m6A,m1A` 精确指定竞争类，**mods与targets二选一**；targets不自动补对照。

常用标签含义如下；是否支持以各方法清单为准，不能把本表所有标签传给任意配置。

| 标签 | 修饰名称 |
|---|---|
| m6A / m1A | N6-甲基腺苷 / N1-甲基腺苷 |
| Inosine | 肌苷（A-to-I相关类别） |
| Am / Gm / Um | 2′-O-甲基腺苷 / 2′-O-甲基鸟苷 / 2′-O-甲基尿苷 |
| m5C / hm5C / FormylC | 5-甲基胞苷 / 5-羟甲基胞苷 / 5-甲酰基胞苷 |
| ac4C | N4-乙酰胞苷 |
| m7G / o8G | N7-甲基鸟苷 / 8-氧代鸟苷 |
| pseudoU / m1psi / m5U | 假尿苷（Ψ）/ N1-甲基假尿苷 / 5-甲基尿苷 |

### 不用sequence81时怎样融合

上面四条融合命令均可替换二代来源，其余chemistry/task/mods不变：

| 已有输入 | 替换 `--sequence-column sequence81` 的参数 | 还需SecondGen权重吗 |
|---|---|---|
| 独立二代概率表 `second.tsv` | `--second-probabilities second.tsv --join-on record_key` | 不需要 |
| 主表内嵌二代概率 | 提供全部 `second_prob_<label>` 列，不传序列/独立表参数 | 不需要 |
| 对应参考FASTA及正确中心坐标 | `--context-fasta transcripts.fa --fasta-id-column reference_name --center-column center_pos0 --coordinate-base 0`；按FASTA实际方向决定是否加 `--strand-column strand` | 需要 |

两种概率输入都必须有完整15类二代轴；独立表按唯一身份键匹配，不能只填待测修饰概率或按行号拼接。FASTA中心坐标属于所指定记录，工具不自动转换基因组/转录本坐标，也不自动从BAM解读剪接路径。已按转录方向的FASTA不要再次反向互补。详见[输入合同](#inputs)。

### 结果可以做什么、如何筛选

| 文件（相对于指定output） | 用途 |
|---|---|
| `thirdgen/s2/predictions.tsv` | 纯三代的逐read×候选预测及概率 |
| `fused/s2/predictions.tsv` | 融合预测及证据来源 |
| `secondgen/none/predictions.tsv` | 二代支持类别的分支结果，需请求该输出 |
| 各路线的 `selection_audit.tsv` | 全部输入记录、保留/拒绝状态及原因 |
| `run_manifest.json` | 参数、输入/模型hash、类别、计数和可选评价 |

S1输出在相应s1目录；`--strategy both`分别写s1/s2。查看 `status`、`predicted_label`、`candidate_label`、`confidence`、类别概率和 `evidence_source`，可比较纯三代/二代/融合，并按样本、read、候选身份进行后续分析。

可按验证数据预先确定 `--min-confidence`、`--min-selected-mass`、`--min-second-mass`（默认均0，严格低于阈值拒绝）；`--output-filter modified`只输出有效修饰预测。阈值默认0仍做格式、中心、方向和有效支持检查。confidence是所选家族类别中的条件概率，不是位点修饰率。**CLI不自动做motif/MAPQ/coverage筛选或read-to-site修饰率汇总**；这些需要显式上游/下游步骤。

`--fusion-scope extended`（默认）保留三代独有类别并标记证据；`--fusion-scope common`仅比较共同支持的类别。完整计算和所有参数见[参数表](#options)与[公式](#calculation)。

<a id="quickstart"></a>

## 首次检查：无需权重的人工示例

在上述源码目录和已安装环境中运行，不需要下载模型：

```bash
nanoremix predict --input examples/minimal_binary_probabilities.csv \
  --output results/quickstart --chemistry RNA004 --task binary --mods m6A \
  --input-kind probabilities --mode thirdgen --strategy s2
```

预期两行结果分别为m6A（prob_m6A=0.90）和A（prob_A=0.80），输出在 `results/quickstart/thirdgen/s2/predictions.tsv`。人工概率用于软件检查，不代表模型准确率。

## 详细参考

下面保留完整输入、预处理版本、全部41个predict/validate长选项、公式及模型布局。v0.3.1没有改变预测数学算法或v0.3.0权重。

<a id="configurations"></a>

## 1. 版本数量与名称

NanoRemix 统一推理工具版本为 **v0.3.1**；继续使用v0.3.0的原始权重。这里有 **4 套第三代 chemistry/task 模型配置**，并非 4 个不同的软件安装包：

| 编号 | chemistry | task | 三代模型输入 | 三代预训练对象 |
|---|---|---|---|---:|
| 1 | RNA002 | Binary | 9 nt 序列上下文 + 20 个信号数值 | 18 个 `.pkl` |
| 2 | RNA002 | Multi | 9 nt 序列上下文 + 20 个信号数值 | 6 个 `.pkl`，另有 S1/S2 manifests |
| 3 | RNA004 | Binary | 17 nt 序列上下文 + 36 个信号数值 | 30 个 `.pkl` |
| 4 | RNA004 | Multi | 17 nt 序列上下文 + 36 个信号数值 | 10 个 `.pkl`，另有 S1/S2 manifests |

四套三代配置都能运行纯三代或二代+三代融合，共 **8 条主要预测路线**。每条可选 S1 或 S2，所以有 16 种 route/strategy 配置；`--strategy both` 是一次分别跑 S1 和 S2，不计算两者共识。另有共享的二代序列模型，可单独输出二代结果。所有路线由同一个 `nanoremix` CLI 管理。

“二代/三代”是项目内部模型分支名称，不代表测序技术代际。二代使用 81 nt 连续序列；三代使用相应 chemistry 的序列和 Nanopore 电流特征。**Binary/Multi** 是分类任务类型，**RNA002/RNA004** 是 RNA 化学版本，**S1/S2** 是模型策略；这些维度不能互相替代。

旧模型文件名 `strategy1` 对应 CLI `s1`，`strategy3` 对应 CLI `s2`。CLI 接受旧别名 `--strategy s3` 并将其解释为 `s2`。不应把 S1/S2 称为二代/三代融合；两代融合是另一层独立操作。

## 2. 四套配置与完整类别的对应

四种方式的用途、修饰、输入、上游版本和成对运行命令已分别列在上方：[RNA002 Binary](#rna002-binary)、[RNA002 Multi](#rna002-multi)、[RNA004 Binary](#rna004-binary)、[RNA004 Multi](#rna004-multi)。同一chemistry的Binary/Multi共享prepared格式，不能共用任务权重。

共享二代网络的15类固定为 `A,C,G,U,Am,m1A,m6A,Inosine,ac4C,m5C,Gm,m7G,m5U,pseudoU,Um`。默认 `--fusion-scope extended` 保留三代独有类别：RNA002的hm5C、FormylC、o8G、m1psi；RNA004的hm5C、o8G、m1psi。这些类别只有三代证据。`common`只使用共有轴（RNA002 Multi为14类，RNA004 Multi为13类，均包括A/C/G/U）。`Cm`不在正式模型范围；`nanoremix catalog`可查询原生标签及别名。

## 3. 纯三代、二代+三代融合和纯二代

| 路线 | 主输入 | 第二数据来源 | 必需模型资产 |
|---|---|---|---|
| 纯三代 | 对应 RNA002/RNA004 的 prepared 特征表，或三代预计算概率 | 无 | 对应 chemistry/task 的三代包；若主输入是三代概率，则无需加载三代模型 |
| 二代+三代融合（现场跑二代） | 同一条三代 prepared 特征/概率记录 | 同候选、同方向的连续 81 nt 序列，或能按身份和坐标抽取该序列的 FASTA | 对应三代包 + 共享 SecondGen 包（25 个 checkpoint） |
| 二代+三代融合（已有二代概率） | 同一条三代 prepared 特征/概率记录 | 完整 15 类二代概率，内嵌或按唯一身份键连接 | 对应三代包；不用下载 SecondGen checkpoint |
| 纯二代 | 中心参考碱基与连续 81 nt，或完整二代概率 | 无 | 现场推理时用 SecondGen；二代概率输入时不需 checkpoint |

生物样本的默认路线为融合，必须能取得两代证据；工具不会悄悄降级成纯三代。生物样本纯三代请显式使用 `--mode thirdgen`。`--data-type synthetic` 表示合成数据常用的纯三代路由，不要为规避二代输入把真实样本标成 synthetic。

使用 `--data-type biological` 的融合运行可以加 `--also-thirdgen --also-secondgen` 同时保存两条组成分支；不能把这两个选项加在显式 `--mode fused/both` 后。`--mode both` 表示分别跑纯三代和融合（仍须二代证据）；不包含纯二代模式。`--mode secondgen` 单独运行二代。`--data-type` 与 `--mode` 互斥。

纯二代路线不读取三代信号特征，不需要三代权重；输入要有中心碱基与连续 81 nt 序列，或可替代的完整二代概率。其类别仅限共享 15 类二代轴。二代-only 常规调用使用 `--task multi`；需要同时评估三代和二代时，更常见的用法是融合路线。

<a id="inputs"></a>

## 4. 输入文件合同

主表为 UTF-8 CSV/TSV，可使用 `.csv.gz`/`.tsv.gz`。一行对应一个 read 覆盖的候选中心，不是已经汇总的位点。推荐给 `record_key` 一个全表唯一、非空值；额外保留 `sample_id`、`read_id`、`site_id`、参考名称、坐标和链方向方便追溯。不同样本中可能重复的 read/site ID 应加样本前缀。跨表连接按明确的 `record_key` 或复合键，不按行号连接。

### 4.1 特征输入

`--input-kind features` 是默认输入类型。

| 列 | RNA002 | RNA004 |
|---|---|---|
| `9mer_final` | 9 个碱基 token，中心第 5 位 | 17 个碱基 token，中心第 9 位；历史列名仍叫 `9mer_final` |
| `mean_output1..K` | K=5 | K=9 |
| `median_output1..K` | K=5 | K=9 |
| `std_output1..K` | K=5 | K=9 |
| `length_output1..K` | K=5 | K=9 |
| 信号数值合计 | 20 | 36 |

Binary 和 Multi 在相同 chemistry 下使用相同 prepared 列合同。可提供 `center_ref_base` 并要求它与上下文中心一致。序列上下文由正常参考碱基构成。数值必须有限，std 和 length 非负；使用匹配 chemistry 原定义的归一化，不能随意再做标准化。length 表示事件/dwell 长度统计，不是核酸长度。每行信号属于该候选中心，不能先跨多个 read 取均值后输入。

### 4.2 三代预计算概率输入

`--input-kind probabilities` 可跳过三代模型特征推理，但不能跳过任务对应的概率格式检查。

| Task | 必需列 |
|---|---|
| Binary | 所选修饰与策略的 pair 概率，如 `s1_prob_m6A`、`s2_prob_m6A`；值表示专用 Binary pair 中的 P(mod)，不是 Multi 原始分数 |
| RNA002 Multi | 每种策略完整 18 类原生轴，如 `s2_prob_A` 加全部其他 RNA002 标签 |
| RNA004 Multi | 每种策略完整 16 类原生轴 |

Multi 完整轴每行概率和应为 1。`--strategy both` 时须提供两种策略的概率列。输入仍需中心正常碱基信息；融合还需一种有效的二代来源。

### 4.3 二代序列或概率

现场二代推理可在主表添加 `sequence81`，并传 `--sequence-column sequence81`；也可使用 `--context-fasta`，同时指定 `--fasta-id-column`、`--center-column`、`--coordinate-base 0|1`，需要时再指定 `--strand-column`。工具从指定 FASTA 记录的中心取前后各 40 nt。81 nt 必须是连续 A/C/G/T（U 转为 T）、中心第 41 位，不能含 N；方向须与候选一致。RNA002 的中央 9 nt 必须与 prepared 上下文一致；RNA004 的中央 17 nt 必须一致。跨剪接边界时需提供成熟转录本或沿正确 BAM CIGAR 路径提取的上下文；普通基因组连续窗口不会自动跳过内含子。

已有二代概率时可在主表提供完整的 `second_prob_<label>` 列，也可用 `--second-probabilities second.tsv --join-on record_key` 读取独立表。需要完整二代 15 类轴、有限非负值且总和约为 1；独立表键非空且唯一。提供序列列时可用于检查身份/上下文，但不会重新运行二代网络。

<a id="preprocessing"></a>

## 5. 原始信号到 NanoRemix 输入的准备流程

NanoRemix 不从 FAST5/POD5 直接 basecall，也不直接提取电流特征。原始信号预处理与 NanoRemix 现代推理环境分开运行。下面的软件版本来自超算已有配置、环境元数据和日志；不能推断为每个历史作业都具有完整锁定环境。

### 5.1 RNA002：multi-FAST5 到 9 nt/20 值

1. 使用 **Guppy 6.3.8+d9e0f64**、`rna_r9.4.1_70bps_hac.cfg`、FLO-MIN106/SQK-RNA002 basecall，并保留 basecalled FAST5。
2. 用 **ont-fast5-api 4.1.1** 的 `multi_to_single_fast5` 拆成 single-read FAST5。已有配置环境是 Python 3.6.15。
3. 用 **Tombo 1.0** 执行 resquiggle，配合 **GraphMap 0.5.2** 和匹配的参考 cDNA；Tombo 是独立 Python 2.7 环境。检查 corrected Events/Alignment 和失败记录，不能只依据退出码。
4. 运行上游 `extract_tombo_step4.py`，提取连续 5 个事件的 5-mer 信号统计和同步 `.coordinates.tsv` sidecar。已有特征环境记录为 Python 3.6.9、NumPy 1.19.2、h5py 2.10.0。现存作业使用 `--clip 10`；脚本的 `--rna` 是提取器参数，不是 Tombo 参数。
5. 得到 Step4 TSV 后运行 `nanoremix prepare --input STEP4.tsv --format step4 --chemistry RNA002 --output prepared.csv`。每 5 个连续重叠的 5-mer 构成 9 nt；信号取中间行的 20 个值，不对 5 行再次平均。不要预先删掉邻接事件，否则候选窗口断裂。
6. 纯三代可直接用 prepared。融合时需用 sidecar 的 `id` 显式连接 read、链、转录本 ID 和中心坐标，再从匹配的参考转录本或完整 corrected-event 序列得到同中心连续 81 nt。event index 不是转录本/基因组坐标；不能把剪切后的片段序列与旧 event index 混用。

RNA002 软件环境需拆开：Guppy/FAST5、Tombo、特征提取和 NanoRemix 推理依赖彼此较老/不同，不能把 Tombo 1.0 与现代模型环境强行装在同一个 Python 环境。

### 5.2 RNA004：POD5 到 17 nt/36 值

1. 使用 **pod5 0.3.15**（现存环境 Python 3.9.20）读取 POD5。
2. 用 **Dorado 1.3.0** 和 `rna004_130bps_hac@v5.2.0` basecall，输出 move table（`--emit-moves`）。
3. 用相同参考的 Dorado aligner 比对，并用 samtools 过滤、排序和索引。现存脚本的 `samtools view -F 2308 -q 10` 排除 unmapped、secondary、supplementary alignment 并要求 MAPQ≥10；这是上游脚本设置，不是 NanoRemix CLI 隐式过滤。历史配置没有锁定 samtools 版本，新流程需记录 `samtools --version`。
4. 用 **ont-remora 3.3.0** 的 `SigMapRefiner` 和匹配的 `9mer_levels_v1.txt` 做 refined signal mapping/9-mer 特征，保留 NPZ 与逐行坐标 sidecar。现存特征环境记录 pysam 0.22.1、NumPy 1.24.4、pandas 2.3.3、h5py 3.11.0、torch 2.4.1。
5. 对齐核实 NPZ 与 sidecar 行数后由 `merge_step6.py` 汇总 Step6。中心基因组坐标为 `ref_start+4`。Remora 的数值字段映射为 mean、std、median、length；输出还应保留 read、reference、start/end、strand 和候选身份。
6. 构造中心正确、RNA 方向一致的 17 nt 上下文并保留中央 Step6 的 36 个信号数值。融合时另需相同候选的 81 nt 连续成熟 RNA 上下文；跨内含子应沿剪接路径提取。

**RNA004 Step6转换限制：**旧v0.3.0转换器按参考坐标升序拼接相邻9-mer，对RNA方向负链可能丢候选。v0.3.1的 `nanoremix prepare --format step6` 要求显式 `strand` 列，值为+/-；若包含任意负链行（包括正负链混合），在写出结果前明确报错。同一read若声明多个参考也会报错，避免跨比对拼接。不要改写链标签或删除负链来冒充完整处理；此保护没有实现通用负链适配。已有 `prepare_model_input.py` 和 splice-aware 上下文脚本依赖已知 catalog/truth 字段，也不能原样用于未知样本。新样本需按 `ref_start+4` 构造上下文、校验中央序列、同步反转序列与数值槽位（需要反向时），并记录每条拒绝原因；当前包没有声称包含经过验证的任意未知样本全转录组适配器。

RNA004 的参考序列、regions、BAM 和 FASTA 必须对应同一参考版本。9-mer region 文件格式和区间是候选 9-mer 左端起点，不等同于任意 site BED；细节见 [RNA004 pipeline 模板](docs/guide_20261104/PIPELINE_命令与输入衔接.md)。

### 5.3 推理环境

原始权重兼容环境建议 Python 3.10/3.11，scikit-learn 1.3.0、XGBoost 2.1.4。项目的 `constraints-models.txt` 给出安装约束。basecalling/信号处理与推理使用独立环境；包内没有包含 Guppy、Dorado、Tombo、Remora、参考基因组或数据。

<a id="options"></a>

## 6. 全部predict / validate选项与条件

`--mods` 与 `--targets` 二选一必填；`--data-type` 与 `--mode` 二选一必填，不能都省略。下面包括帮助在内共41个长选项。

`validate`与`predict`共享以下解析选项，但通常不传output，且不保存预测。完整原始help见 `nanoremix predict --help` 或 `docs/guide_20261104/CLI_help.txt`。

| 参数 | 默认/允许值 | 用法 |
|---|---|---|
| `--input` | 必需 | CSV/TSV，按后缀分隔，可gzip |
| `--output` | predict必需 | 新输出目录；已存在则报错 |
| `--task` | 必需；binary/multi | 任务类型 |
| `--chemistry` | auto；RNA002/RNA004 | 建议显式与真实实验一致 |
| `--mods` | 与targets二选一必填 | 逗号分隔修饰，或all |
| `--targets` | Multi专用 | 精确类别列表，与mods互斥 |
| `--canonical` | matching/all，默认matching | mods路线的正常对照 |
| `--data-type` | biological/synthetic | 与mode二选一必填；biological默认融合 |
| `--mode` | thirdgen/fused/both/secondgen | 显式路由；both为纯三代＋融合，不含纯二代 |
| `--also-thirdgen` | 默认关 | 仅data-type biological可用 |
| `--also-secondgen` | 默认关 | 同上；只输出二代支持的比较 |
| `--strategy` | s1/s2/both，默认s2 | s3可作s2旧名；both分别输出 |
| `--input-kind` | features/probabilities，默认features | 特征或预计算三代概率 |
| `--condition` | full/mean/std/length，默认full | ablation仅Binary；完整ZIP只带full权重，其他条件另需对应权重 |
| `--model-root` | 无自动默认路径 | 包含RNA002、RNA004、secondgen的根；特征推理需给 |
| `--model-checksums` | 默认无 | 建议 `models/SHA256.json`，校验实际加载资产 |
| `--trust-model-files` | 默认关 | 实际加载已核实的pickle、checkpoint和模型Python代码所必需 |
| `--second-model-root` | 默认model-root/secondgen | 覆盖二代根；启用hash检查时须仍在model-root下 |
| `--second-probabilities` | 默认无 | 独立完整15类二代表 |
| `--join-on` | 独立表时必需 | 逗号分隔身份键 |
| `--sequence-column` | 默认无 | 主表81 nt列；与独立二代表并用时验证独立表的序列 |
| `--context-fasta` | 默认无 | 完整参考/事件序列FASTA，可gzip |
| `--fasta-id-column` | FASTA时必需 | 匹配FASTA header第一个token的输入列 |
| `--center-column` | FASTA时必需 | 所命名FASTA记录中的中心坐标列 |
| `--coordinate-base` | FASTA时必需；0/1 | 明确0基或1基 |
| `--strand-column` | 默认不反向互补 | 有则值必须为+/-，负向窗口反向互补 |
| `--fusion-scope` | extended/common，默认extended | Multi保留三代独有类或仅共有类 |
| `--second-weight` | 1 | 二代融合指数；此CLI要求有限且严格大于0 |
| `--third-weight` | 1 | 三代融合指数；同上，不要求两者和为1 |
| `--binary-bias` | 0 | 仅Binary融合，log-odds偏置；Multi或非融合非零值报错 |
| `--min-confidence` | 0 | [0,1]，最高条件概率下限 |
| `--min-selected-mass` | 0 | [0,1]，所选家族质量下限 |
| `--min-second-mass` | 0 | [0,1]，二代所选家族质量下限 |
| `--output-filter` | all/retained/modified，默认all | 全记录/仅status=ok/仅ok且预测为修饰 |
| `--on-invalid` | error/mark，默认error | 支持的行级错误报错或标记；schema/身份/模型错误仍中止 |
| `--truth-column` | 默认无 | 仅评价；不参与模型或融合参数拟合 |
| `--chunk-size` | 50000 | predict主表分块；键集合/FASTA/独立二代表等仍驻内存，不是全流程内存上限 |
| `--batch-size` | 256 | 二代序列推理批大小 |
| `--device` | cpu/cuda/mps，默认cpu | 二代设备；需本机兼容，不代表所有设备均经过验证 |
| `--keep-features` | 默认关 | 保留额外输入列，排除各类输入概率列；默认仅常用元数据 |
| `-h/--help` | — | 命令帮助 |

完整性注意：`--second-weight 0` 或 `--third-weight 0` 在底层某些数学函数可计算，但**本CLI会拒绝**，因此不用零权重模拟纯代；请明确使用mode。

`--condition mean/std/length` 是消融模型条件，并非额外的一套chemistry、策略或数据类型。只删除列而继续用full权重是不正确的；本ZIP没有附带这三组消融权重。

`--also-thirdgen`、`--also-secondgen` 只用于 `--data-type biological`；`--strategy both` 不生成S1/S2共识。实际参数是 `--condition`，没有 `--feature-set`。

先按所选标签和中心参考碱基家族条件化、再归一化和argmax。例如只选m6A/m5C时，U/G中心通常标为 `outside_selected_family`；`--canonical all` 额外保留其他正常对照。Multi的targets精确指定竞争类，不自动补正常碱基。

三个质量阈值默认0，范围[0,1]，严格低于阈值才拒绝；默认0仍执行格式、中心/方向和有效支持检查。confidence是获胜类别的条件概率，不一定是P(mod)，也不是位点修饰比例；selected mass不是read覆盖度。遗漏真实修饰可能把概率重新分配给选中的类别。

工具不自动做motif、MAPQ、coverage、位点阈值、S1/S2共识或read-to-site修饰率汇总。上游/下游筛选应明确执行和记录。

<a id="calculation"></a>

## 7. 主要计算方法

### Binary 三代与融合

S1 的多个校准 pair 分数按 logit 相加：

```text
p3 = sigmoid(sum_j(logit(p3_j)))
```

S2 使用单个联合模型分数。二代完整 15 类概率先按目标正常碱基/修饰 pair 条件化，再按默认融合权重和偏置合并：

```text
pf = sigmoid(w2*logit(p2) + w3*logit(p3) + beta)
```

默认 `w2=w3=1, beta=0`；例如 p2=0.8、p3=0.6，结果为 0.857142857。logit 概率会在计算时进行数值边界保护。

### Multi 三代内部聚合与两代融合

三代 S1/S2 的多个组件概率按模型 manifest 中的 class axis、家族权重和 temperature 聚合。随包实际 manifest 使用冻结的 geometric/gate 参数；不要把它改写成简单等权平均。具体值见 `guide_20261104/frozen_multi_fusion.json`。

融合时先在研究标签和参考碱基家族内条件化概率。对二代与三代共有类别，按默认权重 1/1 做归一化概率幂乘积；`extended` 将融合结果乘回共有类别原来占有的三代质量，三代独有类别保持原三代概率。完全没有有效共有正概率支持时会标为无效，不造均匀分布。`common` 先限制到共有轴再融合。二代不支持的类别必须通过 `evidence_source` 识别为三代证据。

### 阈值与决策顺序

先检查 schema、chemistry、中心/方向、数值和身份，再确定类别/家族并计算概率；接着计算 selected mass/second mass、argmax 与 confidence；只有 status 为 `ok` 才有接受的 `predicted_label`。阈值拒绝通常保留 `candidate_label` 和审计原因，但 accepted label 为空。精确公式和数值案例见 [详细计算章节](docs/guide_20261104/README_版本输入预处理参数与计算.md)。

<a id="models"></a>

## 8. 五个权重 Release 附件

从本仓库 [Releases](https://github.com/iris982022/nanoremix_tool/releases) 下载所需附件；如果页面没有附件，则相应模型的公开发布尚未完成。纯三代只取对应的一包；现场执行序列型二代融合时再取共享 SecondGen 包；已有二代概率表不需要 SecondGen。

| Release asset | 包含内容 | 用途 |
|---|---|---|
| `NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip` | 18 个 `.pkl` | RNA002 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip` | 6 个 `.pkl`、S1/S2 manifests | RNA002 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip` | 30 个 `.pkl` | RNA004 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip` | 10 个 `.pkl`、S1/S2 manifests | RNA004 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_SecondGen_20261104.zip` | 网络定义及共享 25 个 `.pt` checkpoint | 二代推理；RNA002/RNA004 共用 |

四组三代共 64 个 `.pkl`，共享二代共 25 个 `.pt`，完整集合为 89 个权重。Binary 消融模型未包含。模型包解压到项目根目录后组成 `models/`。包内源码没有模型权重，不需要重新训练。

<a id="examples"></a>

## 9. 环境、预检和辅助命令

模型推理依赖约束在 `constraints-models.txt`：numpy1.26.4、pandas2.3.1、scipy1.11.4、scikit-learn1.3.0、XGBoost2.1.4、joblib1.5.1、torch2.7.1。安装命令见[安装与模型](#start-install)；四套配置的纯三代/融合命令见[逐一使用](#start-run)。

```bash
nanoremix --version
nanoremix catalog --chemistry RNA002
nanoremix catalog --chemistry RNA004
nanoremix predict --help
nanoremix validate --help
nanoremix doctor --help
nanoremix prepare --help
```

| 命令 | 功能 |
|---|---|
| `catalog` | 查询chemistry、修饰和类别别名 |
| `doctor` | 预检依赖、文件、模型manifest和hash，不运行真实模型；只有data-type选项，没有mode |
| `prepare` | 转换已有Step4/extracted/Step6特征为prepared，不执行basecall/比对/信号提取 |
| `validate` | 使用predict相同输入/模型参数执行推理检查，通常不传output，不保存预测表 |
| `predict` | 保存逐候选预测、筛选审计、运行清单及可选truth评价 |

`prepare`必需input、output、chemistry和format；RNA002支持step4/extracted，RNA004仅支持满足显式正链要求的step6，输出须为新文件。`--chunk-size`只控制predict主表批量，不是prepare或全流程的内存上限；大表按完整read分片并保持全局身份唯一。

<a id="references"></a>

## 10. 资料导航和适用范围

- [详细中文方法手册](docs/guide_20261104/README_版本输入预处理参数与计算.md)：全部 41 个 predict/validate 选项、字段、版本证据、完整筛选顺序和公式例子。
- [RNA002/RNA004 上游命令模板](docs/guide_20261104/PIPELINE_命令与输入衔接.md)：命令、Step4/Step6 输入衔接和坐标示例。
- [输入列定义](docs/INPUTS.md)、[模型布局与依赖](docs/MODELS.md)、[CLI 选项说明](docs/OPTIONS.md)、[限制说明](docs/KNOWN_LIMITATIONS.zh-CN.md)。
- [人工输入示例目录](docs/guide_20261104/examples/README.md)：所有示例均为人工数据，不是准确率 benchmark。
- [原始上游副本索引](upstream_reference/README.md) 与 [`model_path_mapping.tsv`](docs/guide_20261104/model_path_mapping.tsv)：上游脚本/权重路径对应。

包完整性通过 `python verification/verify_distribution.py` 检查只代表文件哈希、清单和缺失文件通过；它不运行模型、不重做原始数据处理，也不证明新样本准确率。源码许可状态见 [LICENSE](LICENSE)，模型来源和分发元数据见 [NOTICE.md](NOTICE.md)。
