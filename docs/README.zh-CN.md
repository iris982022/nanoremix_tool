# NanoRemix 完整使用手册

NanoRemix v0.3.1 使用 v0.3.0 发布的模型权重，支持 RNA002 和 RNA004 的逐 read 候选修饰预测。以下步骤从已有 FAST5/POD5 开始。手头已有 prepared 特征表或概率表时，可直接从[第 3 步](#inputs)开始。

[返回仓库首页](../README.md) · [下载发布权重](https://github.com/iris982022/nanoremix_tool/releases/tag/v0.3.0) · [已知限制](KNOWN_LIMITATIONS.zh-CN.md)

1. [选定配置和预测路线](#choose)
2. [完成预处理](#preprocessing)
3. [整理输入文件](#inputs)
4. [安装并检查模型](#start-install)
5. [分路线运行预测](#start-run)
6. [检查和分析结果](#outputs)
7. [所有 CLI 参数](#options)
8. [计算方法与参考资料](#calculation)

<a id="choose"></a>
<a id="configurations"></a>
<a id="modifications"></a>

## 1. 选定配置和预测路线

### 1.1 按实验试剂和研究问题选配置

| 配置 | 试剂 | 修饰范围 | 三代 prepared 输入 |
| --- | --- | --- | --- |
| [RNA002 Binary](#rna002-binary) | SQK-RNA002 | m6A、m5C、pseudoU | 9 nt + 20 个信号值 |
| [RNA002 Multi](#rna002-multi) | SQK-RNA002 | 14 种修饰 + 4 个未修饰类 | 9 nt + 20 个信号值 |
| [RNA004 Binary](#rna004-binary) | SQK-RNA004 | m6A、m5C、pseudoU | 17 nt + 36 个信号值 |
| [RNA004 Multi](#rna004-multi) | SQK-RNA004 | 12 种修饰 + 4 个未修饰类 | 17 nt + 36 个信号值 |

Binary 分别比较 A/m6A、C/m5C、U/pseudoU 三组配对。可同时请求多个配对，候选仅进入参考中心对应的配对。Multi 在每个候选的参考碱基家族内比较所选类别。

| 参考中心 | RNA002 Multi 的修饰 | RNA004 Multi 的修饰 |
| --- | --- | --- |
| A | m6A、m1A、Inosine、Am | m6A、m1A、Am |
| C | m5C、hm5C、ac4C、FormylC | m5C、hm5C、ac4C |
| G | m7G、o8G | o8G、Gm |
| U | pseudoU、m1psi、m5U、Um | pseudoU、m1psi、m5U、Um |

用 `--mods m6A,m5C` 提供修饰列表，或用 `--mods all` 请求当前路线支持的全部修饰。默认加入同家族的未修饰对照。Multi 也可改用 `--targets A,m6A,m1A` 精确列出参与竞争的类别。

| 标签 | 修饰名称 |
|---|---|
| m6A / m1A | N6-甲基腺苷 / N1-甲基腺苷 |
| Inosine | 肌苷（A-to-I相关类别） |
| Am / Gm / Um | 2′-O-甲基腺苷 / 2′-O-甲基鸟苷 / 2′-O-甲基尿苷 |
| m5C / hm5C / FormylC | 5-甲基胞苷 / 5-羟甲基胞苷 / 5-甲酰基胞苷 |
| ac4C | N4-乙酰胞苷 |
| m7G / o8G | N7-甲基鸟苷 / 8-氧代鸟苷 |
| pseudoU / m1psi / m5U | 假尿苷（Ψ）/ N1-甲基假尿苷 / 5-甲基尿苷 |

### 1.2 选路线和策略

“二代”和“三代”指 NanoRemix 内部的两个模型分支。二代使用 81 nt 序列，三代使用序列和 Nanopore 电流特征。

| 路线 | 主输入 | 第二数据来源 | 必需模型资产 |
|---|---|---|---|
| 纯三代 | 对应 RNA002/RNA004 的 prepared 特征表，或三代预计算概率 | 无 | 对应 chemistry/task 的三代包；若主输入是三代概率，则无需加载三代模型 |
| 二代+三代融合（现场跑二代） | 同一条三代 prepared 特征/概率记录 | 同候选、同方向的连续 81 nt 序列，或能按身份和坐标抽取该序列的 FASTA | 特征输入需对应三代包，另需共享 SecondGen 包（25 个 checkpoint） |
| 二代+三代融合（已有二代概率） | 同一条三代 prepared 特征/概率记录 | 完整 15 类二代概率，内嵌或按唯一身份键连接 | 特征输入需对应三代包；预计算概率可直接使用 |
| 纯二代 | 中心参考碱基与连续 81 nt，或完整二代概率 | 无 | 现场推理时用 SecondGen；二代概率输入时不需 checkpoint |

`--data-type biological` 选择融合并要求两代证据，`--data-type synthetic` 选择纯三代。生物样本只做三代时，显式指定 `--mode thirdgen`。

要随融合保存两条组成分支，可用 `--data-type biological --also-thirdgen --also-secondgen`。`--mode both` 输出纯三代和融合，`--mode secondgen` 单独输出二代。`--data-type` 与 `--mode` 互斥。

纯二代的 Multi 类别来自共享 15 类二代轴，Binary 按三组专用配对选择。主表需提供中心参考碱基或短上下文，供对应家族筛选。

四套配置均可运行纯三代和融合，各可选 S1 或 S2。`--strategy s2` 是默认值，`s3` 是 S2 的旧名称。`--strategy both` 分别输出 S1 和 S2。纯二代只输出一次，策略目录为 `none/`。

<a id="preprocessing"></a>

## 2. 完成预处理

| 原始数据 | 上游产物 | NanoRemix 后续入口 |
| --- | --- | --- |
| RNA002 FAST5 | Step4 特征 TSV 与坐标侧表 | `prepare --chemistry RNA002 --format step4` |
| RNA004 POD5 | Step6 特征 CSV，含方向和坐标 | 显式正链旧 Step6 用 `prepare --chemistry RNA004 --format step6`。已正确处理方向的 prepared 表可直接用于预测。 |

NanoRemix 不从 FAST5/POD5 直接 basecall，也不直接提取电流特征。原始信号预处理与 NanoRemix 现代推理环境分开运行。下面的软件版本来自超算已有配置、环境元数据和日志；不能推断为每个历史作业都具有完整锁定环境。

### 2.1 RNA002 从 multi-read FAST5 到 9 nt 与 20 个信号值

1. 使用 **Guppy 6.3.8+d9e0f64**、`rna_r9.4.1_70bps_hac.cfg`、FLO-MIN106/SQK-RNA002 basecall，并保留 basecalled FAST5。
2. 用 **ont-fast5-api 4.1.1** 的 `multi_to_single_fast5` 拆成 single-read FAST5。已有配置环境是 Python 3.6.15。
3. 用 **Tombo 1.0** 执行 resquiggle，配合 **GraphMap 0.5.2** 和匹配的参考 cDNA；Tombo 是独立 Python 2.7 环境。检查 corrected Events/Alignment 和失败记录，不能只依据退出码。
4. 运行上游 `extract_tombo_step4.py`，提取连续 5 个事件的 5-mer 信号统计和同步 `.coordinates.tsv` sidecar。已有特征环境记录为 Python 3.6.9、NumPy 1.19.2、h5py 2.10.0。现存作业使用 `--clip 10`；脚本的 `--rna` 是提取器参数，不是 Tombo 参数。
5. 将 Step4 TSV 交给下方的 `prepare` 命令。每 5 个连续重叠的 5-mer 构成 9 nt，信号取中间行的 20 个值。完整保留邻接事件供窗口构造。
6. 纯三代可直接用 prepared。融合时需用 sidecar 的 `id` 显式连接 read、链、转录本 ID 和中心坐标，再从匹配的参考转录本或完整 corrected-event 序列得到同中心连续 81 nt。event index 不是转录本/基因组坐标；不能把剪切后的片段序列与旧 event index 混用。

RNA002 软件环境需拆开：Guppy/FAST5、Tombo、特征提取和 NanoRemix 推理依赖彼此较老/不同，不能把 Tombo 1.0 与现代模型环境强行装在同一个 Python 环境。

| 上游参数 | 含义 |
| --- | --- |
| Guppy `--fast5_out` | 保留带 basecalling 信息的 FAST5，供拆分及 Tombo 使用。 |
| Tombo `--basecall-group`、`--basecall-subgroups` | 定位 FAST5 内的 basecalling 数据。命令模板用 `Basecall_1D_000` 和 `BaseCalled_template`。 |
| 提取器 `--fast5-list`、`--output-prefix` | 指定单 read FAST5 清单和输出文件前缀。 |
| 提取器 `--rna`、`--clip 10` | 按历史直接 RNA 流程反转原始信号，并在 5-mer 窗口所需侧翼之外，从两端各裁去 10 个事件。 |

### 2.2 RNA004 从 POD5 到 17 nt 与 36 个信号值

1. 使用 **pod5 0.3.15**（现存环境 Python 3.9.20）读取 POD5。
2. 用 **Dorado 1.3.0** 和 `rna004_130bps_hac@v5.2.0` basecall，输出 move table（`--emit-moves`）。
3. 用相同参考的 Dorado aligner 比对，并用 samtools 过滤、排序和索引。现存脚本的 `samtools view -F 2308 -q 10` 排除 unmapped、secondary、supplementary alignment 并要求 MAPQ≥10；这是上游脚本设置，不是 NanoRemix CLI 隐式过滤。历史配置没有锁定 samtools 版本，新流程需记录 `samtools --version`。
4. 用 **ont-remora 3.3.0** 的 `SigMapRefiner` 和匹配的 `9mer_levels_v1.txt` 做 refined signal mapping/9-mer 特征，保留 NPZ 与逐行坐标 sidecar。现存特征环境记录 pysam 0.22.1、NumPy 1.24.4、pandas 2.3.3、h5py 3.11.0、torch 2.4.1。
5. 对齐核实 NPZ 与 sidecar 行数后由 `merge_step6.py` 汇总 Step6。中心基因组坐标为 `ref_start+4`。Remora 的数值字段映射为 mean、std、median、length；输出还应保留 read、reference、start/end、strand 和候选身份。
6. 构造中心正确、RNA 方向一致的 17 nt 上下文并保留中央 Step6 的 36 个信号数值。融合时另需相同候选的 81 nt 连续成熟 RNA 上下文；跨内含子应沿剪接路径提取。

**RNA004 Step6转换限制：**旧v0.3.0转换器按参考坐标升序拼接相邻9-mer，对RNA方向负链可能丢候选。v0.3.1的 `nanoremix prepare --format step6` 要求显式 `strand` 列，值为+/-；若包含任意负链行（包括正负链混合），在写出结果前明确报错。同一read若声明多个参考也会报错，避免跨比对拼接。不要改写链标签或删除负链来冒充完整处理；此保护没有实现通用负链适配。已有 `prepare_model_input.py` 和 splice-aware 上下文脚本依赖已知 catalog/truth 字段，也不能原样用于未知样本。新样本需按 `ref_start+4` 构造上下文、校验中央序列、同步反转序列与数值槽位（需要反向时），并记录每条拒绝原因；当前包没有声称包含经过验证的任意未知样本全转录组适配器。

RNA004 的参考序列、regions、BAM 和 FASTA 必须对应同一参考版本。9-mer region 文件格式和区间是候选 9-mer 左端起点，不等同于任意 site BED；细节见 [RNA004 pipeline 模板](guide_20261104/PIPELINE_命令与输入衔接.md)。

| 上游参数 | 含义 |
| --- | --- |
| Dorado `--emit-moves` | 将 move table 写入 BAM，后续提取器据此定位信号。 |
| samtools `-F 2308 -q 10` | 排除未比对、次要和补充比对，要求 MAPQ 至少为 10。 |
| 特征提取器 `--reg`、`--level` | 指定候选起点区间表和信号水平模型。区间表字段为 `seqnames,start,end,strand`，`[start,end)` 使用 0 基 9-mer 左端坐标。 |
| 特征提取器 `--kmer 9`、`--step 1`、`--resquiggle True` | 以 9 nt 为特征窗口，候选起点步长为 1，并开启信号映射细化。 |

### 2.3 检查上游产物并转换

| 输入格式 | 转换前核对 |
| --- | --- |
| RNA002 `step4` | TSV 的前七列依次为 `id`、`base`、`sequence_input`、`mean`、`std`、`median`、`length`，可附第八列 `kmer_raw`。序列是 5 nt，每组统计用竖线分隔 5 个值。保留用于连接的 `id` 和同步坐标侧表。[Step4 示例](guide_20261104/examples/RNA002_step4.tsv) |
| RNA002 `extracted` | 表头需含 `id`、`sequence_input`、`fast5_prefix`、`fast5_order`，以及 mean、median、std、length 各 5 个已展开的 `*_outputN` 列。 |
| RNA004 `step6` | 表头需含 `read_id`、`ref_start`、`sequence_input`、`strand`，以及 mean、median、std、length 各 9 个 `*_outputN` 列。序列是 9 nt，旧转换要求每行 `strand=+`。建议保留 `ref_name` 和 `ref_end` 便于核对。[正链 Step6 示例](guide_20261104/examples/RNA004_step6_positive.csv) |

在仓库根目录运行对应命令，并将输入路径换成自己的文件。

```bash
nanoremix prepare --input sample002.step4.feature.tsv --output prepared002.csv \
  --chemistry RNA002 --format step4

nanoremix prepare --input sample004.positive.step6.csv --output prepared004.csv \
  --chemistry RNA004 --format step6
```

上游命令模板包含每一步的读入路径、产物和核心参数，应按参考版本和本地工具环境填写。[FAST5/POD5 上游命令与坐标衔接](guide_20261104/PIPELINE_命令与输入衔接.md)

<a id="start-inputs"></a>
<a id="inputs"></a>

## 3. 整理输入文件

### 3.1 确认主表的文件格式和行身份

预测主表使用带表头的 UTF-8 CSV 或 TSV，支持 `.csv.gz` 和 `.tsv.gz`。一行对应一条 read 在一个候选中心的观测。跨表连接使用全表唯一、非空的 `record_key` 或由样本、read、位置组成的复合键。

### 3.2 特征输入

格式示例均为人工输入。[RNA002 prepared CSV](guide_20261104/examples/RNA002_prepared_with81.csv) · [RNA004 prepared CSV](guide_20261104/examples/RNA004_prepared_with81.csv)

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

常用元数据列如下。

| 列 | 含义 |
| --- | --- |
| `record_key` | 推荐提供全表唯一、非空的观测键。省略时工具按输入文件哈希和行号生成。跨表匹配时优先显式指定稳定键。 |
| `center_ref_base` | 未修饰参考碱基 A/C/G/U，T 会转换为 U。若同时提供短上下文，中心必须一致。 |
| `chemistry` | 可选，全表使用同一试剂版本。自动判断依据为此列或短上下文长度。 |
| `condition`、`data_type` | 可选。提供 `condition` 时需匹配 `--condition`。使用 `--data-type` 且表中有 `data_type` 时，两者必须一致。 |
| `sample_id`、`read_id`、`site_id`、`reference_name`、`chrom`、`pos`、`strand` | 可选，供定位、连接和分组。自行记录参考版本和坐标基准。默认输出保留这些列中输入时已存在的列。 |

### 3.3 三代预计算概率输入

`--input-kind probabilities` 可跳过三代模型特征推理，但不能跳过任务对应的概率格式检查。

| Task | 必需列 |
|---|---|
| Binary | 所选修饰与策略的 pair 概率，如 `s1_prob_m6A`、`s2_prob_m6A`；值表示专用 Binary pair 中的 P(mod)，不是 Multi 原始分数 |
| RNA002 Multi | 每种策略完整 18 类原生轴，如 `s2_prob_A` 加全部其他 RNA002 标签 |
| RNA004 Multi | 每种策略完整 16 类原生轴 |

所有概率须为 [0, 1] 内的有限值。Multi 完整轴每行概率和须在绝对误差 1e-6 内等于 1。`--strategy both` 时须提供两种策略的概率列。S2 也接受旧版 `s3_prob_*` 列，同一类别只能提交其中一种列名。输入还需 `9mer_final` 或 `center_ref_base`；只提供后者时，请明确指定 chemistry 或添加 `chemistry` 列。融合还需有效的二代来源。[概率表示例](guide_20261104/examples/README.md)

### 3.4 二代序列或概率

现场二代推理可在主表添加 `sequence81`，并传 `--sequence-column sequence81`；也可使用 `--context-fasta`，同时指定 `--fasta-id-column`、`--center-column`、`--coordinate-base 0|1`，需要时再指定 `--strand-column`。工具从指定 FASTA 记录的中心取前后各 40 nt。81 nt 必须是连续 A/C/G/T（U 转为 T）、中心第 41 位，不能含 N；方向须与候选一致。RNA002 的中央 9 nt 必须与 prepared 上下文一致；RNA004 的中央 17 nt 必须一致。跨剪接边界时需提供成熟转录本或沿正确 BAM CIGAR 路径提取的上下文；普通基因组连续窗口不会自动跳过内含子。

已有二代概率时可在主表提供完整的 `second_prob_<label>` 列，也可用 `--second-probabilities second.tsv --join-on record_key` 读取独立表。也接受原导出字段如 `probability_background_A` 或 `probability_PseudoU`，同一类别只提交一个字段。二代 15 类的概率必须有限且位于 [0, 1]，每行总和在绝对误差 1e-6 内等于 1。独立表键非空且唯一。与独立概率表同时指定 `--sequence-column` 时，工具校验的是独立表中的序列列，概率仍来自独立表。

使用序列做二代推理时，安装共享 SecondGen 模型。同一个融合命令可使用下列来源。

| 已有输入 | 替换 `--sequence-column sequence81` 的参数 | 还需SecondGen权重吗 |
|---|---|---|
| 独立二代概率表 `second.tsv` | `--second-probabilities second.tsv --join-on record_key` | 不需要 |
| 主表内嵌二代概率 | 提供全部 `second_prob_<label>` 列，不传序列/独立表参数 | 不需要 |
| 对应参考FASTA及正确中心坐标 | `--context-fasta transcripts.fa --fasta-id-column reference_name --center-column center_pos0 --coordinate-base 0`；按FASTA实际方向决定是否加 `--strand-column strand` | 需要 |

两种概率输入都必须有完整15类二代轴；独立表按唯一身份键匹配，不能只填待测修饰概率或按行号拼接。FASTA中心坐标属于所指定记录，工具不自动转换基因组/转录本坐标，也不自动从BAM解读剪接路径。已按转录方向的FASTA不要再次反向互补。详见[输入合同](#inputs)。

共享二代概率轴依次包含 `A,C,G,U,Am,m1A,m6A,Inosine,ac4C,m5C,Gm,m7G,m5U,pseudoU,Um`。默认扩展融合中，RNA002 的 hm5C、FormylC、o8G、m1psi 和 RNA004 的 hm5C、o8G、m1psi 只保留三代证据。共有模式使用两代都支持的类别。

<a id="start-install"></a>
<a id="models"></a>

## 4. 安装并检查模型

在 Python 3.10 或 3.11 环境中安装源码。以下命令安装两代的推理依赖。

```bash
git clone https://github.com/iris982022/nanoremix_tool.git
cd nanoremix_tool
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[models,second]' -c constraints-models.txt
nanoremix --version
```

也可通过 **Code → Download ZIP** 获取源码。解压后在源码根目录创建环境。Windows PowerShell 激活命令为 `.\.venv\Scripts\Activate.ps1`。

只从三代特征推理，或三代特征加已有二代概率时，可将安装命令中的 `'.[models,second]'` 改为 `'.[models]'`。三代和二代都用预计算概率时，运行 `python -m pip install -e .` 即可。

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

### 4.1 核对所下载的权重

| Release asset | 包含内容 | 用途 |
|---|---|---|
| `NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip` | 18 个 `.pkl` | RNA002 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip` | 6 个 `.pkl`、S1/S2 manifests | RNA002 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip` | 30 个 `.pkl` | RNA004 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip` | 10 个 `.pkl`、S1/S2 manifests | RNA004 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_SecondGen_20261104.zip` | 网络定义及共享 25 个 `.pt` checkpoint | 二代推理；RNA002/RNA004 共用 |

四组三代共 64 个 `.pkl`，共享二代共 25 个 `.pt`，完整集合为 89 个权重。Binary 消融模型未包含。模型包解压到项目根目录后组成 `models/`。包内源码没有模型权重，不需要重新训练。

单包预检时应按实际需要限定 `--chemistry`、`--task`、`--mods` 和 `--strategy`。若融合已有二代概率，把 `--data-type` 设为 `biological` 并加 `--second-source probabilities`，这样只检查所选三代模型。`doctor` 检查文件、清单、哈希和依赖，`validate` 运行实际输入校验与推理。两者的适用范围与报告含义见[模型说明](MODELS.md)。

<a id="start-run"></a>

## 5. 分路线运行预测

在已安装 NanoRemix 的仓库根目录运行。以下 `prepared002.csv`、`prepared004.csv` 是按第 2–3 步取得的输入表，`*_with81.csv` 额外包含同中心、同方向的 `sequence81`。每次使用新的输出目录。

先确定要运行的配置和命令。可把相同预测命令中的 `predict` 改为 `validate` 并省略 `--output`，对首批数据做一次实际推理检查。校验通过后再运行 `predict`。

<a id="rna002-binary"></a>

### 5.1 RNA002 Binary

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

### 5.2 RNA002 Multi

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

### 5.3 RNA004 Binary

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

### 5.4 RNA004 Multi

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

<a id="quickstart"></a>
<a id="examples"></a>

### 5.5 用人工概率快速检查安装

在已安装 NanoRemix 的仓库根目录运行，不需要下载模型：

```bash
nanoremix predict --input examples/minimal_binary_probabilities.csv \
  --output results/quickstart --chemistry RNA004 --task binary --mods m6A \
  --input-kind probabilities --mode thirdgen --strategy s2
```

预期两行结果分别为m6A（prob_m6A=0.90）和A（prob_A=0.80），输出在 `results/quickstart/thirdgen/s2/predictions.tsv`。人工概率用于软件检查，不代表模型准确率。

### 5.6 用预计算概率跑融合

此人工表同时提供 RNA004 S2 的二分类概率和完整二代概率。先做校验，再输出融合及其组成分支。两步均可在基础安装下运行。

```bash
nanoremix validate --input docs/guide_20261104/examples/RNA004_binary_both_probabilities.csv \
  --chemistry RNA004 --task binary --mods m6A,m5C,pseudoU \
  --input-kind probabilities --data-type biological --strategy s2

nanoremix predict --input docs/guide_20261104/examples/RNA004_binary_both_probabilities.csv \
  --chemistry RNA004 --task binary --mods m6A,m5C,pseudoU \
  --input-kind probabilities --data-type biological --strategy s2 \
  --also-thirdgen --also-secondgen --output results/probability_fused
```

<a id="outputs"></a>

## 6. 检查和分析结果

### 6.1 确认实际运行的路线和参数

先打开输出目录的 `run_manifest.json`。核对 `chemistry`、`task`、`output_modes`、`strategy` 和 `selected_labels_by_mode`，再看 `input_rows`、各路线的处理计数与模型哈希。使用预计算概率的路线可能没有模型哈希。

| 输出文件 | 内容 |
| --- | --- |
| `thirdgen/s2/predictions.tsv` | 纯三代预测及所选类别概率 |
| `fused/s2/predictions.tsv` | 融合预测、所选类别概率及证据来源 |
| `secondgen/none/predictions.tsv` | 纯二代预测，仅在请求对应路线时生成 |
| 同目录的 `selection_audit.tsv` | 全部候选的状态、候选标签和筛选原因 |
| `run_manifest.json` | 输入、模型哈希、参数、标签范围、统计和可选评价 |

### 6.2 先看审计，再读预测标签

1. 在 `selection_audit.tsv` 按 `status` 计数。`ok` 表示接受，`outside_selected_family` 表示参考中心不在选定类别范围内，`below_confidence` 等表示阈值筛除。搜索 `invalid_*`、`missing_*`、`insufficient_*` 和 `81nt_*`，按原因回查输入、身份键、坐标和参考方向。
2. 在 `predictions.tsv` 读 `predicted_label`、`prob_<label>`、`confidence`，融合时再看 `evidence_source`。`candidate_label` 可以保留被阈值拒绝的最高分候选，只有 `status=ok` 才写入接受的预测标签。
3. 需要对比三代和融合时，按 `record_key` 连表，并固定 chemistry、任务、策略和候选集合。纯二代可能覆盖更少的类别，比较前核对 `selected_labels_by_mode` 和审计状态。`--strategy both` 的文件分别位于 `s1/` 和 `s2/`。

以下命令统计第 5.6 节融合示例的状态。运行自己的数据后，将路径替换为对应路线的审计表。

```bash
python - <<'PY'
from collections import Counter
import csv

with open("results/probability_fused/fused/s2/selection_audit.tsv", encoding="utf-8", newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t")
    for status, count in Counter(row["status"] for row in rows).most_common():
        print(status, count)
PY
```

### 6.3 根据目的汇总或评价

- 筛选接受的记录可用 `--output-filter retained`，仅保存接受的修饰预测可用 `modified`。无论选何种过滤，`selection_audit.tsv` 仍保留全部审计记录。
- 输出以 read×候选为单位。若需位点汇总，先在原输入保留样本、参考、链和中心坐标，必要时用 `--keep-features` 保存自定义列。随后按位点和样本分组，明确覆盖度、有效观测数和拒绝观测的处理，再计算相应类别的比例。`confidence` 是所选类别内获胜类别的条件概率。
- 若有可靠真值，可用 `--truth-column` 输出混淆计数、保留行的条件准确率及把范围内拒绝行计作错误的准确率。不同路线应同时报告各自的保留数和评价范围。

<a id="options"></a>

## 7. 所有 CLI 参数

本节覆盖 `nanoremix` 根命令和全部五个子命令。每张表列出必填项、默认值、允许值及生效条件。开关参数直接填写选项名，省略时关闭。

### 7.1 根命令

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--version` | 开关 | 显示已安装的 NanoRemix 版本并退出。 |
| `--help`（`-h`） | 开关 | 显示用法与子命令。将它放在子命令后可查看该命令帮助，如 `nanoremix prepare --help`。 |

### 7.2 prepare 转换特征表

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--input` | 必填，文件路径 | 第 2 步中的旧格式特征表。Step4 用 TSV，extracted 或 Step6 用 CSV，也可用 `.tsv` 或 `.tab` 保存其制表符版本。 |
| `--output` | 必填，新文件路径 | 输出 prepared 表。`.tsv` 或 `.tab` 使用制表符，其余使用逗号。 |
| `--chemistry` | 必填，`RNA002` 或 `RNA004` | 决定上下文构造器和信号列数。 |
| `--format` | 必填，`step4`、`extracted` 或 `step6` | `step4` 和 `extracted` 只用于 RNA002。`step6` 只用于 RNA004，要求显式正链。 |
| `--help`（`-h`） | 开关 | 显示转换命令帮助。 |

### 7.3 catalog 查询类别

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--chemistry` | 默认显示两套，或指定 `RNA002`、`RNA004` | 限定输出的试剂版本。JSON 给出原生类别、融合共有类别、三代独有类别和二代类别。 |
| `--help`（`-h`） | 开关 | 显示类别查询帮助。 |

### 7.4 doctor 检查模型与依赖

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--model-root` | `models` | 待检查的模型目录。 |
| `--chemistry` | 默认检查两套，或指定 `RNA002`、`RNA004` | 限定试剂版本。 |
| `--task` | 默认检查两种，或指定 `binary`、`multi` | 限定任务。只解压一套模型包时请显式指定。 |
| `--mods` | `all`，或逗号分隔的修饰列表 | 限定选定配置的修饰。Binary 按配对检查权重，Multi 仍检查完整的组件。 |
| `--strategy` | `both`，可选 `s1`、`s2`、`both` | 检查一种或两种三代策略，旧名 `s3` 等同 `s2`。 |
| `--data-type` | `biological` 或 `synthetic`，默认 `biological` | `synthetic` 只检查三代。`biological` 根据二代来源决定是否同时检查共享二代权重。 |
| `--second-source` | `sequence` 或 `probabilities`，默认 `sequence` | 与 `biological` 一起使用。`sequence` 检查共享二代网络和权重，`probabilities` 跳过这部分。二代概率表在预测阶段验证。 |
| `--condition` | `full`、`mean`、`std`、`length`，默认 `full` | 对应全特征或单信号消融权重。消融需要显式指定 `--task binary` 并另外准备匹配模型。 |
| `--model-checksums` | 默认省略，接受 SHA256 索引 JSON | 额外核对所需文件的路径和摘要，如 `models/SHA256.json`。 |
| `--strict` | 开关，默认关闭 | 预检失败时退出状态为 2。省略时应从 JSON 的 `preflight_passed` 判断报告结果。参数或文件读取错误仍可能直接失败。 |
| `--help`（`-h`） | 开关 | 显示预检命令帮助。 |

<a id="predict-options"></a>

### 7.5 predict / validate 的 41 个长选项

两者使用相同的参数。`predict` 保存预测及审计文件，`validate` 执行相同的输入校验与推理，可省略 `--output`。`--mods` 与 `--targets` 必选其一，`--data-type` 与 `--mode` 也必选其一。[CLI 帮助全文](guide_20261104/CLI_help.txt)

#### 7.5.1 主表、类别与运行路线

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--input` | 必填，文件路径 | 主表，UTF-8 CSV 或 TSV，可 gzip 压缩。[输入格式](#inputs) |
| `--output` | `predict` 必填，`validate` 可省略 | 新输出目录。`validate` 接受该选项并检查路径，验证后仅返回摘要。 |
| `--task` | 必填，`binary` 或 `multi` | 专用修饰/未修饰配对，或参考碱基家族内的多类比较。 |
| `--chemistry` | `auto`、`RNA002`、`RNA004`，默认 `auto` | `auto` 根据统一的 `chemistry` 列或短上下文长度判断。主表只有参考中心碱基而无短上下文时，请提供试剂版本。 |
| `--mods` | `all` 或逗号分隔的修饰名 | 自动按 `--canonical` 加入未修饰对照。`all` 请求当前路线支持的全部修饰。 |
| `--targets` | 仅 Multi，逗号分隔的精确类别名 | 逐一列出参与竞争的类别，需要未修饰对照时显式写入。单个修饰类在同家族内须有竞争类。 |
| `--canonical` | `matching` 或 `all`，默认 `matching` | 使用 `--mods` 时，`matching` 加入目标修饰的同家族未修饰类，`all` 加入 A/C/G/U。 |
| `--data-type` | `biological` 或 `synthetic` | `biological` 运行融合，`synthetic` 运行纯三代。主表若带 `data_type`，需与所传值一致。 |
| `--mode` | `thirdgen`、`fused`、`both`、`secondgen` | 显式指定纯三代、融合、同时输出纯三代与融合、或纯二代。生物样本只做三代时用 `thirdgen`。 |
| `--also-thirdgen` | 开关，默认关闭 | 与 `--data-type biological` 一起使用，在融合结果外另存纯三代。 |
| `--also-secondgen` | 开关，默认关闭 | 与 `--data-type biological` 一起使用，在融合结果外另存二代可支持的比较。 |
| `--strategy` | `s1`、`s2`、`both`，默认 `s2` | 选择三代策略。`both` 分别保存 S1 和 S2，旧名 `s3` 等同 `s2`。纯二代策略目录为 `none/`。 |
| `--input-kind` | `features` 或 `probabilities`，默认 `features` | 选择三代 prepared 信号特征或预计算三代概率。 |
| `--condition` | `full`、`mean`、`std`、`length`，默认 `full` | 选模型的信号条件。单信号条件只用于 Binary，需匹配的消融权重；发布包提供 `full` 权重。 |

#### 7.5.2 模型路径和二代输入

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--model-root` | 默认省略，目录路径 | 对应模型目录的父目录，例如 `models`。在特征三代推理中必需。 |
| `--model-checksums` | 默认省略，SHA256 索引 JSON | 建议用 `models/SHA256.json` 核对实际加载的文件。加载资产时须提供 `--model-root`。 |
| `--trust-model-files` | 开关，默认关闭 | 加载模型 pickle、checkpoint 或模型 Python 文件时必须提供。用于用户已经核实来源的本地模型。 |
| `--second-model-root` | 未指定时使用 `<model-root>/secondgen` | 覆盖二代模型目录。加载并校验哈希时，该目录需位于 `--model-root` 内。 |
| `--second-probabilities` | 默认省略，独立 CSV/TSV 路径 | 提供完整 15 类二代概率，与 `--join-on` 一起使用。 |
| `--join-on` | 默认省略，逗号分隔的列名 | 独立二代概率表需要显式指定匹配键，如 `record_key` 或 `sample_id,read_id,site_id`。独立表键需非空且唯一。 |
| `--sequence-column` | 默认省略，列名 | 从主表的指定 81 nt 列运行二代推理。与 `--second-probabilities` 同用时，改为校验独立表中的 81 nt 列。 |
| `--context-fasta` | 默认省略，FASTA 路径，可 gzip | 从指定记录按坐标提取 81 nt，需同时提供下面三项必填坐标选项。 |
| `--fasta-id-column` | FASTA 时必填，主表列名 | 该列的值匹配 FASTA 标题的第一个字段。 |
| `--center-column` | FASTA 时必填，主表列名 | 指定 FASTA 记录中的候选中心坐标。 |
| `--coordinate-base` | FASTA 时必填，`0` 或 `1` | 坐标从 0 或 1 起始，须与坐标列的定义一致。 |
| `--strand-column` | FASTA 时可选，主表列名 | 值为 `+` 或 `-`。`-` 对提取窗口取反向互补，省略时直接使用 FASTA 方向。 |

#### 7.5.3 融合、筛选、评价与资源

| 参数 | 默认值或允许值 | 含义 |
| --- | --- | --- |
| `--fusion-scope` | `extended` 或 `common`，默认 `extended` | Multi 融合中，`extended` 保留三代独有类别，`common` 仅选择两代共有类别。与纯三代一起运行时，共有模式的类别选择也适用于该次运行。 |
| `--second-weight` | `1.0`，有限正数 | 二代在融合公式中的权重。 |
| `--third-weight` | `1.0`，有限正数 | 三代在融合公式中的权重。两项权重可独立设置。 |
| `--binary-bias` | `0.0`，有限实数 | Binary 融合 log-odds 相加时的固定偏置。只有 Binary 融合接受非零值。 |
| `--min-confidence` | `0.0`，[0, 1] | 获胜类别的条件概率下限。 |
| `--min-selected-mass` | `0.0`，[0, 1] | 所选类别在参考碱基家族中的概率质量占比下限。纯二代 Multi 使用二代家族，三代或融合 Multi 使用三代家族。 |
| `--min-second-mass` | `0.0`，[0, 1] | 在有二代结果的路线中，二代所选类别的家族概率质量占比下限。 |
| `--output-filter` | `all`、`retained`、`modified`，默认 `all` | 分别保存全部候选、接受候选、接受且判为修饰的候选。审计表始终写出全部候选。 |
| `--on-invalid` | `error` 或 `mark`，默认 `error` | 支持的行级问题触发中止或留下明确状态。缺列、重复身份键和模型错误会中止。 |
| `--truth-column` | 默认省略，列名 | 指定独立的真值列，仅用于在运行清单中评价。 |
| `--chunk-size` | `50000`，正整数 | 主表每批读入的行数。身份键集合、FASTA 和独立二代概率表仍可能驻内存。 |
| `--batch-size` | `256`，正整数 | 二代序列网络的批大小。 |
| `--device` | `cpu`、`cuda`、`mps`，默认 `cpu` | 二代序列网络的运行设备，须与本地 PyTorch 环境相容。 |
| `--keep-features` | 开关，默认关闭 | 写出额外的输入列，跳过输入概率列和下划线开头的内部列。 |
| `--help`（`-h`） | 开关 | 显示相应子命令的帮助。 |

三项筛选阈值只在对应指标严格低于设定值时拒绝该候选。默认阈值为 0 时仍执行格式、方向和概率支持检查。融合权重必须大于 0，运行某个模型分支时直接用相应 `--mode`。[计算方法](#calculation)

<a id="calculation"></a>

## 8. 计算方法与参考资料

### 8.1 Binary 三代与融合

S1 的多个校准 pair 分数按 logit 相加：

```text
p3 = sigmoid(sum_j(logit(p3_j)))
```

S2 使用单个联合模型分数。二代完整 15 类概率先按目标正常碱基/修饰 pair 条件化，再按默认融合权重和偏置合并：

```text
pf = sigmoid(w2*logit(p2) + w3*logit(p3) + beta)
```

默认 `w2=w3=1, beta=0`；例如 p2=0.8、p3=0.6，结果为 0.857142857。logit 概率会在计算时进行数值边界保护。

### 8.2 Multi 三代内部聚合与两代融合

三代 S1/S2 的多个组件概率按模型 manifest 中的 class axis、家族权重和 temperature 聚合。随包实际 manifest 使用冻结的 geometric/gate 参数；不要把它改写成简单等权平均。具体值见 `guide_20261104/frozen_multi_fusion.json`。

融合时先在研究标签和参考碱基家族内条件化概率。对二代与三代共有类别，按默认权重 1/1 做归一化概率幂乘积；`extended` 将融合结果乘回共有类别原来占有的三代质量，三代独有类别保持原三代概率。完全没有有效共有正概率支持时会标为无效，不造均匀分布。`common` 先限制到共有轴再融合。二代不支持的类别必须通过 `evidence_source` 识别为三代证据。

### 8.3 阈值与决策顺序

先检查 schema、chemistry、中心/方向、数值和身份，再确定类别/家族并计算概率；接着计算 selected mass/second mass、argmax 与 confidence；只有 status 为 `ok` 才有接受的 `predicted_label`。阈值拒绝通常保留 `candidate_label` 和审计原因，但 accepted label 为空。精确公式和数值案例见 [详细计算章节](guide_20261104/README_版本输入预处理参数与计算.md)。

<a id="references"></a>

### 8.4 参考资料

- [详细中文方法手册](guide_20261104/README_版本输入预处理参数与计算.md)：全部 41 个 predict/validate 选项、字段、版本证据、完整筛选顺序和公式例子。
- [RNA002/RNA004 上游命令模板](guide_20261104/PIPELINE_命令与输入衔接.md)：命令、Step4/Step6 输入衔接和坐标示例。
- [输入列定义](INPUTS.md)、[模型布局与依赖](MODELS.md)、[CLI 选项说明](OPTIONS.md)、[限制说明](KNOWN_LIMITATIONS.zh-CN.md)。
- [人工输入示例目录](guide_20261104/examples/README.md)：所有示例均为人工数据，不是准确率 benchmark。
- [原始上游副本索引](../upstream_reference/README.md) 与 [`model_path_mapping.tsv`](guide_20261104/model_path_mapping.tsv)：上游脚本/权重路径对应。

源码包可在仓库根目录用 `python verification/verify_distribution.py` 校验文件哈希与清单。源码许可与权重来源见 [LICENSE](../LICENSE) 和 [NOTICE.md](../NOTICE.md)。
