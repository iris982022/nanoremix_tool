# NanoRemix v0.3.0：版本、输入、预处理、参数与计算说明

> **2026-11-05维护版说明：**下文v0.3.0描述保留原始模型和历史流程依据。v0.3.1不改变模型预测与权重，修复分发并保护旧Step6入口：RNA004必须显式声明strand（+/-），包含负链、缺失/无效方向或同read多个参考会报错且不写输出。历史负链样例用于说明限制，不代表v0.3.1已经提供通用负链适配。


核对日期：2026-11-04。方法核对所依据的原始 v0.3.0 基线包名为：

```text
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/NanoRemix_complete_v0.3.0.zip
```

该原 ZIP 是本次分析时的可安装源码和权重基线；发布权重已拆分成更小的 Release ZIP，原 ZIP 不再作为当前上传文件。本文以原包内的 `src/nanoremix/` 实现为准。上游预处理则依据超算现存脚本、软件环境和日志说明。**没有重新运行原始数据预处理或证明跨流程生物学等效性。** 附带的人工示例用于格式和计算核验，不可用于评价准确率。核验结果、完整 CLI help、精确类别清单、冻结的 Multi 聚合参数、模型路径对应表见同目录附件。

## 1. 应该使用哪个包

| 文件 | 包含内容 | 是否有模型 |
|---|---|---|
| `NanoRemix_complete_v0.3.0.zip`，404,858,144 bytes | 可安装的统一命令 `nanoremix`、源码、文档、测试及推理模型 | **有：64 个三代 `.pkl`，25 个二代 `.pt`** |
| `NanoRemix_final_code_bundle.zip`，46,491,918 bytes | 原始四套工具的源码，以及论文绘图、验证、benchmark 与来源审计代码 | **没有模型权重** |

原始完整包的89个权重已逐个与当时正式 `Binary/`、`Multi/` 目录中的对应文件进行 SHA256 对照，全部一致。`complete` 包中的 `NOTICE.md` 也准确引用了 `final_code_bundle.zip` 的 SHA256。两者当时分别是统一推理包与独立源码分析归档；当前小型上传集合不包含论文分析归档，它不是运行依赖。

完整包是2026-10-09的快照；当前 `FINAL_TOOL/` 下11月更新的转换器、运行脚本、文档没有自动进入该 ZIP。尤其不要将旧四包入口的参数直接传给 v0.3.0。

本文“二代”指81 nt序列模型，“三代”指Nanopore序列与信号特征模型。这两个名称在这里用于区分模型分支；二代分支不直接读取Illumina FASTQ/BAM。

## 2. 有几个版本和运行组合

### 2.1 四套任务、八种主要工作流

chemistry × task 一共有4套：RNA002 Binary、RNA002 Multi、RNA004 Binary、RNA004 Multi。每套都可以运行纯三代或二代＋三代融合，所以有**8种主要工作流**。

| chemistry | task | 输出模式 | 必需主输入 | 额外二代输入 | 所用权重 |
|---|---|---|---|---|---|
| RNA002 | Binary | 纯三代 | RNA002 prepared：9 nt＋20数值列 | 不需要 | RNA002 Binary |
| RNA002 | Binary | 融合 | 同上 | 81 nt序列/FASTA，或完整15类二代概率 | RNA002 Binary＋共享二代模型；已有二代概率时不加载二代权重 |
| RNA002 | Multi | 纯三代 | RNA002 prepared：9 nt＋20数值列 | 不需要 | RNA002 Multi |
| RNA002 | Multi | 融合 | 同上 | 同上 | RNA002 Multi＋共享二代模型，或二代概率 |
| RNA004 | Binary | 纯三代 | RNA004 prepared：17 nt＋36数值列 | 不需要 | RNA004 Binary |
| RNA004 | Binary | 融合 | 同上 | 同上 | RNA004 Binary＋共享二代模型，或二代概率 |
| RNA004 | Multi | 纯三代 | RNA004 prepared：17 nt＋36数值列 | 不需要 | RNA004 Multi |
| RNA004 | Multi | 融合 | 同上 | 同上 | RNA004 Multi＋共享二代模型，或二代概率 |

以上每种含 `s1`、`s2` 两套三代策略，因此可细分为16种策略配置。`--strategy both` 一次输出两套结果，**不额外产生一套共识模型**。四套任务共享同一组25个二代网络，不是四组独立重训的二代权重。

此外可以单独输出纯二代。`--mode secondgen --task multi` 的可选类别来自共享15类轴，无需三代数值特征和三代权重；仍需明确参考中心碱基、身份键、序列/二代概率来源。

### 2.2 S1/S2与旧工具名称的对应

| 论文/本 v0.3.0 CLI | 旧四包代码名称 | 含义 |
|---|---|---|
| `s1` | `s1` / `strategy1` | 独立位移视图；RNA002为5个组件，RNA004为9个组件 |
| `s2`（默认） | `s3` / `strategy3` | 联合多位置特征的单个 `m_all` 组件 |
| `both` | 此处只表示分别运行上述两套 | 保存独立的 `s1/`、`s2/` 结果 |

本 CLI 接受 `--strategy s3` 作为 `s2` 的旧名，但不是另一套模型。Multi `models/.../multi/s2/manifest.json` 内仍写 `strategy: "s3"`，这是保留的原始模型标识，不能手工改写。

**旧四包 Binary 中 `both` 可组合S1和S3概率；本统一 v0.3.0 的 `--strategy both` 不做该组合。** 相同权重不意味着两种入口在所有配置下的最终输出相同。

### 2.3 支持修饰

Binary两种chemistry均仅支持3套专用二分类：`A/m6A`、`C/m5C`、`U/pseudoU`。

| 家族 | RNA002 Multi三代 | RNA004 Multi三代 |
|---|---|---|
| A | A、m6A、m1A、Inosine、Am | A、m6A、m1A、Am |
| C | C、m5C、hm5C、ac4C、FormylC | C、m5C、hm5C、ac4C |
| G | G、m7G、o8G | G、o8G、Gm |
| U | U、pseudoU、m1psi、m5U、Um | U、pseudoU、m1psi、m5U、Um |
| 总计 | 18类，含4个正常碱基 | 16类，含4个正常碱基 |

共享二代15类轴，**按下列顺序**：

```text
A,C,G,U,Am,m1A,m6A,Inosine,ac4C,m5C,Gm,m7G,m5U,pseudoU,Um
```

融合的 `--fusion-scope common` 仅支持两代共有类别：RNA002为14类，RNA004为13类，均包含ACGU。默认 `extended` 可保留完整三代18/16类：RNA002的 `hm5C,FormylC,o8G,m1psi`，RNA004的 `hm5C,o8G,m1psi` 沿用三代证据。不能将这些类别描述为获得了独立二代证据。

`Cm` 不在上述正式模型范围。RNA002三代不含Gm；RNA004三代不含Inosine、m7G。别名如 `RNA-editing→Inosine`、`PseudoU→pseudoU`、`T→U` 可被规范化，但不会扩大模型范围。

## 3. 各模式到底需要什么输入文件

### 3.1 主输入 prepared 表

`nanoremix predict --input-kind features`（默认）接收UTF-8 CSV或TSV，可用 `.csv.gz` / `.tsv.gz`。每行是一个read覆盖一个候选中心的记录。不能先把各read的信号取平均再当成一个site输入。

| 字段 | RNA002 | RNA004 |
|---|---|---|
| `9mer_final` | 9个碱基token，中心第5位 | 17个token，中心第9位 |
| `mean_output1..K` | K=5 | K=9 |
| `median_output1..K` | K=5 | K=9 |
| `std_output1..K` | K=5 | K=9 |
| `length_output1..K` | K=5 | K=9 |
| 推荐身份字段 | `record_key,sample_id,read_id,site_id` | 同左 |
| 可选中心碱基 | `center_ref_base`，须与上下文中心一致 | 同左 |
| 可选参考坐标 | 同批坐标表中的转录本ID和中心0基位置 | `reference_name,center_pos0,strand` 等 |

数值列建议按位置写为 `mean_output1,median_output1,std_output1,length_output1,...`。工具按列名对齐，不要求文件中恰好按这一顺序。所有值必须有限；std与length非负。length是电流采样点数/dwell统计，不是序列长度。不要把mean电流、标准差或dwell任意标准化后再送原权重。

推荐用正常参考碱基，如：

```text
RNA002: A_C_G_U_A_C_G_U_A                 # 9 token
RNA004: A_C_G_U_A_C_G_U_A_C_G_U_A_C_G_U_A # 17 token
```

字段名 `9mer_final` 沿用历史名称，RNA004仍叫这个名字。内部使用的 `9mer` 中间列可能含81个拼接token；这**不等于连续81 nt序列**。

`record_key` 须全文件非空且唯一。未提供时，预测阶段生成“源文件SHA256＋0基输入行号”的键。一个site有多条read时，不应直接用不唯一的site_id作为record_key。跨文件连接最好事先提供稳定的身份键。真值只需在评价时额外提供，不参与预测。

同一输入不要混RNA002/RNA004，也不要混多个样本中可能重名的read。`--chemistry auto` 只根据chemistry列或上下文长度判断，不能根据文件名证明实验化学版本。

### 3.2 纯三代

只需相应prepared表和三代权重，不需要81 nt、二代概率或二代权重。

```bash
# 生物样本也可以明确选择纯三代；不要为了纯三代把真实样本标成synthetic。
nanoremix predict --input sample002.prepared.csv --output results/002_multi_third \
  --chemistry RNA002 --task multi --mods all --mode thirdgen --strategy s2 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

`--data-type synthetic` 是“合成数据、仅三代”的常规路由；生物样本的纯三代可用上面的显式 `--mode thirdgen`。有 `data_type` 列时，用 `--data-type` 会检查每行与声明一致。`--data-type` 与 `--mode` 互斥。

### 3.3 融合：prepared之外提供以下一种来源

| 二代来源 | 具体输入与参数 | 是否加载25个二代网络 |
|---|---|---|
| 表内连续81 nt | prepared增加 `sequence81`（或自定义列名），用 `--sequence-column sequence81` | 是 |
| FASTA＋明确中心坐标 | `--context-fasta`，同时提供 `--fasta-id-column`、`--center-column`、`--coordinate-base`；必要时加 `--strand-column` | 是 |
| 表内完整15类概率 | 每类一列 `second_prob_<label>`，例如 `second_prob_A`、`second_prob_m6A` | 否 |
| 独立完整15类概率表 | `--second-probabilities second.tsv --join-on record_key`，或明确的复合键 | 否 |

后两种都是预计算二代概率。每行须是完整15类分布，有限、非负且总和与1相差不超过1e-6。不能只给所选的两三类，也不能把Multi三代概率冒充二代概率。

原导出器列名也可识别，如 `probability_background_A`、`probability_background_T`、`probability_RNA-editing`、`probability_PseudoU` 等。每个类只允许一种列名；标准列与别名同时出现会因歧义报错。

独立二代表的连接键必须非空且唯一。主表允许多条read连接同一个明确的序列/site概率，但键须覆盖所有必要维度，例如 `sample_id,site_id,strand`；工具不会替用户判断一个site是否跨样本重名。禁止按行号拼接。若独立概率表同时指定 `--sequence-column`，该列用于核对独立表中的序列，**不会重新运行二代模型**。

81 nt必须满足：

1. 连续正常碱基，仅ACGT；输入U会转成T，不能含N。
2. 中心为0基40，即第41位；与参考中心碱基一致。
3. RNA002中央9 nt为 `seq[36:45]`，须与prepared上下文一致。
4. RNA004中央17 nt为 `seq[32:49]`，须与prepared上下文一致。
5. 序列、数值槽位、参考中心必须在同一方向。FASTA已按转录方向排列时，不应再因原始基因组负链而反向互补一次。

FASTA坐标是**所命名FASTA记录内的中心坐标**。例如 `fast5_order` 只有在完整事件序列FASTA且编号仍匹配时才是正确中心；不能拿它直接定位基因组/转录组FASTA。工具取 `[center-40:center+41]`，不足40 nt侧翼时不填N。

### 3.4 预计算三代概率替代特征

用 `--input-kind probabilities`。这一路不加载三代模型，也不要求20/36个数值列。仍需 `center_ref_base`，并明确 `--chemistry` 或一致的chemistry列。

| task | 必需三代概率列 |
|---|---|
| Binary | 每个所选修饰、所选策略的pair概率，如 `s2_prob_m6A`；表示专用A/m6A二分类中的P(m6A)，不是Multi的原始m6A分数 |
| Multi | 每个所选策略的**完整原生类别轴**：RNA002 18列，RNA004 16列，如 `s2_prob_A,...`；每行总和为1 |

`--strategy both` 要同时提供s1和s2。旧 `s3_prob_*` 可以代替s2，但同类同时给s2和s3会报歧义。融合时仍需第3.3节的二代来源。

### 3.5 可直接参考的现存真实prepared文件

以下文件已存在且表头与所述格式一致；它们是经过catalog筛选的验证数据，不代表无筛选全转录组输入：

```text
RNA002:
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA002/HEK293_WT/catalog_prepared.csv
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA002/HEK293_WT/binary_catalog_prepared.csv

RNA004:
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/HeLa/fused_model_results/input/multi_human_fused_input.csv.gz
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/HEK293T/fused_model_results/input/multi_human_fused_input.csv.gz
```

这些表的81 nt列名为 `reference_sequence_81bp`，所以统一工具需要显式 `--sequence-column reference_sequence_81bp`，不会因为该列名存在而自动选择它。Binary仅对所选A/C/U家族使用专用pair，表内其他修饰真值不会自动变成正确的二分类benchmark标签。

## 4. 原始数据怎样生成输入：RNA002

### 4.1 已有工具和版本

| 阶段 | 软件/版本 | 版本依据 |
|---|---|---|
| RNA002 basecall | Guppy **6.3.8+d9e0f64**；`rna_r9.4.1_70bps_hac.cfg`；FLO-MIN106 / SQK-RNA002 | HEK293_WT的Guppy日志直接记录 |
| multi→single FAST5 | ont-fast5-api **4.1.1**；现环境Python3.6.15 | 2026-11-04读取配置所指环境的包元数据 |
| 事件校正与比对 | **Tombo 1.0**，环境Python2.7.15；**GraphMap v0.5.2** | Tombo包元数据，GraphMap二进制help |
| 5-mer特征 | 本地 `extract_tombo_step4.py`；现环境Python3.6.9、NumPy1.19.2、h5py2.10.0 | 脚本及配置所指环境 |
| prepared转换 | 本v0.3.0的 `nanoremix prepare` | 包内代码 |
| 模型推理 | Python3.10/3.11、sklearn1.3.0、XGBoost2.1.4 | 包依赖及原权重兼容要求 |

环境元数据表示当前安装版本，不将其冒充历史每个作业的完整软件锁定记录。Guppy的6.3.8有实际运行日志支持。各阶段使用独立环境；不要将Tombo1.0与现代推理依赖装进同一个Python环境。

### 4.2 已执行流程与产物

```text
RNA002 原始 multi-FAST5
  → Guppy RNA002 basecall，保留FAST5
  → multi_to_single_fast5
  → Tombo resquiggle＋GraphMap，对参考cDNA比对并写corrected Events
  → extract_tombo_step4.py：5-event统计＋同步坐标
  → Step4 TSV
  → nanoremix prepare --format step4
  → prepared CSV（9 nt、20数值列）
  ├→ 纯三代 Binary/Multi
  └→ 连接正确的参考/事件81 nt → 二代＋三代融合
```

现存脚本目录：

```text
/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/
```

| 顺序 | 现存脚本 | 关键内容/产物 |
|---|---|---|
| 1 | `slurm/01_basecall_rna002.sbatch` | `--flowcell FLO-MIN106 --kit SQK-RNA002 --fast5_out`；生成 `01_guppy_rna002/workspace/` |
| 2 | `slurm/02_make_single_fast5.sbatch` | `multi_to_single_fast5 --recursive`；生成 `02_single_fast5/` |
| 3 | `slurm/03_tombo_resquiggle.sbatch` | `tombo resquiggle`，明确basecall组，GraphMap、`--include-event-stdev`；抽样检查成功Events/Alignment |
| 4 | `slurm/04_extract_step4.sbatch` | 调 `python/extract_tombo_step4.py --rna --clip 10`；生成 `.feature.tsv` 和 `.coordinates.tsv` |
| 5 | 本发布包 `nanoremix prepare` | 生成通用prepared，不做真值/catalog筛选 |
| 6 | 本发布包 `nanoremix predict` | 指定Binary/Multi、策略、模式和研究修饰 |

当前参考转录组为 `/gpfs/work/bio/yuxinzhang17/DirectMultiRM/Homo_sapiens.GRCh38.cdna.all.fa`。换样本时应使用匹配物种、参考版本和转录本ID的FASTA。

Tombo会写corrected组，因此现存脚本处理新输出目录中的single-FAST5副本。`--rna`在这里是**特征提取脚本**的参数；不要把它不加核对地写成Tombo1.0命令参数。实际 `.sbatch` 中的Tombo命令未传 `--rna`。

旧 `submit_pipeline.sh` 会继续运行历史m6A-only推理、DRACH聚合与绘图。若要本统一v0.3.0的通用Binary/Multi，取其Step1–4作为上游，随后转到本工具，不应把旧05–07步骤当作新通用入口。提交资源、账号、数组与依赖按现存脚本和自己的样本配置；本文没有提交Slurm作业。

### 4.3 Step4与prepared是不同层级

Step4是TSV，前7/8列为：

```text
id  base  sequence_input  mean  std  median  length  [kmer_raw]
```

其中 `sequence_input` 长5，每个统计字段由5个 `|` 分隔的值构成。id形如 `read_id|event_index|strand`。一个Step4行对应中心事件及其周围5个事件的统计。

转换示例，在解压后的NanoRemix根目录运行：

```bash
nanoremix prepare --input /path/to/sample.step4.feature.tsv \
  --format step4 --chemistry RNA002 --output /path/to/sample002.prepared.csv

# 若已经是展开的5-mer extracted表：
nanoremix prepare --input /path/to/extracted002.csv \
  --format extracted --chemistry RNA002 --output /path/to/sample002.from_extracted.csv
```

转换器按read分组、按event index排序，要求连续5行编号每次加1，且相邻5-mer的4 nt重叠一致。这5行形成9 nt上下文；取**中间那行**的20个数值，不对5行信号再平均。无有效窗口时会报错。

不要先按修饰、中心碱基或truth删掉Step4邻居，否则会破坏候选构造。raw id应真实携带event index；旧读入器对缺失编号有按行补号逻辑，不能将补号视为已经验证的坐标。

`prepare` 是全表读入和按read构造，`predict --chunk-size` 不能控制其内存。亿级Step4应在上游按完整read分片，每片保留完整相邻行；不能任意每N行切开后忽略边界候选。保持样本分开，避免重复read ID。

### 4.4 信号数值如何计算

RNA raw signal按原提取定义先反转。设反转后的ADC数组为x，U为其**唯一值集合**：

\[
c=\operatorname{median}(U),\quad d=\operatorname{median}(|U-c|)/0.6744897501960817,
\qquad x'_i=(x_i-c)/d.
\]

以Tombo corrected事件边界切分x′，每个事件段S输出：

\[
\mathrm{mean}=\frac1n\sum S_i,\quad
\mathrm{std}=\sqrt{\frac1n\sum(S_i-\overline S)^2},\quad
\mathrm{median}=\operatorname{median}(S),\quad \mathrm{length}=n.
\]

std使用population SD（`np.std` 默认ddof=0）。每行取连续5事件。每端clip10事件，另留2事件构成5-mer。MAD为0、事件越界、异常碱基、失败的Tombo记录会被拒绝或记录错误。

### 4.5 RNA002融合所需81 nt如何得到

两条可用语义：

1. **参考转录组窗口**：用同批 `.coordinates.tsv` 中的 `id` 与prepared中间行的 `id` 显式连接。获得转录本ID、中心 `transcript_position0` 和链。将转录本ID列传给 `--fasta-id-column`，该0基中心传给 `--center-column --coordinate-base 0`；若需要相对FASTA反向互补则给正确链列。提取的81 nt必须通过中央9 nt检查。sidecar读入和连接时核实唯一性，不按行号。
2. **完整corrected事件序列FASTA**：从每个single-FAST5的完整 `RawGenomeCorrected_000/BaseCalled_template/Events.base` 按事件顺序导出，header精确匹配fast5_prefix，中心为原始0基fast5_order。不能只导出clip后的序列却仍用旧下标。现存Step4提取器仅输出feature/coordinate/error/summary，**不自动导出该完整FASTA**；不能在流程图中虚构已有 `.feature.fa` 产物。

第一条用参考序列，第二条用校正事件序列；应记录来源，不能把事件上下文叫作额外进行过的二代测序。若二代已有相同身份的完整15类概率，可直接走概率连接路线。

## 5. 原始数据怎样生成输入：RNA004

### 5.1 已有工具和版本

| 阶段 | 工具/配置 | 依据与边界 |
|---|---|---|
| POD5读取 | pod5 **0.3.15**，Python3.9.20 | 当前remora-env元数据 |
| Basecall | **Dorado1.3.0**，模型 **`rna004_130bps_hac@v5.2.0`** | 现存HeLa/HEK293T配置；需 `--emit-moves` |
| 比对 | 同套Dorado aligner＋samtools | 脚本用module加载samtools，未锁定其历史版本；实际运行应另存 `samtools --version`，本文不猜版本 |
| 信号对齐和9-mer统计 | **ont-remora3.3.0**，`SigMapRefiner`，`9mer_levels_v1.txt` | 当前remora-env元数据及提取代码 |
| BAM/FASTA坐标 | pysam **0.22.1** | 当前remora-env元数据 |
| 特征环境 | NumPy1.24.4、pandas2.3.3、h5py3.11.0、torch2.4.1 | 当前安装记录；与推理环境分开 |
| 推理 | 与RNA002相同的独立模型兼容环境 | 第6节 |

电流水平表的具体文件是：

```text
/gpfs/work/bio/yuxinzhang17/DirectMultiRM/ONT/9mer_levels_v1.txt
```

换新预处理实现、basecaller模型或level表时，需重新核对输入分布及独立数据表现。RNA004特征裁成5个槽位也不能当作RNA002正式输入。

### 5.2 已执行的流程

```text
RNA004 POD5
  → Dorado1.3.0 / RNA004 v5.2.0 basecall，输出含move table的BAM
  → Dorado aligner＋samtools过滤、排序、索引
  → Remora refined signal mapping，按预先指定region提取9-mer统计
  → NPZ＋逐行坐标CSV
  → merge_step6.py → Step6 data_df.csv
  → 构造同中心、同方向的17 nt参考上下文和36个信号值 → prepared
  ├→ 纯三代
  └→ 准备同中心81 nt（跨剪接用BAM路径） → 融合
```

脚本目录：

```text
HeLa:
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/
HEK293T:
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_HEK293T_pipeline/
通用9-mer提取器:
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_pipeline/feature_extraction_rna004.py
```

| 顺序 | 脚本 | 实际参数/产物 |
|---|---|---|
| 0 | `00_extract.sbatch` | 如原始为tar归档，解包并按split分片；归档不是模型输入 |
| 1 | `01_basecall.sbatch` | `basecaller --emit-moves --recursive --batchsize 512 --chunksize 2000`，给定RNA004模型；输出split BAM |
| 2 | `02_align.sbatch` | `dorado aligner REF BAM`；`samtools view -bh -F 2308 -q 10`；sort/index |
| 3 | `03_features.sbatch` | Remora提取器：`--kmer 9 --step 1 --resquiggle True`，指定POD5/BAM/level/regions；输出 `splitN.npz` 和 `splitN.csv` |
| 4 | `04_step6.sbatch` / `merge_step6.py` | 核对NPZ与sidecar行数，保留read/reference/start/end/strand；生成Step6 CSV |
| 5 | `05_prepare_model.sbatch` / `prepare_model_input.py` | **现有benchmark专用**：连接reserved catalog、检查参考方向、构造17 nt、筛选catalog及truth |
| 6 | `06_fused_context.sbatch` | BAM CIGAR路径提取剪接感知81 nt；当前实现仍依赖catalog字段 |
| 7 | 本v0.3.0 `predict` | 使用对应prepared和可选81 nt列 |

`-F 2308` 排除unmapped(4)、secondary(256)、supplementary(2048)，并保留MAPQ≥10。此为**上游脚本过滤规则**，不是NanoRemix CLI隐藏阈值。是否过滤duplicate、低质量read等，应在自己的分析设计中明确，不能声称已有CLI自动处理。

现存人类参考为 `/gpfs/work/bio/yuxinzhang17/DirectMultiRM/hg38.fa`。BAM、region、catalog/自定义候选和FASTA必须来自同一参考版本。

### 5.3 region文件不是任意BED的同义词

提取器需要带表头CSV：

```csv
seqnames,start,end,strand
chr1,9992,10001,+
```

代码把 `[start,end)` 中每个整数当成一个**9-mer左端起点**，不是把该区间当作一个待测site本身。对于0基目标中心p，要覆盖9个相邻9-mer起点以重建17 nt，可给起点范围 `[p-8,p+1)`；其中央9-mer起点为p-4，最终中心为p。必须保证起点非负、完整9-mer在参考内。

若只提取一条中央9-mer，可用外部可信参考构造17 nt；但不能期待旧相邻9行转换器从一行Step6重建其余侧翼。跨剪接边界需使用相应参考路径，不能把包含内含子的连续窗口视为成熟RNA序列。

当前提取器会收集指定region涉及的read并处理相应POD5；这是一条有region约束的链。`--step 1` 是此脚本的候选起点采样设置，稀疏采样可能使旧连续窗口构造失败。

### 5.4 Step6数值与方向

Step6至少包含：

```text
read_id,ref_start,sequence_input,
mean_output1..9,std_output1..9,median_output1..9,length_output1..9
```

实际合并器还保留 `ref_name,ref_end,kmer_length,strand`，这些字段在后续定位中心和方向时应保留。`sequence_input` 是9 nt，Step6还不是17 nt prepared。

Remora在refined mapping上调用 `dwell_mean_sd`，另计算quantiles。NPZ stat每位置为：

```text
0 mean, 1 sd, 2 dwell, 3 min, 4 q1, 5 median, 6 q3, 7 max
```

`merge_step6.py` 对应取0→mean、1→std、5→median、2→length。正常预测只使用每位置的这4项，level和basecall-error等NPZ字段不是额外直接传入CLI的特征。

SigMapRefiner参数为 `do_rough_rescale=True, scale_iters=0, do_fix_guage=True`。这里的mean/SD在Remora归一化信号尺度上计算；不能套用RNA002 unique-ADC MAD公式再归一遍。

原提取器将负链序列和所有数值位置保持在RNA方向。将窗口反向互补时，数值槽位也必须同步反序，不能只改序列。对当前Step6，中心基因组坐标为 **`ref_start+4`**，`ref_start`本身是中央9-mer的左边界。

### 5.5 v0.3.0旧prepare转换器的实际限制

命令确实存在：

```bash
nanoremix prepare --input /path/to/data_df.csv \
  --format step6 --chemistry RNA004 --output /path/to/sample004.prepared.csv
```

但它按read分组、按 `ref_start` **升序**找9条连续记录，再按字符串正向检查9-mer重叠。该算法适用于与这一遍历方向一致的9-mer。**对于当前提取器输出的RNA方向负链9-mer，直接套用可能没有有效窗口或遗漏负链候选。** 本次人工正负链例子复现：正链生成1个候选；对应负链报 `No valid 9-row windows were generated from the input.`。

因此，不将该命令描述为当前RNA004原始数据的全链条通用转换方案，也不能用只留下正链来冒充完整处理。当前 `FINAL_TOOL/prepare_features.py` 同样调用该核心构造器，并非已修复的负链替代品。

两条正确的衔接原则：

- **复查已有catalog验证结果**：使用前述已验证prepared表。现有 `prepare_model_input.py` 按 `ref_start+4` 连接catalog，核对中央9-mer与正/反向参考，并同步数值方向，用可信的参考上下文构造17 nt；不依赖这个升序字符串窗口假设。
- **处理未知样本**：用同一参考/剪接路径，针对每条Step6的中心 `ref_start+4` 构造RNA方向17 nt，核对其中间9 nt与Step6 `sequence_input`，保留该中央Step6行的36数值；如需反向必须同步反序数值。去掉catalog truth筛选，按自己的候选region定位并记录未通过的行。或适配相邻9行构造器，按目标链方向正确遍历，联合校验参考、坐标、序列和数值。当前ZIP没有附带经过独立验证的任意未知样本全转录组适配器。

上面是处理新样本所需的数据合同和适配步骤，不声称已实现并验证了新的通用原始数据pipeline。本文不更改现有工具与模型代码。

### 5.6 RNA004融合的81 nt

推荐有序转录本FASTA＋转录本坐标，或基因组BAM的read特异成熟RNA参考路径。直接在hg38上连续取81 nt只适用于该窗口不跨剪接边界且方向正确的情况。

现存脚本：

```text
/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/mouse_model_pipeline/add_splice_aware_context.py
```

它沿BAM CIGAR保留参考消耗的 `M/D/=/X`，跳过 `N` 内含子，忽略不消耗参考的 `I/S/H/P`，取中心前后各40个参考碱基；负链反向并互补。输出是参考路径窗口，不是含插入/测序错误的原始read窗口。

但该脚本要求 `reference_sequence_21bp,true_label,site_id` 等catalog字段，并用21 nt核验；其Binary输出还按truth范围筛选。**不应原样用于未知样本，更不能给未知候选伪造truth。** 处理新样本要保留其正确的CIGAR路径逻辑，同时用与候选相关的可信上下文核验，移除truth依赖，保留审计。

统一工具本身的 `--context-fasta` 只做连续FASTA切片，不解析BAM、不自动跨剪接。已有正确 `reference_sequence_81bp` 时直接传该列最清楚。

## 6. 安装、八种主要组合与实用命令

### 6.1 安装与检查

在新的安装目录解压后进入 `NanoRemix/`：

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[models,second]' -c constraints-models.txt
nanoremix --version
nanoremix catalog
nanoremix doctor --model-root models --model-checksums models/SHA256.json --strict
```

源码依赖numpy/pandas；三代需要scipy、scikit-learn1.3.0、XGBoost2.1.4、joblib、NumPy<2；序列二代另需torch。原包建议Python3.10/3.11。当前超算模型环境可参考：

```text
/gpfs/work/bio/zhixingwu16/.conda/envs/nanoremix_env/bin/python
Python3.10.18 / numpy1.24.3 / pandas2.0.3 / scipy1.10.1 /
scikit-learn1.3.0 / xgboost2.1.4 / torch2.2.2+cu118 / joblib1.5.2
```

上面是本次核验使用的现存超算环境，和ZIP内推荐锁定文件不完全相同。ZIP的 `constraints-models.txt` 明确固定：

```text
numpy1.26.4 / pandas2.3.1 / scipy1.11.4 / scikit-learn1.3.0 /
xgboost2.1.4 / joblib1.5.1 / torch2.7.1
```

新安装按该constraints文件执行；本次格式/概率核验不构成对所有依赖组合的逐概率等效性证明。

`doctor`检查文件、metadata、hash和依赖/native runtime，不运行权重。`validate`接受predict同样的输入参数，实际构造特征/加载所需权重/推理，但不保存预测。仅有CLI help或静态检查通过，不等于真实数据准确性验证。

序列推理和原始预处理应放在合适的计算节点执行；不要在登录节点运行全量数据。

### 6.2 八种主要组合的命令选择

统一模板如下。`--mods`按下表选，输出目录必须是新的，模型目录不得用作输出目录。

```bash
# 纯三代模板
nanoremix predict --input /path/to/PREPARED.csv --output /path/to/NEW_RESULTS \
  --chemistry RNA002 --task multi --mods all --mode thirdgen --strategy s2 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files

# 融合模板：表内已有81 nt
nanoremix predict --input /path/to/PREPARED_WITH81.csv --output /path/to/NEW_RESULTS \
  --chemistry RNA002 --task multi --mods all --data-type biological --strategy s2 \
  --sequence-column reference_sequence_81bp \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

| 组合 | `--chemistry` | `--task` | 推荐明确的 `--mods` | 模式与输入 |
|---|---|---|---|---|
| 002 Binary纯三代 | RNA002 | binary | `m6A,m5C,pseudoU`或子集 | `--mode thirdgen`；002 prepared |
| 002 Binary融合 | RNA002 | binary | 同上 | `--data-type biological`＋81 nt/二代概率；002 prepared |
| 002 Multi纯三代 | RNA002 | multi | `all`或支持修饰子集 | `--mode thirdgen`；002 prepared |
| 002 Multi融合 | RNA002 | multi | 同上 | `--data-type biological`＋81 nt/二代概率；002 prepared |
| 004 Binary纯三代 | RNA004 | binary | `m6A,m5C,pseudoU`或子集 | `--mode thirdgen`；004 prepared |
| 004 Binary融合 | RNA004 | binary | 同上 | `--data-type biological`＋81 nt/二代概率；004 prepared |
| 004 Multi纯三代 | RNA004 | multi | `all`或支持修饰子集 | `--mode thirdgen`；004 prepared |
| 004 Multi融合 | RNA004 | multi | 同上 | `--data-type biological`＋81 nt/二代概率；004 prepared |

为了比较同一研究范围的融合与两代分量：

```bash
nanoremix predict --input sample004.with81.csv --output results/004_multi_compare \
  --chemistry RNA004 --task multi --mods m6A,m5C,pseudoU --data-type biological \
  --strategy both --sequence-column sequence81 --also-thirdgen --also-secondgen \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

该命令写5张表：融合s1/s2、纯三代s1/s2、纯二代none，不做S1/S2共识。若希望比较全原生三代与common融合，需注意同一次运行的附加三代仍使用该次统一选定的类别；要全原生类别应另行运行纯三代。

### 6.3 其他二代来源的命令

FASTA已有转录本方向、坐标为0基中心时：

```bash
nanoremix predict --input sample004.prepared.csv --output results/004_fasta_fused \
  --chemistry RNA004 --task multi --mods m6A,m5C --data-type biological \
  --context-fasta transcripts.fa --fasta-id-column reference_name \
  --center-column transcript_center0 --coordinate-base 0 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

如果FASTA方向需要反向互补，加正确的 `--strand-column`。这里 `transcript_center0` 须是真实的转录本中心坐标，不能填基因组坐标。

已算二代完整概率：

```bash
nanoremix predict --input sample002.prepared.csv --output results/002_probability_fused \
  --chemistry RNA002 --task binary --mods m6A,m5C,pseudoU --data-type biological \
  --second-probabilities second.tsv --join-on record_key --strategy s2 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

上式仍加载三代模型，但不加载二代网络。若三代也已算成完整/对应pair概率，改 `--input-kind probabilities` 并使用符合第3.4节的主表，则两代均可不加载权重。

纯二代序列预测（不需要三代信号值）：

```bash
nanoremix predict --input sequence_candidates.csv --output results/sequence_only \
  --chemistry RNA002 --task multi --mods all --mode secondgen \
  --sequence-column sequence81 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

`sequence_candidates.csv` 至少含有效 `center_ref_base,sequence81`，建议含record_key。这里chemistry显式声明用于接口，二代网络本身为共享序列模型；提供短上下文时仍检查其长度与对应chemistry一致。

### 6.4 一条现存真实输入的衔接示例

以下是可供在计算节点使用的模板，本文没有执行全量样本：

```bash
nanoremix predict \
  --input /gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/HeLa/fused_model_results/input/multi_human_fused_input.csv.gz \
  --output /path/to/NEW_HeLa_v030_results \
  --chemistry RNA004 --task multi --mods m6A,m5C,pseudoU \
  --data-type biological --strategy s2 --sequence-column reference_sequence_81bp \
  --also-thirdgen --also-secondgen --on-invalid mark \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

## 7. 条件、筛选与全部predict/validate参数

### 7.1 筛选发生的顺序

```text
格式/chemistry/身份核对
  → 正常参考上下文＋有效数值
  → 确定研究类别、限制参考碱基家族
  → 运行所需三代及二代来源
  → 在研究集合内条件化与归一化，或融合
  → 检查所选质量和二代质量
  → argmax得到候选，检查最高概率confidence
  → status=ok才给predicted_label
  → 按output-filter选择写入预测表的行；完整审计另写
```

`--mods m6A,m5C` 默认加A和C，只允许A家族或C家族候选参与相应竞争；G/U中心记录为 `outside_selected_family`。`--canonical all` 另保留G/U正常对照，这些没有选中修饰的家族只会有正常类，confidence=1不代表评估了其所有修饰。

Multi `--targets A,m6A,m1A` 指精确竞争类别，不自动补正常碱基；不能与mods并用。一个家族只选一个修饰而无其他竞争类会报错；可以选择同一家族两个修饰作明确的条件比较，但这时没有正常对照。`--canonical`只影响mods路线。

**先选类别、再归一化、再argmax**。它不是全类别argmax后再删掉未研究的标签。遗漏真实修饰可能被重新分配给所选类别，故confidence是条件概率，不能直接解释成绝对修饰率。

### 7.2 阈值定义与默认值

| 参数 | 默认 | 实际含义 |
|---|---|---|
| `--min-confidence` | 0 | 获胜类别的条件概率下限，可为正常类；不是专门的P(mod)阈值 |
| `--min-selected-mass` | 0 | 在该参考碱基家族内，所选类别保留了原生分布多少质量 |
| `--min-second-mass` | 0 | 二代相应所选集合占二代该家族的质量下限 |

三个阈值都在[0,1]，只有**严格低于**阈值才拒绝。默认0没有额外置信度筛选，但仍受格式、有效支持、家族和序列一致性检查。

设原生概率p，参考家族b、研究集合T，则：

\[
\mathrm{selected\_family\_mass}=
\frac{\sum_{c\in T\cap F_b}p_c}{\sum_{c\in F_b}p_c}.
\]

Multi纯三代/融合用三代原生分布计算上式；纯二代用二代分布。Binary的三代pair为自身完整比较，所以 `selected_family_mass=1`；其 `second_selected_family_mass` 为二代选中pair占二代整个参考家族的质量。不是pair占全部15类的未条件化质量，也不是样本修饰比例。

Multi融合中二代不能提供有效竞争的家族，输出会注明三代证据；部分未参与共有概率更新的行的二代mass指标为实现默认值1。不能把这解释成“二代模型对全部类别支持充分”。应同时查看 `evidence_source`。

示例：C家族原分布 `C=.2,m5C=.2,hm5C=.6`，只选C/m5C后为 `.5/.5`，但保留质量仅.4；`--min-selected-mass 0.5`拒绝这一记录。阈值应由与最终测试独立的验证数据预先确定，本文不提供未经验证的通用最优值。

### 7.3 全部predict/validate选项

`validate`与`predict`共享以下解析选项，但通常不传output，且不保存预测。完整原始help在 `CLI_help.txt`。

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

### 7.4 错误与输出过滤

常见status：`outside_selected_family`、`invalid_context_or_reference`、`invalid_signal`、`missing_second_match`、`invalid_second_probability`、`insufficient_81nt_flanks`、`81nt_center_mismatch`、`81nt_context_mismatch`、`invalid_or_zero_probability_support`、`below_selected_mass`、`below_second_mass`、`below_confidence`、`unsupported_by_second_generation`。

`--on-invalid mark` 可保留实现支持的行级失败并输出原因；重复record_key、缺必需列、未知类、二代连接键不唯一、模型哈希/兼容性错误等仍会中止。阈值拒绝和家族外记录属于普通审计排除，不等同于运行失败。

阈值拒绝时可能保留 `candidate_label`、confidence及概率，但 `predicted_label` 为空。`--output-filter modified` 仅留下已接受的修饰预测，不改变用于评价的完整运行记录。每个输出目录都有完整 `selection_audit.tsv`，但该审计表不保存全部prob列。

CLI没有 `--motif`、`--min-coverage`、`--mapq`、`--site-threshold`、`--top10` 选项，不自动筛DRACH、不自动做site覆盖度过滤，也不自动把read结果聚合成位点修饰率。论文的S1/S2共识、top10%、confidence/target-mass队列规则不能不加说明套成此CLI默认行为。

## 8. 具体怎么算

### 8.1 三代特征与模型

序列位置以正常ACGU one-hot编码，不把输入中修饰标签当作真值。Binary中心类别已由专用pair固定，主要编码侧翼；Multi额外编码中心参考家族。

| 模型 | S1每组件的输入 | S1组件数 | S2的输入 |
|---|---|---|---|
| RNA002 Binary | 4个侧翼位置one-hot＋对应槽位4个统计，共20维 | 每pair5个 | 8侧翼位置one-hot＋20数值，共52维 |
| RNA002 Multi | 上述位置再加中心参考one-hot，共24维 | 5个 | 加中心参考one-hot，共56维 |
| RNA004 Binary | 8侧翼位置one-hot＋对应槽位4统计，共36维 | 每pair9个 | 16侧翼位置one-hot＋36数值，共100维 |
| RNA004 Multi | 加中心参考one-hot，共40维 | 9个 | 加中心参考one-hot，共104维 |

S1不同组件使用不同的滑移位置子集及相应统计槽位。S2使用全部规定侧翼和数值。实际特征顺序从schema/model feature names核对，不自行排序列后硬传裸数组。

### 8.2 Binary三代

每个校准组件输出P(mod)；从 `classes_` 找阳性类，不猜概率列顺序。S1将组件pair概率做归一化乘积，数值实现为：

\[
p_3=\operatorname{sigmoid}\!\left(\sum_j\operatorname{logit}(p_{3,j})\right).
\]

S2直接用对应 `m_all_calibrated.pkl` 的pair概率。logit前将概率clip到 `[10^{-15},1-10^{-15}]`。本v0.3.0的 `strategy both` 分别计算两套，不把它们再乘一次。

### 8.3 Multi三代

组件校准后的原生类概率，先在参考家族内取值并归一化，记 \(p_{m,c}\)。各家族有冻结的聚合权重 \(w_{b,m}\)、温度 \(T_b\) 和method，来自所选strategy的manifest。

几何聚合：

\[
z_c=\sum_m w_{b,m}\log p_{m,c},\qquad
r_c=\operatorname{softmax}_{c\in F_b}(z_c/T_b).
\]

算术聚合时改为 \(z_c=\log(\sum_mw_{b,m}p_{m,c})\)。本包四份实际manifest均为geometric，gate=true；精确参数见 `frozen_multi_fusion.json`，不是统一写死为等权乘积。

例如A家族的冻结温度为：RNA002 S1≈0.610095、S2≈1.349353；RNA004 S1≈0.946026、S2≈1.515072。S2虽只有一个组件，也仍要执行其冻结temperature。不要以为单组件便可忽略manifest。

这些是**三代模型内部聚合参数**，与CLI二代/三代融合的1/1权重是两层不同的参数。

### 8.4 共享二代网络

81 nt经原网络定义编码，再由5 repeat×5 fold共25个CNN-GRU-attention模型推理。每个checkpoint有自己的正温度 \(T_j\)：

\[
p_{2,c}=\frac1{25}\sum_{j=1}^{25}
\operatorname{softmax}(\mathrm{logits}_j/T_j)_c.
\]

顺序是**各自temperature-softmax→25份概率平均**，不是先平均logits、也不是只取一个fold。重复81 nt会在处理批次内复用推理结果；这不改变输入记录身份。模型结构、类顺序、81 nt长度、checkpoint数量和state_dict均核对。

### 8.5 Binary二代pair与融合

对目标家族的正常类b和修饰m，先从完整二代分布取pair并归一化：

\[
p_2=\frac{P_2(m)}{P_2(b)+P_2(m)}.
\]

其余修饰概率**不并入正常类**。随后：

\[
p_f=\operatorname{sigmoid}\bigl(w_2\operatorname{logit}p_2+
w_3\operatorname{logit}p_3+\beta\bigr).
\]

默认w2=w3=1、bias=0，等价于：

\[
p_f=\frac{p_2p_3}{p_2p_3+(1-p_2)(1-p_3)}.
\]

例：p2=.8、p3=.6，融合为 `.48/(.48+.08)=.857142857`；最终输出正常类.142857143、修饰.857142857。这是归一化概率乘积，不能简写成未归一化的 `.8×.6=.48`。

### 8.6 Multi选择、共有类别与扩展融合

先在选中集合T且参考家族b内条件化三代概率，得到r。设该家族所选共有类别为S，\(M=\sum_{c\in S}r_c\)。将二代在S内归一化为q2，共有部分计算：

\[
u_c=\frac{(r_c/M)^{w_3}q_{2,c}^{w_2}}
{\sum_{d\in S}(r_d/M)^{w_3}q_{2,d}^{w_2}},\quad c\in S.
\]

extended输出：

\[
P_f(c)=M u_c\ (c\in S),\qquad P_f(c)=r_c\ (c\notin S).
\]

因此三代独有类别保留原质量，共有部分只重新分配其自身的M。common模式先限制到共有类别，M=1。如果只剩一个正常共有类而修饰全为三代独有，这个家族没有可由二代提供的有效修饰对照，证据标签会注明三代来源。两个分布完全没有共同正支持时不输出人为均匀分布，而是标记无效。

数值例：选择 `C,m5C,hm5C`，三代在此家族分布为 `.2,.5,.3`；二代在C/m5C内为 `.4,.6`。M=.7，共有部分乘积归一化后为 `4/19,15/19`，extended最终：

```text
C     = .7 × 4/19  = .147368421
m5C   = .7 × 15/19 = .552631579
hm5C  = .3
```

hm5C没有获得二代概率修正。最高值confidence=.552631579；若min-confidence=.6则拒绝该记录。

### 8.7 最终决策、阈值与评价

跨参考家族的概率为0。最终candidate是所选概率的argmax；精确平局按固定类别顺序取首项，普通Binary pair中因此正常类优先。之后仅status=ok给接受标签。

`truth-column`只用于统计，绝不拟合模型、temperature、融合权重或阈值。输出的评价为所选truth范围内已保留记录的混淆矩阵、accuracy、precision/recall/F1；同时记范围外truth和拒绝数，并报告“范围内拒绝计错”的accuracy。它不是无条件全修饰性能；本CLI不自动计算AUROC/AP和论文的benchmark汇总。

## 9. 输出文件、核验与阅读顺序

```text
NEW_RESULTS/
  run_manifest.json
  fused/s1/predictions.tsv               # 按所选模式/策略出现
  fused/s1/selection_audit.tsv
  fused/s2/...
  thirdgen/s1/...                        # 纯三代或附加输出
  thirdgen/s2/...
  secondgen/none/...                     # 纯二代只写一次
```

主要字段：身份/坐标、chemistry、task、mode、strategy、status、candidate_label、predicted_label、confidence、selected_family_mass、second_selected_family_mass、evidence_source和所选prob列。

`run_manifest.json`包含输入hash、实际加载资产hash、CLI参数、软件版本、策略映射、研究范围、各模式类轴、阈值、固定融合系数、计数和可选评价。`model_retrained=false`；推理不会自动重训或“寻找最优参数”。

建议先读第2–3节确定版本和输入，再按第4/5节核对原始特征来源，第6节运行，第7–8节解释条件概率、筛选与输出。

本次附件：

- `docs/OVERVIEW.zh-CN.md`：四套模型、八种主要工作流、输入、准备流程、参数、筛选、计算和发布权重的中文总览。
- `CLI_help.txt`：从此ZIP的解析器生成的原始命令帮助。
- `PIPELINE_命令与输入衔接.md`：原始数据处理的命令模板及RNA002坐标连接示例；明确RNA004通用转换的适用边界。
- `catalog.json`：精确原生类轴、共有/独有范围。
- `frozen_multi_fusion.json`：四组manifest的完整内部聚合权重/温度。
- `model_path_mapping.tsv`：89个权重在ZIP与当前正式目录间的对应关系和hash核验。
- `provenance.json`：被核对ZIP及上游脚本hash、版本记录。
- `examples/`：明确标注为人工构造的CSV/TSV/FASTA，用于输入合同和参数演示。
- `verification/README.md`：精简包中保留的分发完整性检查方法和验证范围说明；自动生成的运行报告不随包分发。

## 10. 核对所依据的代码

ZIP内部：

| 主题 | 代码 |
|---|---|
| CLI参数、prepare | `src/nanoremix/cli.py` |
| 类范围、别名、s2↔s3 | `src/nanoremix/labels.py` |
| 输入、坐标、81 nt核验 | `src/nanoremix/tables.py` |
| 分支模型加载、二代25成员 | `src/nanoremix/engines.py` |
| 概率条件化、pair与Multi融合 | `src/nanoremix/probability.py` |
| 模式路由、筛选、输出、评价 | `src/nanoremix/pipeline.py` |
| 5行/9行候选构造 | `src/nanoremix/_vendor/binary5/features.py`、`binary9/features.py` |
| 原始Step4/Step6格式 | `src/nanoremix/_vendor/binary5/io.py`、`binary9/io.py` |
| Multi特征、聚合与metadata | `src/nanoremix/_vendor/multi5/`、`multi9/`及对应manifest |
| 二代网络定义 | `models/secondgen/scripts/train.py` |

现存超算上游：

- RNA002：`Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/config.sh`、Step1–4 sbatch和 `python/extract_tombo_step4.py`。
- RNA004：`NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/` 的config、01–06脚本、`merge_step6.py`、`prepare_model_input.py`。
- RNA004原始特征：`validation_RNA004/preprocess_pipeline/feature_extraction_rna004.py`。
- 剪接参考窗口：`validation_RNA004/mouse_model_pipeline/add_splice_aware_context.py`。

所有相对超算路径均以 `/gpfs/work/bio/zhixingwu16/` 为根。历史结果目录名、文件名中的“final”不替代模型hash、输入定义和证据范围核对。
