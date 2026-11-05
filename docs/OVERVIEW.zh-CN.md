# NanoRemix v0.3.0 完整功能与使用说明

本文件是源码包内的中文总览，集中说明 NanoRemix 的模型配置、运行路线、输入文件、上游数据准备、命令参数、筛选和计算。逐项 CLI help、复杂坐标规则和完整推导仍以本目录链接的专题文档及实际代码为准。

## 1. 版本数量与名称

NanoRemix 统一推理工具版本为 **v0.3.0**。这里有 **4 套第三代 chemistry/task 模型配置**，并非 4 个不同的软件安装包：

| 编号 | chemistry | task | 三代模型输入 | 三代预训练对象 |
|---|---|---|---|---:|
| 1 | RNA002 | Binary | 9 nt 序列上下文 + 20 个信号数值 | 18 个 `.pkl` |
| 2 | RNA002 | Multi | 9 nt 序列上下文 + 20 个信号数值 | 6 个 `.pkl`，另有 S1/S2 manifests |
| 3 | RNA004 | Binary | 17 nt 序列上下文 + 36 个信号数值 | 30 个 `.pkl` |
| 4 | RNA004 | Multi | 17 nt 序列上下文 + 36 个信号数值 | 10 个 `.pkl`，另有 S1/S2 manifests |

四套三代配置都能运行纯三代或二代+三代融合，共 **8 条主要预测路线**。每条可选 S1 或 S2，所以有 16 种 route/strategy 配置；`--strategy both` 是一次分别跑 S1 和 S2，不计算两者共识。另有共享的二代序列模型，可单独输出二代结果。所有路线由同一个 `nanoremix` CLI 管理。

“二代/三代”是项目内部模型分支名称，不代表测序技术代际。二代使用 81 nt 连续序列；三代使用相应 chemistry 的序列和 Nanopore 电流特征。**Binary/Multi** 是分类任务类型，**RNA002/RNA004** 是 RNA 化学版本，**S1/S2** 是模型策略；这些维度不能互相替代。

旧模型文件名 `strategy1` 对应 CLI `s1`，`strategy3` 对应 CLI `s2`。CLI 接受旧别名 `--strategy s3` 并将其解释为 `s2`。不应把 S1/S2 称为二代/三代融合；两代融合是另一层独立操作。

## 2. 四套配置具体预测什么

### 2.1 RNA002 Binary

将每个支持的修饰作为一个专用二分类：A 对 m6A、C 对 m5C、U 对 pseudoU。`--mods` 可选择这些类别的子集或 `all`。预测竞争发生在对应正常碱基与该专用修饰之间；不是 Multi 的全类别分数。

### 2.2 RNA002 Multi

在参考碱基家族内进行多类竞争。原生三代类别共 18 类（包括 A/C/G/U）：

- A 家族：A、m6A、m1A、Inosine、Am。
- C 家族：C、m5C、hm5C、ac4C、FormylC。
- G 家族：G、m7G、o8G。
- U 家族：U、pseudoU、m1psi、m5U、Um。

### 2.3 RNA004 Binary

与 RNA002 Binary 支持同一组三种专用 pair：A/m6A、C/m5C、U/pseudoU；输入特征长度和对应权重换为 RNA004。

### 2.4 RNA004 Multi

在参考碱基家族内进行多类竞争。原生三代类别共 16 类（包括 A/C/G/U）：

- A 家族：A、m6A、m1A、Am。
- C 家族：C、m5C、hm5C、ac4C。
- G 家族：G、o8G、Gm。
- U 家族：U、pseudoU、m1psi、m5U、Um。

共享二代网络的 15 类固定为 `A,C,G,U,Am,m1A,m6A,Inosine,ac4C,m5C,Gm,m7G,m5U,pseudoU,Um`。默认 `--fusion-scope extended` 会融合共有类别，并保留三代独有类别：RNA002 的 hm5C、FormylC、o8G、m1psi；RNA004 的 hm5C、o8G、m1psi。这些独有类别沿用三代证据，不能说有二代支持。`--fusion-scope common` 只报告共有轴（RNA002 为 14 类，RNA004 为 13 类）。`Cm` 不属于随包正式支持范围。

## 3. 纯三代、二代+三代融合和纯二代

| 路线 | 主输入 | 第二数据来源 | 必需模型资产 |
|---|---|---|---|
| 纯三代 | 对应 RNA002/RNA004 的 prepared 特征表，或三代预计算概率 | 无 | 对应 chemistry/task 的三代包；若主输入是三代概率，则无需加载三代模型 |
| 二代+三代融合（现场跑二代） | 同一条三代 prepared 特征/概率记录 | 同候选、同方向的连续 81 nt 序列，或能按身份和坐标抽取该序列的 FASTA | 对应三代包 + 共享 SecondGen 包（25 个 checkpoint） |
| 二代+三代融合（已有二代概率） | 同一条三代 prepared 特征/概率记录 | 完整 15 类二代概率，内嵌或按唯一身份键连接 | 对应三代包；不用下载 SecondGen checkpoint |
| 纯二代 | 中心参考碱基与连续 81 nt，或完整二代概率 | 无 | 现场推理时用 SecondGen；二代概率输入时不需 checkpoint |

生物样本的默认路线为融合，必须能取得两代证据；工具不会悄悄降级成纯三代。生物样本纯三代请显式使用 `--mode thirdgen`。`--data-type synthetic` 表示合成数据常用的纯三代路由，不要为规避二代输入把真实样本标成 synthetic。

融合时可以加 `--also-thirdgen --also-secondgen` 同时保存两条组成分支。`--mode both` 表示分别跑纯三代和融合（仍须二代证据）；不包含纯二代模式。`--mode secondgen` 单独运行二代。`--data-type` 与 `--mode` 互斥。

纯二代路线不读取三代信号特征，不需要三代权重；输入要有中心碱基与连续 81 nt 序列，或可替代的完整二代概率。其类别仅限共享 15 类二代轴。二代-only 常规调用使用 `--task multi`；需要同时评估三代和二代时，更常见的用法是融合路线。

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

**负链限制必须注意：**现有 v0.3.0 `nanoremix prepare --format step6` 按参考坐标升序检查相邻 9-mer 重叠；当前 Remora 已把负链序列和信号槽位按 RNA 方向保存，直接调用旧转换器可能漏掉负链候选。不能把它称为已验证的正负链通用入口。已有 `prepare_model_input.py` 和 splice-aware 上下文脚本依赖已知 catalog/truth 字段，也不能原样用于未知样本。新样本需按 `ref_start+4` 构造上下文、校验中央序列、同步反转序列与数值槽位（需要反向时），并记录每条拒绝原因；当前包没有声称包含经过验证的任意未知样本全转录组适配器。

RNA004 的参考序列、regions、BAM 和 FASTA 必须对应同一参考版本。9-mer region 文件格式和区间是候选 9-mer 左端起点，不等同于任意 site BED；细节见 [RNA004 pipeline 模板](guide_20261104/PIPELINE_命令与输入衔接.md)。

### 5.3 推理环境

原始权重兼容环境建议 Python 3.10/3.11，scikit-learn 1.3.0、XGBoost 2.1.4。项目的 `constraints-models.txt` 给出安装约束。basecalling/信号处理与推理使用独立环境；包内没有包含 Guppy、Dorado、Tombo、Remora、参考基因组或数据。

## 6. CLI 常用选项与条件

调用形式：`nanoremix predict --input TABLE --output NEW_DIR ...`。每次使用新的输出目录。完整参数以 `CLI_help.txt` 和 `nanoremix predict --help` 为准。

| 选项 | 作用和限制 |
|---|---|
| `--chemistry RNA002|RNA004` | 选择 RNA 化学版本，应与真实实验和模型一致 |
| `--task binary|multi` | 选择四套配置之一 |
| `--mods m6A,m5C` / `--mods all` | 选择支持修饰集合或子集；家族过滤依据类别发生 |
| `--targets A,m6A,m1A` | Multi 精确选择竞争类；和 `--mods` 互斥，且不会自动补正常碱基 |
| `--canonical matching|all` | `mods` 路线的正常碱基对照，默认 matching |
| `--mode thirdgen|fused|both|secondgen` | 显式选择模型分支；both=分别三代与融合，不含纯二代 |
| `--data-type biological|synthetic` | 未显式 mode 时 biological 默认融合、synthetic 默认三代；与 mode 互斥 |
| `--strategy s1|s2|both` | 默认 s2；both 单独输出两种策略，不做共识 |
| `--input-kind features|probabilities` | 主表是 prepared 特征还是三代概率 |
| `--sequence-column COL` | 指定 81 nt 序列列 |
| `--context-fasta FILE` | 使用参考/事件 FASTA；必须同时给 ID、中心列和坐标基准，必要时给链列 |
| `--second-probabilities FILE --join-on KEY` | 读取独立完整二代概率表并按身份键连接 |
| `--fusion-scope extended|common` | Multi 是否保留三代独有类，默认 extended |
| `--second-weight W --third-weight W` | 融合指数权重，必须有限且大于 0，默认均 1 |
| `--binary-bias B` | Binary 融合 log-odds 偏置，默认 0；不用于 Multi |
| `--min-confidence X` | 获胜类别条件概率下限，默认 0 |
| `--min-selected-mass X` | 选中类在对应参考家族原生分布中保留质量的下限，默认 0 |
| `--min-second-mass X` | 二代选中集合在对应参考家族中保留质量的下限，默认 0 |
| `--output-filter all|retained|modified` | 写出全部记录、通过筛选记录或已接受修饰记录，默认 all |
| `--on-invalid error|mark` | 行级无效记录报错或标状态，默认 error；全局 schema/key 错误仍会中止 |
| `--also-thirdgen --also-secondgen` | 融合运行时可额外保存相应分支结果 |
| `--model-root models` | 模型根目录；解压 Release 权重后使用 |
| `--model-checksums models/SHA256.json` | 校验模型文件索引 |
| `--trust-model-files` | 明确允许加载可信的 pickle/checkpoint 和二代模型代码 |
| `--chunk-size N` / `--batch-size N` / `--device cpu|cuda|mps` | 控制主表分块、二代 batch 与运行设备；不是整个工作流内存上限 |
| `--truth-column COL` | 只做评价统计，不用于拟合或预测 |

`mods` 先限定标签及参考碱基家族，再条件化、融合和 argmax。比如只选 m6A/m5C 时，U/G 候选通常标为 `outside_selected_family`；`--canonical all` 会额外保留未选择修饰的正常碱基记录。Multi 的 `--targets` 是精确竞争集合，不会暗中添加 canonical；只有一个修饰且无同家族对照时会被拒绝为无效设置。

三个质量阈值均须在 0–1 范围内，严格低于阈值才拒绝。默认 0 不附加概率阈值，但输入格式、中心/方向一致性、类别支持及模型完整性检查始终执行。`confidence` 是选定类别集合中获胜类别的条件概率，不是绝对 P(mod) 或位点修饰率；遗漏真实修饰可能令概率质量重分配给已选标签。selected mass 也不是 read coverage。

输出包含 `predictions.tsv`、对应的 `selection_audit.tsv` 和 `run_manifest.json`。audit 记录每行 status/拒绝原因，run manifest 记录选项、输入/模型 hash 和数量。CLI 不自动筛 DRACH/motif、MAPQ、coverage、site threshold，不聚合 read-level 预测成 site-level 修饰率。

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

先检查 schema、chemistry、中心/方向、数值和身份，再确定类别/家族并计算概率；接着计算 selected mass/second mass、argmax 与 confidence；只有 status 为 `ok` 才有接受的 `predicted_label`。阈值拒绝通常保留 `candidate_label` 和审计原因，但 accepted label 为空。精确公式和数值案例见 [详细计算章节](guide_20261104/README_版本输入预处理参数与计算.md)。

## 8. 五个权重 Release 附件

纯三代只取对应的一包；现场执行序列型二代融合时再取共享 SecondGen 包；已有二代概率表不需要 SecondGen。

| Release asset | 包含内容 | 用途 |
|---|---|---|
| `NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip` | 18 个 `.pkl` | RNA002 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip` | 6 个 `.pkl`、S1/S2 manifests | RNA002 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip` | 30 个 `.pkl` | RNA004 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip` | 10 个 `.pkl`、S1/S2 manifests | RNA004 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_SecondGen_20261104.zip` | 网络定义及共享 25 个 `.pt` checkpoint | 二代推理；RNA002/RNA004 共用 |

四组三代共 64 个 `.pkl`，共享二代共 25 个 `.pt`，完整集合为 89 个权重。Binary 消融模型未包含。模型包解压到项目根目录后组成 `models/`。包内源码没有模型权重，不需要重新训练。

## 9. 安装和命令示例

在解压得到的 `NanoRemix/` 项目根目录使用 Python 3.10/3.11：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[models,second]' -c constraints-models.txt
python verification/verify_distribution.py
```

RNA002 Binary 纯三代特征输入：

```bash
nanoremix predict --input prepared002.csv --output results/run_002_binary \
  --chemistry RNA002 --task binary --mods m6A --mode thirdgen --strategy s2 \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

RNA004 Multi 融合，主表有连续 81 nt 的 `sequence81`：

```bash
nanoremix predict --input prepared004.csv --output results/run_004_multi_fused \
  --chemistry RNA004 --task multi --mods m6A,m5C,pseudoU \
  --data-type biological --strategy both --sequence-column sequence81 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

换用 RNA002 Multi 或 RNA004 Binary 时，改变 `--chemistry`、`--task`、支持的标签和对应三代权重。概率输入可使用 `--input-kind probabilities`，不加载被预计算概率替代的那一代模型。每次预测指定新的 `--output` 目录。

## 10. 资料导航和适用范围

- [详细中文方法手册](guide_20261104/README_版本输入预处理参数与计算.md)：全部 41 个 predict/validate 选项、字段、版本证据、完整筛选顺序和公式例子。
- [RNA002/RNA004 上游命令模板](guide_20261104/PIPELINE_命令与输入衔接.md)：命令、Step4/Step6 输入衔接和坐标示例。
- [输入列定义](INPUTS.md)、[模型布局与依赖](MODELS.md)、[CLI 选项说明](OPTIONS.md)、[限制说明](KNOWN_LIMITATIONS.zh-CN.md)。
- [人工输入示例目录](guide_20261104/examples/README.md)：所有示例均为人工数据，不是准确率 benchmark。
- [原始上游副本索引](../upstream_reference/README.md) 与 [`model_path_mapping.tsv`](guide_20261104/model_path_mapping.tsv)：上游脚本/权重路径对应。

包完整性通过 `python verification/verify_distribution.py` 检查只代表文件哈希、清单和缺失文件通过；它不运行模型、不重做原始数据处理，也不证明新样本准确率。原始 [LICENSE](../LICENSE) 标为待权利人确认，公开源码、模型和参考材料前需要明确分发许可。
