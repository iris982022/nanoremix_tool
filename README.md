# NanoRemix

**RNA002 / RNA004 RNA修饰预测：Binary、Multi、纯三代及二代＋三代融合。**

工具版本 **v0.3.1** · Python **3.10/3.11** · 主页说明更新 **2026-11-05**

[快速运行](#quickstart) · [四套配置](#configurations) · [输入](#inputs) · [上游流程](#preprocessing) · [全部参数](#options) · [计算](#calculation) · [权重](#models) · [运行与输出](#examples)

源码与权重分别分发。源码不含模型，运行特征推理前需从 [Releases](https://github.com/iris982022/nanoremix_tool/releases) 下载对应附件。首次可先跑无需权重的人工概率示例。本次维护版恢复缺失文件，修复源码打包和参数说明，并为旧RNA004 Step6转换加输入保护。预测算法与权重不变；模型继续下载v0.3.0 Release的五个ZIP。



### 本次维护版修复

- 恢复 `.gitignore` 与自动测试workflow，并在CI检查分发完整性。
- 统一必填参数、融合权重必须严格大于0及实际参数名 `--condition`。
- 源码构建保留README引用的详细手册、上游脚本和示例；生成当前文件的完整性清单，精简源码仍为99个文件。
- RNA004旧Step6转换遇到负链、缺失/无效方向或同read多个参考时明确中止，避免静默输出不完整结果。已有经核实的负链prepared表仍可用于predict。

版本关系：GitHub当前发布的模型为v0.3.0；本维护源码包为v0.3.1，权重和预测数学算法沿用原版。更新源码不需要重新训练或重传五个模型包。

<a id="quickstart"></a>

## 快速运行：无需权重的示例

获取保留目录结构的源码，进入含 `pyproject.toml` 和 `src/` 的根目录：

```bash
git clone https://github.com/iris982022/nanoremix_tool.git
cd nanoremix_tool

python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
nanoremix --version

nanoremix predict \
  --input examples/minimal_binary_probabilities.csv \
  --output results/quickstart \
  --chemistry RNA004 --task binary --mods m6A \
  --input-kind probabilities --mode thirdgen --strategy s2
```

也可 **Code → Download ZIP**，解压后进入项目根目录，从创建环境开始。Windows PowerShell激活命令为 `.venv\Scripts\Activate.ps1`。每次使用新的输出目录。

在 `results/quickstart/thirdgen/s2/predictions.tsv` 中预期得到：

| record_key | predicted_label | prob_A | prob_m6A | confidence |
|---|---|---:|---:|---:|
| demo-001 | m6A | 0.10 | 0.90 | 0.90 |
| demo-002 | A | 0.80 | 0.20 | 0.80 |

这是人工概率示例，验证软件输入输出路径，不加载真实模型，也不代表生物学准确率。

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

## 9. 安装和四套配置的运行示例

在含 `pyproject.toml` 和 `src/` 的项目根目录，安装真实模型的兼容依赖：

```bash
python -m pip install -e '.[models,second]' -c constraints-models.txt
```

模型依赖约束为numpy 1.26.4、pandas 2.3.1、scipy 1.11.4、scikit-learn 1.3.0、XGBoost 2.1.4、joblib 1.5.1、torch 2.7.1。上游预处理用其各自环境，不能与推理环境混装。

四套纯三代示例均需要自己的prepared表和匹配权重。生物样本也可明确用 `--mode thirdgen` 做纯三代，无需改称synthetic。每次输出到新目录。

### RNA002 Binary：纯三代

```bash
nanoremix predict --input prepared002.csv --output results/002_binary_third \
  --chemistry RNA002 --task binary --mods m6A,m5C,pseudoU \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

### RNA002 Multi：纯三代

```bash
nanoremix predict --input prepared002.csv --output results/002_multi_third \
  --chemistry RNA002 --task multi --mods all \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

### RNA004 Binary：纯三代

```bash
nanoremix predict --input prepared004.csv --output results/004_binary_third \
  --chemistry RNA004 --task binary --mods m6A,m5C,pseudoU \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

### RNA004 Multi：纯三代

```bash
nanoremix predict --input prepared004.csv --output results/004_multi_third \
  --chemistry RNA004 --task multi --mods all \
  --mode thirdgen --strategy s2 --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

### 上述任一种改为二代＋三代融合

将 `--mode thirdgen` 替换为 `--data-type biological --sequence-column sequence81`，使用含正确81 nt的对应prepared表，另准备SecondGen模型。下面以RNA004 Multi为例，并额外保存两代分支结果：

```bash
nanoremix predict --input prepared004_with81.csv --output results/004_multi_fused \
  --chemistry RNA004 --task multi --mods m6A,m5C,pseudoU \
  --data-type biological --sequence-column sequence81 --strategy s2 \
  --also-thirdgen --also-secondgen --model-root models \
  --model-checksums models/SHA256.json --trust-model-files
```

若已有独立二代概率表，用 `--second-probabilities second.tsv --join-on record_key` 替换sequence-column；不需要加载二代checkpoint。

### 解压模型与预检

在项目根目录解压选中的模型ZIP。例如RNA004 Multi融合：

```bash
unzip -o NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip -d .
unzip -o NanoRemix_v0.3.0_models_SecondGen_20261104.zip -d .

nanoremix doctor --model-root models --model-checksums models/SHA256.json \
  --chemistry RNA004 --task multi --mods all --strategy s2 \
  --data-type biological --strict
```

各包合并成 `models/`，其中SHA256索引相同。doctor检查文件、manifest、hash与依赖，不执行真实模型。只预检纯三代资产时，doctor用 `--data-type synthetic` 选择范围；doctor没有 `--mode`。若融合将使用预计算二代概率，给doctor加 `--second-source probabilities`。

### 其他命令

```bash
nanoremix --version
nanoremix catalog --chemistry RNA004
nanoremix predict --help
nanoremix validate --help
nanoremix doctor --help
nanoremix prepare --help
```

`validate`运行所给输入/模型推理，但不保存预测表；使用与predict相同的输入、类别、路由和模型参数，通常省略output。

`prepare`必需 `--input`、`--output`、`--chemistry RNA002|RNA004`、`--format step4|extracted|step6`。RNA002支持step4/extracted，RNA004支持显式声明为正链的step6，负链或方向缺失会报错；输出须为新文件。它只转换已有特征，不执行basecall、比对或信号提取。

### 结果文件

| 文件 | 内容 |
|---|---|
| `thirdgen/s2/predictions.tsv` | 纯三代S2预测 |
| `fused/s2/predictions.tsv` | 融合S2预测 |
| `secondgen/none/predictions.tsv` | 纯二代预测，只保存一份 |
| 各路线旁 `selection_audit.tsv` | 全部输入记录和拒绝原因 |
| `run_manifest.json` | 选项、标签、输入/模型hash和计数 |

S1或both写对应s1目录。请同时查看 `status`、`predicted_label`、`candidate_label`、`confidence`、概率质量和 `evidence_source`。

<a id="references"></a>

## 10. 资料导航和适用范围

- [详细中文方法手册](docs/guide_20261104/README_版本输入预处理参数与计算.md)：全部 41 个 predict/validate 选项、字段、版本证据、完整筛选顺序和公式例子。
- [RNA002/RNA004 上游命令模板](docs/guide_20261104/PIPELINE_命令与输入衔接.md)：命令、Step4/Step6 输入衔接和坐标示例。
- [输入列定义](docs/INPUTS.md)、[模型布局与依赖](docs/MODELS.md)、[CLI 选项说明](docs/OPTIONS.md)、[限制说明](docs/KNOWN_LIMITATIONS.zh-CN.md)。
- [人工输入示例目录](docs/guide_20261104/examples/README.md)：所有示例均为人工数据，不是准确率 benchmark。
- [原始上游副本索引](upstream_reference/README.md) 与 [`model_path_mapping.tsv`](docs/guide_20261104/model_path_mapping.tsv)：上游脚本/权重路径对应。

包完整性通过 `python verification/verify_distribution.py` 检查只代表文件哈希、清单和缺失文件通过；它不运行模型、不重做原始数据处理，也不证明新样本准确率。源码许可状态见 [LICENSE](LICENSE)，模型来源和分发元数据见 [NOTICE.md](NOTICE.md)。
