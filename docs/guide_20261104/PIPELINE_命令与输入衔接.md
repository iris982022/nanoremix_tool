# NanoRemix v0.3.0 上游处理命令与输入衔接

本文件补充主说明的第4–6节。命令是供核对的模板，没有在本次工作中执行原始数据计算或提交作业。先替换 `/path/to/`、物种参考和样本配置，给每个样本使用全新的输出目录，在具有所需CPU/GPU/内存的计算节点分阶段运行。各阶段使用独立环境；并非将所有工具装进NanoRemix环境。

## A. RNA002：原始FAST5→Step4

以下保留现有HEK/HMEC链的核心参数；脚本使用包内原样副本，外部可执行程序改为待替换路径。原超算路径见主说明与来源记录。新样本使用自己的RNA002数据和匹配参考。

```bash
set -euo pipefail

TASK_PACKAGE_ROOT=/path/to/NanoRemix
TASK002_RAW=/path/to/RNA002_multi_fast5
TASK002_RUN=/path/to/NEW_RNA002_run
TASK002_REF=/path/to/matching_transcriptome.fa
TASK002_PIPE="$TASK_PACKAGE_ROOT/upstream_reference/RNA002"
TASK002_GUPPY=/path/to/guppy-6.3.8/bin/guppy_basecaller
TASK002_FAST5=/path/to/ont-fast5-api-4.1.1-env/bin/multi_to_single_fast5
TASK002_TOMBO=/path/to/tombo1.0-env/bin/tombo
TASK002_GRAPHMAP=/path/to/graphmap-0.5.2/bin/Linux-x64/graphmap
TASK002_PY=/path/to/RNA002-feature-env/bin/python

test -d "$TASK002_RAW"
test -s "$TASK002_REF"
test ! -e "$TASK002_RUN"
mkdir -p "$TASK002_RUN/01_guppy_rna002"

# 1. Guppy6.3.8，GPU；RNA002参数，保留basecalled FAST5。
"$TASK002_GUPPY" -i "$TASK002_RAW" -s "$TASK002_RUN/01_guppy_rna002" \
  --flowcell FLO-MIN106 --kit SQK-RNA002 --recursive --fast5_out \
  --device cuda:0 --num_callers 4 --cpu_threads_per_caller 1

# 2. ont-fast5-api4.1.1，拆分新的single-read FAST5。
mkdir -p "$TASK002_RUN/02_single_fast5"
"$TASK002_FAST5" -i "$TASK002_RUN/01_guppy_rna002/workspace" \
  -s "$TASK002_RUN/02_single_fast5" -t 8 --recursive

# 3. Tombo1.0＋GraphMap0.5.2。
# 私有FASTA副本供GraphMap建索引；Tombo只修改本run内的single-FAST5副本。
mkdir -p "$TASK002_RUN/03_reference"
cp "$TASK002_REF" "$TASK002_RUN/03_reference/transcriptome.fa"
"$TASK002_TOMBO" resquiggle \
  --basecall-group Basecall_1D_000 --basecall-subgroups BaseCalled_template \
  --failed-reads-filename "$TASK002_RUN/03_failed_reads.tsv" \
  "$TASK002_RUN/02_single_fast5" "$TASK002_RUN/03_reference/transcriptome.fa" \
  --graphmap-executable "$TASK002_GRAPHMAP" --processes 8 --align-processes 8 \
  --overwrite --include-event-stdev

# 不单凭退出码0判断成功，核对corrected Events/Alignment。
"$TASK002_PY" "$TASK002_PIPE/python/validate_tombo_fast5.py" \
  --fast5-dir "$TASK002_RUN/02_single_fast5" --sample-size 2000 \
  --min-ready 10 --min-ready-fraction 0.01

# 4. 原定义5-mer统计＋坐标。--rna在此提取器上，非Tombo参数。
mkdir -p "$TASK002_RUN/04_step4"
rg --files "$TASK002_RUN/02_single_fast5" -g '*.fast5' \
  | LC_ALL=C sort > "$TASK002_RUN/04_step4/fast5.list"
"$TASK002_PY" "$TASK002_PIPE/python/extract_tombo_step4.py" \
  --fast5-list "$TASK002_RUN/04_step4/fast5.list" \
  --output-prefix "$TASK002_RUN/04_step4/sample.step4" \
  --workers 8 --clip 10 --progress-every 10000 --rna
```

资源数字仅示范命令接口，需与自己分配的计算资源一致。上面的Tombo抽检下限是现存脚本的最低就绪检查，不代表高resquiggle成功率要求已满足；还要检查失败记录、提取summary及最终候选覆盖。

产物包括：

```text
04_step4/sample.step4.feature.tsv
04_step4/sample.step4.coordinates.tsv
04_step4/sample.step4.errors.tsv
04_step4/sample.step4.summary.json
```

确切error/summary文件名以提取器stdout为准。feature中每read中心事件索引必须保留。进入已安装NanoRemix v0.3.0的环境后：

```bash
nanoremix prepare --input /path/to/NEW_RNA002_run/04_step4/sample.step4.feature.tsv \
  --format step4 --chemistry RNA002 --output /path/to/sample002.prepared.csv
```

该转换全表读入。大数据按完整read分片转换；不要不留边界地任意按行切片，不要先按中心碱基/真值删除邻居。

## B. RNA002：连接同步坐标后取参考81 nt

对于一个内存可容纳的分片，可用如下显式连接逻辑。大表应改用同样身份规则的分片/数据库连接。该示例不创造truth，也不把event index当转录本坐标。

```python
from pathlib import Path
import pandas as pd

prepared = pd.read_csv('/path/to/sample002.prepared.csv')
coords = pd.read_csv('/path/to/sample.step4.coordinates.tsv', sep='\t')
columns = ['id', 'read_id', 'mapped_strand', 'transcript_id', 'transcript_position0']
if coords['id'].isna().any() or coords['id'].duplicated().any():
    raise ValueError('coordinate id must be nonempty and unique')

# prepare默认保留的中间行id是连接键。侧表字段使用新列名，避免覆盖原有元数据。
payload = coords[columns].rename(columns={
    'read_id': 'context_read_id',
    'mapped_strand': 'context_strand',
    'transcript_id': 'context_reference_name',
    'transcript_position0': 'context_center0',
})
joined = prepared.merge(payload, on='id', how='left', validate='many_to_one', indicator=True)
if not joined['_merge'].eq('both').all():
    raise ValueError('unmatched prepared rows')
joined = joined.drop(columns='_merge')
joined['record_key'] = 'sample002:' + joined['id'].astype(str)
if joined['record_key'].duplicated().any():
    raise ValueError('duplicate record_key')
output = Path('/path/to/sample002.prepared_with_coordinates.csv')
if output.exists():
    raise FileExistsError(output)
joined.to_csv(output, index=False)
```

参考FASTA必须与Tombo使用的是同一份序列/ID，坐标已由sidecar从事件位置转换为参考位置：

```bash
nanoremix predict --input /path/to/sample002.prepared_with_coordinates.csv \
  --output /path/to/NEW_002_fused \
  --task multi --chemistry RNA002 --mods all --data-type biological --strategy s2 \
  --context-fasta /path/to/matching_transcriptome.fa \
  --fasta-id-column context_reference_name --center-column context_center0 \
  --coordinate-base 0 --strand-column context_strand \
  --model-root models --model-checksums models/SHA256.json --trust-model-files \
  --on-invalid mark
```

0基中心前后40 nt侧翼不足、N、方向或中央9 nt不符会被标记。人工核对拒绝原因后再解释结果。对已按转录方向的事件FASTA，应使用事件ID和event index，并且不要再套用基因组负链方向。

## C. RNA004：POD5→Step6的核心命令

下面给出单个 `split0` 的处理模板。多分片对每个split运行basecall/align/features，再统一merge。脚本约定POD5在 `pod5/splitN/`，BAM为 `aligned/splitN_sorted.bam`；将自己已有POD5按此约定组织。

```bash
set -euo pipefail

TASK_PACKAGE_ROOT=/path/to/NanoRemix
TASK004_RUN=/path/to/NEW_RNA004_run
TASK004_REF=/path/to/matching_genome.fa
TASK004_REGIONS=/path/to/candidate_kmer_start_regions.csv
TASK004_DORADO=/path/to/dorado-1.3.0-linux-x64/bin/dorado
TASK004_MODEL=/path/to/dorado-models/rna004_130bps_hac@v5.2.0
TASK004_LEVELS="$TASK_PACKAGE_ROOT/upstream_reference/RNA004/assets/9mer_levels_v1.txt"
TASK004_PY=/path/to/remora-3.3.0-env/bin/python
TASK004_FEATURE="$TASK_PACKAGE_ROOT/upstream_reference/RNA004/preprocess_pipeline/feature_extraction_rna004.py"
TASK004_MERGE="$TASK_PACKAGE_ROOT/upstream_reference/RNA004/preprocess_human_pipeline/merge_step6.py"

# 假定原始POD5已放入本次新run的pod5/split0；此处不复用旧计算输出。
test -d "$TASK004_RUN/pod5/split0"
test -s "$TASK004_REF"
test -s "$TASK004_REGIONS"
test ! -e "$TASK004_RUN/basecall"
mkdir -p "$TASK004_RUN/basecall" "$TASK004_RUN/aligned" "$TASK004_RUN/features" "$TASK004_RUN/step6"

# samtools需由当前作业环境/module提供，实际保存版本号。
samtools --version > "$TASK004_RUN/samtools.version.txt"

# 1. Dorado1.3.0 / RNA004 v5.2.0，必须保留moves。
"$TASK004_DORADO" basecaller --device cuda:0 --batchsize 512 --chunksize 2000 \
  --emit-moves --recursive "$TASK004_MODEL" "$TASK004_RUN/pod5/split0" \
  > "$TASK004_RUN/basecall/split0.bam"
samtools quickcheck -v -u "$TASK004_RUN/basecall/split0.bam"

# 2. 匹配参考的主比对，MAPQ≥10。
"$TASK004_DORADO" aligner "$TASK004_REF" "$TASK004_RUN/basecall/split0.bam" \
  | samtools view -@ 8 -bh -F 2308 -q 10 - \
  | samtools sort -@ 8 -o "$TASK004_RUN/aligned/split0_sorted.bam" -
samtools quickcheck -v "$TASK004_RUN/aligned/split0_sorted.bam"
samtools index "$TASK004_RUN/aligned/split0_sorted.bam"

# 3. Remora3.3.0，9-mer，候选起点step=1。
"$TASK004_PY" "$TASK004_FEATURE" \
  --pod5_dir "$TASK004_RUN/pod5" --bam "$TASK004_RUN/aligned" \
  --reg "$TASK004_REGIONS" --level "$TASK004_LEVELS" \
  --output "$TASK004_RUN/features" --threads 8 --run-id sample004 \
  --splits 0 --kmer 9 --step 1 --resquiggle True

# 4. NPZ与sidecar严格对应，stat字段按原定义映射。
"$TASK004_PY" "$TASK004_MERGE" --input-dir "$TASK004_RUN/features" \
  --output "$TASK004_RUN/step6/data_df.csv"
```

reg是CSV，列名必须有 `seqnames,start,end,strand`。代码将 `[start,end)` 当作9-mer**左端起点集合**。0基目标p若需9行重叠重建17 nt，提取起点范围为 `[p-8,p+1)`，而不是简单把单碱基位点p照搬成region；详见主说明。

不要把Step6叫作prepared。其序列为9 nt，下一步还需17 nt参考上下文和正确的方向/中心，融合另需81 nt。

## D. RNA004：Step6→prepared→融合的正确边界

**现有catalog验证链**已经执行了以下步骤：

1. `prepare_model_input.py` 用 `center_pos0=ref_start+4` 连接同参考、同链catalog；核对9 nt并同步调整数值槽位；从可信21 nt参考上下文构造17 nt。
2. `add_splice_aware_context.py` 用read的BAM CIGAR路径得到RNA方向81 nt，跳过N内含子，核对同中心21 nt。
3. 保存带 `9mer_final`、36个数值和 `reference_sequence_81bp` 的prepared表。
4. 统一v0.3.0使用 `--sequence-column reference_sequence_81bp`；不要套用旧Binary/Multi脚本参数。

可直接衔接已存在的HeLa验证输入（在解压后的模型根运行）：

```bash
nanoremix predict \
  --input /gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/HeLa/fused_model_results/input/multi_human_fused_input.csv.gz \
  --output /path/to/NEW_HeLa_multi_fused \
  --task multi --chemistry RNA004 --mods all --data-type biological --strategy s2 \
  --sequence-column reference_sequence_81bp --also-thirdgen --also-secondgen \
  --model-root models --model-checksums models/SHA256.json --trust-model-files
```

这些现有脚本包含catalog/truth筛选。未知样本需要用自己的候选region/参考路径构造上下文、去掉truth依赖，不应直接执行原catalog `submit.sh` 就声称完成通用预测。

**旧 `nanoremix prepare --format step6` 不能作为负链通用保证**：它升序遍历ref_start并检查正向9-mer重叠，对当前RNA方向负链数据可能丢候选。本说明附带了可复现例子。正确的新样本适配必须按RNA方向联合处理坐标、上下文和数值槽位，或用可信参考/剪接路径直接构造17 nt；当前ZIP没有已经验证的未知样本一键适配器。

这一区分也意味着：本文提供了已执行链的正确介绍、可衔接输入及新样本的适配合同，没有把尚未完成验证的原始数据链包装为已支持的功能。原始模型和代码包均保持原样。
