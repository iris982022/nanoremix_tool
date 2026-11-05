# 上游预处理参考副本

这些是超算现有处理链的原样文件副本，仅补充本地可读代码，不替换统一 `nanoremix` CLI。`SOURCE_MAPPING.json` 记录原路径、包内路径、大小和 SHA256；副本均未修改。命令模板已改为引用这些包内脚本，外部程序和自己的数据路径仍需配置。

RNA002：步骤 1–4 和对应 extractor/validator。原 `config.sh`、Slurm 脚本保留了原集群账户、绝对路径和样本配置，用于追溯，不可原样提交到其他集群。

RNA004：basecall/align/features/merge、catalog prepared 适配及剪接上下文脚本；额外附原使用的 9-mer level 表。Dorado/Guppy/Tombo/Remora 等程序未包含，需独立安装。

**未知样本限制**：`prepare_model_input.py` 与 `add_splice_aware_context.py` 是既有 catalog 验证链，仍依赖 catalog、site_id/类别等信息；后者的 Binary 分支按 truth 过滤。不能把这些条件搬到未知样本预测。负链 Step6 转换及正确 81 nt 取法见[使用边界](../docs/KNOWN_LIMITATIONS.zh-CN.md)。

| 包内文件 | 原路径 |
|---|---|
| [RNA002/config.sh](RNA002/config.sh) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/config.sh` |
| [RNA002/python/extract_tombo_step4.py](RNA002/python/extract_tombo_step4.py) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/python/extract_tombo_step4.py` |
| [RNA002/slurm/01_basecall_rna002.sbatch](RNA002/slurm/01_basecall_rna002.sbatch) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/slurm/01_basecall_rna002.sbatch` |
| [RNA002/slurm/02_make_single_fast5.sbatch](RNA002/slurm/02_make_single_fast5.sbatch) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/slurm/02_make_single_fast5.sbatch` |
| [RNA002/slurm/03_tombo_resquiggle.sbatch](RNA002/slurm/03_tombo_resquiggle.sbatch) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/slurm/03_tombo_resquiggle.sbatch` |
| [RNA002/slurm/04_extract_step4.sbatch](RNA002/slurm/04_extract_step4.sbatch) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/slurm/04_extract_step4.sbatch` |
| [RNA004/preprocess_human_pipeline/config.sh](RNA004/preprocess_human_pipeline/config.sh) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/config.sh` |
| [RNA004/preprocess_human_pipeline/01_basecall.sbatch](RNA004/preprocess_human_pipeline/01_basecall.sbatch) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/01_basecall.sbatch` |
| [RNA004/preprocess_human_pipeline/02_align.sbatch](RNA004/preprocess_human_pipeline/02_align.sbatch) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/02_align.sbatch` |
| [RNA004/preprocess_human_pipeline/03_features.sbatch](RNA004/preprocess_human_pipeline/03_features.sbatch) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/03_features.sbatch` |
| [RNA004/preprocess_human_pipeline/merge_step6.py](RNA004/preprocess_human_pipeline/merge_step6.py) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/merge_step6.py` |
| [RNA004/preprocess_human_pipeline/prepare_model_input.py](RNA004/preprocess_human_pipeline/prepare_model_input.py) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/prepare_model_input.py` |
| [RNA004/preprocess_pipeline/feature_extraction_rna004.py](RNA004/preprocess_pipeline/feature_extraction_rna004.py) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_pipeline/feature_extraction_rna004.py` |
| [RNA004/mouse_model_pipeline/add_splice_aware_context.py](RNA004/mouse_model_pipeline/add_splice_aware_context.py) | `/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/mouse_model_pipeline/add_splice_aware_context.py` |
| [RNA002/python/validate_tombo_fast5.py](RNA002/python/validate_tombo_fast5.py) | `/gpfs/work/bio/zhixingwu16/Nanopore/NanoRemix/12_HEK_HMEC_m6A_pipeline/python/validate_tombo_fast5.py` |
| [RNA004/assets/9mer_levels_v1.txt](RNA004/assets/9mer_levels_v1.txt) | `/gpfs/work/bio/yuxinzhang17/DirectMultiRM/ONT/9mer_levels_v1.txt` |

[完整版本、算法与上游说明](../docs/guide_20261104/README_版本输入预处理参数与计算.md)；[命令模板](../docs/guide_20261104/PIPELINE_命令与输入衔接.md)。

原文件中保留的绝对路径为来源证据，并非打包后可自动访问的远程依赖。level 表来源见映射；原版权/许可未因复制改变，公开分发许可状态见顶层 LICENSE。
