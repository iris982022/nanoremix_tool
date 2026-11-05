# NanoRemix v0.3.0 补充说明

本目录围绕原始 `NanoRemix_complete_v0.3.0.zip` 中的 v0.3.0 推理代码整理，日期为2026-11-04。它区分原始数据处理、候选构造、模型推理和结果筛选；没有重新训练或更换模型。

- [完整中文手册：版本、输入、预处理、参数、筛选和公式](README_版本输入预处理参数与计算.md)
- [全工具功能与四套模型配置中文总览](../OVERVIEW.zh-CN.md)
- [原始处理命令模板与输入衔接](PIPELINE_命令与输入衔接.md)
- [精简包内人工输入示例](examples/README.md)
- [CLI完整命令帮助](CLI_help.txt)
- [精确类别清单](catalog.json)
- [Multi内部冻结的聚合参数](frozen_multi_fusion.json)
- [89个权重的路径与哈希对应表](model_path_mapping.tsv)
- [原始包和上游脚本来源记录](provenance.json)

RNA002/RNA004 × Binary/Multi 共四套chemistry/task组合，各有纯三代与二代加三代融合两条主要路线；共享二代模型由25个81 nt序列checkpoint组成。每条路线可分别使用S1或S2策略。

人工示例用于检查输入格式和概率路径，不是准确率benchmark。精简包移除了自动生成的核验报告与重复样例，但保留9个代表性数据文件、完整CLI选项文档、测试代码和包完整性检查器。查看根目录[START_HERE.zh-CN.md](../../START_HERE.zh-CN.md)了解安装和模型下载。
