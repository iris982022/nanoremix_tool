# 人工格式示例

所有序列、信号值和概率均为人工构造，不来自真实实验，不能评价生物学性能。本精简包保留下列10个数据文件：

- `RNA002_prepared_with81.csv`、`RNA004_prepared_with81.csv`：分别展示9 nt/20值和17 nt/36值prepared输入，并附带融合使用的81 nt序列。
- `RNA002_binary_both_probabilities.csv`、`RNA004_binary_both_probabilities.csv`：Binary三代pair及完整二代概率输入。
- `RNA002_both_probabilities.csv`、`RNA004_both_probabilities.csv`：Multi原生三代类别轴及完整二代概率输入。
- `RNA002_step4.tsv`：RNA002上游Step4格式，可生成9 nt候选。
- `RNA004_step6_positive.csv`：RNA004正链Step6示例。
- `RNA004_step6_negative_limitation.csv`：复现旧prepare对当前RNA方向负链Step6的限制；不代表负链已获全流程支持。
- `contexts.fa`：与prepared元数据中的reference_name和0基center_pos0=40对应的人工参考序列。

概率示例无需模型；特征推理或81 nt现场推理需要安装对应依赖并下载Release权重。具体命令、输入合同、边界条件和公式见上一级[完整中文手册](../README_版本输入预处理参数与计算.md)。
