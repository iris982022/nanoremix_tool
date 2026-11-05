# 使用边界与已知限制

本次保持 v0.3.0 推理源码和所有模型不变。以下内容决定是否可以将自己的输入交给模型，不能用格式校验通过替代实验流程核验。

1. **RNA004 负链转换**：旧 `nanoremix prepare --format step6` 按参考起点升序拼接重叠 9-mer，而现行提取器将负链序列/数值转成 RNA 方向；两者直接衔接可能丢失负链候选。`docs/guide_20261104/examples/RNA004_step6_negative_limitation.csv` 与核验脚本记录了这一现象。本次没有改变转换器，超算 `FINAL_TOOL/prepare_features.py` 也不能仅凭文件较新就视为已解决。
2. **未知样本的 RNA004 适配仍需落实**：随包 `prepare_model_input.py`、`add_splice_aware_context.py` 是已有 catalog 验证链源码。前者使用参考 catalog 上下文和类别，后者还含 catalog/site_id/真值相关逻辑，Binary 结果按 truth 过滤。不能把这些条件当作未知样本的先验。通用适配需去除真值依赖，核对坐标、链、信号槽位和连续上下文。跨剪接的 81 nt 需按转录本或 read 的 CIGAR 路径跳过 N，不能直接取基因组左右 40 bp。
3. **RNA002 坐标**：Step4 的 event index 不是转录本/基因组坐标。利用同步 sidecar 显式连接后才能在同版本参考中取 81 nt；该提取器不直接产生完整事件 FASTA。
4. **候选与输入**：每行是 read×中心观察。先筛掉邻近事件再拼上下文会改变输入；未经校验的新 basecaller/normalization、不同数值定义或任意标准化并不保证与模型训练输入等价。RNA004 region 的 `[start,end)` 表示 9-mer 左起点集合，单碱基 site BED 不能直接等同。
5. **权重覆盖**：完整包只有 Binary `full` 特征权重。CLI 虽提供消融 feature-set 选项，对应权重需另行准备。支持类别以 `nanoremix catalog` 和精确 JSON 为准；别名不增加类别，三代独有类别没有独立二代证据。
6. **阈值与输出**：置信度是当前所选类别中的条件概率。排除真正类别后其质量可能分配给已选类别；按 modified 输出也不等于显著位点。工具没有自动位点聚合、修饰率、覆盖度门槛、motif 过滤、S1/S2 共识或论文特定 cohort 筛选。
7. **规模和环境**：`predict --chunk-size` 不控制 `prepare`、FASTA 和独立二代概率表的全部内存。按完整 read 分片，保证边界与身份唯一性。上游 CUDA basecaller、Tombo/Remora 和推理使用独立环境；本包不是离线软件镜像。
8. **验证范围**：包完整性、原包一致性、回归测试、人工格式/公式例子与实际权重小样本推理各自只覆盖自己的范围。没有在本次重跑原始 FAST5/POD5，也没有验证生物学准确率或论文数值复现。详见[核验说明](../verification/README.md)。

完整输入与算法说明见[中文手册](guide_20261104/README_版本输入预处理参数与计算.md)，上游阶段见[命令模板](guide_20261104/PIPELINE_命令与输入衔接.md)。
