# 环境与版本证据

## 推理

原包推荐 Python 3.10/3.11，并通过 `constraints-models.txt` 固定 numpy 1.26.4、pandas 2.3.1、scipy 1.11.4、scikit-learn 1.3.0、xgboost 2.1.4、joblib 1.5.1、torch 2.7.1。安装命令见[入口](../START_HERE.zh-CN.md)。这些是原包约束，不应误写成本次实际测试环境。

本地曾使用超算已有 `nanoremix_env` 验证推理，观测到 Python 3.10.18、numpy 1.24.3、pandas 2.0.3、scipy 1.10.1、scikit-learn 1.3.0、xgboost 2.1.4、joblib 1.5.2、torch 2.2.2+cu118；真实模型 smoke 在 CPU 上运行。精简上传包不含自动生成的 smoke 报告。不同环境不能仅凭版本范围推断数值逐位相同；不要把 core CI 的 Python 3.12 支持当作所有旧模型依赖支持。

基础概率输入仅需 core numpy/pandas，不加载模型；特征推理另需 models extras，现场二代推理另需 second extras。安装指令会从软件源下载依赖，本包没有把 conda/venv、CUDA 或 wheel 仓库打进 ZIP。

## RNA002 上游

Guppy 6.3.8+d9e0f64（日志核实），RNA002 HAC 配置 `rna_r9.4.1_70bps_hac.cfg`；ont-fast5-api 4.1.1；Tombo 1.0（独立 Python 2.7 环境）；GraphMap 0.5.2；已有特征提取环境 Python 3.6.9、numpy 1.19.2、h5py 2.10.0。Tombo 和提取器不能因脚本同处目录就共用一个 Python。

## RNA004 上游

Dorado 1.3.0 与 `rna004_130bps_hac@v5.2.0`（现有 pipeline 配置），必须保存 moves。Remora 环境观测为 Python 3.9.20、ont-remora 3.3.0、pod5 0.3.15、pysam 0.22.1、numpy 1.24.4、pandas 2.3.3、h5py 3.11.0、torch 2.4.1。原脚本没有固定 samtools module 的历史版本，运行时需记录 `samtools --version`，本说明不补造版本。

RNA004 refiner 配置、9-mer level 表、region 起点定义、read/strand/参考中心和各统计量字段必须同步核对。`upstream_reference/SOURCE_MAPPING.json` 记录副本 SHA256；原路径和观测证据记录见[provenance](../docs/guide_20261104/provenance.json)。

本次没有重新执行上游软件安装、basecalling、resquiggle 或原始数据计算。外部工具及模型须按其来源与许可获取；随包保留了用于核对的脚本和 level 表，不能用它们替代对应工具安装。
