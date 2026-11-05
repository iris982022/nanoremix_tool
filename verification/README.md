# 分发包完整性检查

精简源码 ZIP 保留以下轻量验证材料：

- `verify_distribution.py`：只依赖Python标准库，核对包内 `SHA256SUMS.txt`、`MANIFEST.json`、缺失文件和模型文件计数；不会反序列化模型。
- `../tests/`：可安装依赖后运行的回归测试源代码。
- `../docs/guide_20261104/examples/`：用于说明输入格式的人工构造文件，不是生物学benchmark。

在解压后的 `NanoRemix/` 根目录运行：

```bash
python verification/verify_distribution.py
python -m pytest -q
```

第一条只检查文件完整性，不运行NanoRemix模型或验证预处理结果。第二条需要先安装项目依赖。精简包不包含自动生成的pytest日志、模型smoke报告和文档检查报告；这类运行产物不属于安装或使用工具所需文件。

本次未重新运行原始FAST5/POD5 basecalling或特征提取，也未测定新样本生物学准确率或证明论文结果可复现。实际模型推理还需要从GitHub Release下载对应权重，并按[模型说明](../docs/MODELS.md)执行检查。
