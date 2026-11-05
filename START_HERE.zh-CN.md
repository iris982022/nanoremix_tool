# NanoRemix：请从主页开始

工具仍为v0.3.0；2026-11-05更新主页与分发清单，没有更换推理代码或模型。

功能、四套配置、支持修饰、输入、上游流程、安装、运行、全部41个predict/validate长选项、筛选计算和输出已统一写进根目录 [README.md](README.md)。放在仓库根目录后，GitHub会自动显示在主页。

- [完整技术手册](docs/guide_20261104/README_版本输入预处理参数与计算.md)：字段、版本依据、预处理和公式。
- [上游命令模板](docs/guide_20261104/PIPELINE_命令与输入衔接.md)：RNA002/RNA004分阶段流程。
- [更新已有GitHub仓库](docs/UPLOAD.zh-CN.md)：99个源码文件、主页README和五个Release附件。
- [已知限制](docs/KNOWN_LIMITATIONS.zh-CN.md)：包括RNA004旧Step6负链问题。

新源码包为 `NanoRemix_v0.3.0_GitHub_READY_20261105.zip`。解压后把NanoRemix内部内容放在仓库根目录，保留src/docs等层级；不要平铺。源码不含权重，模型继续使用20261104的五个ZIP。
