# 更新GitHub源码与主页

目标仓库：https://github.com/iris982022/nanoremix_tool 。工具仍是v0.3.0；本次只更新文档和分发清单。

## 1. 上传哪些内容

本次源码包为 `NanoRemix_v0.3.0_GitHub_READY_20261105.zip`，内部共99个文件。解压得到 `NanoRemix/`，把**文件夹内部的内容**放到仓库根目录，保留src、docs等目录；不要平铺所有子文件，也不要再套一层NanoRemix目录。

正确文件树包括：

```text
README.md
pyproject.toml
constraints-models.txt
src/nanoremix/cli.py
docs/...
examples/...
tests/...
upstream_reference/...
verification/...
.github/workflows/tests.yml
```

根目录README已含全部功能、四套配置、支持修饰、输入、上游流程、安装、运行、全部41个参数、筛选计算和限制。GitHub会自动将它显示在主页；另外提供的README副本与ZIP中相同，只需在仓库保留一份README.md。

GitHub不会自动展开ZIP：仅上传ZIP时，不能直接浏览里面的源文件，也不会自动显示ZIP内README。只更新README可先改善主页，但不能修复源码结构或补齐权重。

## 2. 推荐用GitHub Desktop更新已有仓库

此方法不需要在超算终端输入token。

1. 自己的电脑安装GitHub Desktop，通过浏览器登录GitHub账户。
2. 选择 **File → Clone repository → URL**，填 `https://github.com/iris982022/nanoremix_tool.git` 并克隆到本地。
3. 备份这个本地克隆目录。删除其中之前错位的源码文件，保留隐藏的 `.git/` 目录；不要在其他数据目录执行。
4. 解压新源码ZIP，把 `NanoRemix/` **内部全部内容**复制到本地仓库目录，包含 `.gitignore` 和 `.github/`；macOS/Linux需显示隐藏文件。
5. 确认根目录是README、pyproject、src、docs等，不是99个平铺文件，也不是 `nanoremix_tool/NanoRemix/src`。
6. 在Desktop查看删除和新增清单，Summary填写 `Restore package structure and document NanoRemix usage`，点击 **Commit to main → Push origin**。
7. 打开网页确认中文主页显示，并能进入 `src/nanoremix/` 和 `docs/`。

Git没有网页上传器的100文件批次限制。重整结构会同时删除旧路径、添加正确路径，Desktop变更数量可能超过100，不代表新源码包文件数超过99。

## 3. 如果坚持网页上传

上传解压后项目内容并保留子目录。提交前必须在预览中看到 `src/nanoremix/cli.py` 等正确路径；如果只有cli.py等文件名在根目录，先取消上传。网页批次若因删除与新增数量或浏览器行为无法完成，使用上面的Desktop方式。

只先改主页：打开根目录 **README.md → 铅笔 Edit this file**，用本次提供的README全文替换，然后 **Commit changes**。完整功能运行仍需修复目录并补齐模型下载。

## 4. 权重继续作为Release附件

进入 **Releases → Draft a new release**，选择正确main提交，创建或选择 `v0.3.0` 标签，附加下面五个原始ZIP，不要将解压权重放入源码仓库。

| 附件 | 用途 |
|---|---|
| `NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip` | RNA002 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip` | RNA002 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip` | RNA004 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip` | RNA004 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_SecondGen_20261104.zip` | 共享25个二代checkpoint及网络定义 |

权重不变，继续使用20261104的原始附件及对应SHA256文本。纯三代下载对应chemistry/task包，现场序列融合再加SecondGen，已有二代概率无需SecondGen。下载者解压到项目根目录形成models。发布包没有mean/std/length消融权重；RNA004旧Step6负链限制仍存在，详见[已知限制](KNOWN_LIMITATIONS.zh-CN.md)。
