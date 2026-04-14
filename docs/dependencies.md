# 依赖安装说明

本文档详细说明各模块所需的依赖及其安装方法。

## 目录

- [Python 依赖](#python-依赖)
- [R 依赖](#r-依赖)
- [外部命令行工具](#外部命令行工具)
- [系统要求](#系统要求)
- [Docker 环境](#docker-环境)

---

## Python 依赖

### 必需依赖

| 包名 | 说明 |
|------|------|
| `PyYAML` | 配置文件解析 |

### 安装方法

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## R 依赖

### R 版本要求

- 推荐 `R >= 4.3`
- 不建议使用 `R 4.1.x` 及更低版本

> 说明：`SKAT` 等 R 包在较旧版本 R 环境下可能存在兼容性问题。

### 基础 R 包

所有方法都需要的基础包：

| 包名 | 说明 |
|------|------|
| `optparse` | 命令行参数解析 |
| `jsonlite` | JSON 处理 |

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite"), repos="https://cloud.r-project.org")'
```

### 方法对应 R 包

#### Burden / SKAT-O (engine=skat)

需要额外安装：

| 包名 | 说明 |
|------|------|
| `SKAT` | 集合检验核心包 |

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages("SKAT", repos="https://cloud.r-project.org")'
```

#### 连续表型关联分析

- 当前不依赖 R 包

#### GWAS Catalog 基因注释

- 默认不依赖额外命令行工具
- 如需远程下载 TSV 或调用 API，需保证网络访问

---

## 外部命令行工具

### regenie

| 项目 | 说明 |
|------|------|
| 用途 | Burden / SKAT-O 的替代引擎 |
| 安装方式 | 独立命令行工具，不通过 R 包安装 |
| 位置 | 需预先安装并在 PATH 中，或通过配置文件指定路径 |

安装步骤（Linux）：

```bash
mkdir -p ~/tools/regenie
cd ~/tools/regenie
wget https://github.com/rgcgithub/regenie/releases/download/v3.1.2/regenie_v3.1.2_linux_x86_64.zip
unzip regenie_v3.1.2_linux_x86_64.zip
chmod +x regenie
./regenie --help
```

配置文件中指定路径：

```yaml
burden:
  engine: regenie
  regenie_bin: ~/tools/regenie/regenie
```

### plink2

| 项目 | 说明 |
|------|------|
| 用途 | 连续表型关联分析 |
| 安装方式 | 独立命令行工具 |
| 版本建议 | 近期 alpha 构建 |

安装步骤（Linux）：

```bash
mkdir -p ~/tools/plink2
cd ~/tools/plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20260310.zip
unzip plink2_linux_x86_64_20260310.zip
chmod +x plink2
./plink2 --version
```

配置文件中指定路径：

```yaml
quant_assoc:
  plink2_bin: ~/tools/plink2/plink2
```

---

## 系统要求

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip r-base
```

### CentOS / RHEL

```bash
sudo yum install -y python3 python3-pip R
```

### 环境变量持久化

```bash
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
source ~/.bashrc
```

如果使用 `zsh`，写入 `~/.zshrc`。

---

## Docker 环境

项目提供 Docker 镜像，内置所有依赖。

### 镜像包含

- 基于 `rocker/r2u` 的 R 运行环境
- Python 运行环境
- R 包：`optparse`、`jsonlite`、`SKAT`
- 命令行工具：`plink2`、`regenie`
- 项目代码、配置和示例数据

### 构建镜像

```bash
docker build -t genetic-support-tool .
```

### 运行示例

```bash
# 运行 burden
docker run --rm -it \
  -v $(pwd)/output:/work \
  genetic-support-tool \
  burden --config config/default.yaml

# 运行 skato
docker run --rm -it \
  -v $(pwd)/output:/work \
  genetic-support-tool \
  skato --config config/default.yaml

# 运行 quant-assoc
docker run --rm -it \
  -v $(pwd)/output:/work \
  genetic-support-tool \
  quant-assoc --config config/quant_assoc_plink2_example.yaml
```

### 自动构建

仓库已配置 GitHub Actions 自动构建 Docker 镜像：

- 触发条件：push 到 main、PR 合并、发布 tag
- 镜像推送至 GitHub Container Registry (`ghcr.io`)