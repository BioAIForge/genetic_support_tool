# 遗传学支撑分析工具说明文档

## 1. 工具概述

本工具用于构建遗传学支撑分析流程，采用"Python 总控 + R 统计分析后端"的混合架构。

工具定位：
- 统一管理遗传分析方法的输入、参数、执行和结果输出
- 以 Python 作为流程编排入口
- 以 R 包或其他外部工具作为统计分析后端
- 便于在服务器环境中部署、测试和扩展

当前已实现的方法包括：
- `burden`：Burden 集合检验
- `skato`：SKAT-O 集合检验
- `quant_assoc`：连续表型关联分析
- `gwas_gene_catalog`：基因级 GWAS Catalog 证据回填

---

## 2. 工具架构

本工具采用分层设计：

### 2.1 Python 总控层
负责命令行入口、配置文件解析、输入文件检查、运行参数组织、调用 R 脚本或其他后端程序、结果文件整理。

### 2.2 R 统计分析层
负责执行具体统计模型，输出标准结果表。

### 2.3 外部工具扩展层
当前已预留对以下类型后端的接入能力：
- R 包
- Python 脚本
- Java 工具
- 命令行程序，如 `PLINK2`、`regenie`

---

## 3. 项目结构

```
genetic_support_tool/
  config/           # 配置文件
  data/             # 示例数据
  output/           # 输出目录
  scripts/
    python/        # Python 入口脚本
    r/              # R 统计脚本
  docs/             # 详细文档
  requirements.txt  # Python 依赖
```

---

## 4. 快速开始

### 4.1 环境安装

```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip r-base

# 创建虚拟环境并安装依赖
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 安装 R 包
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite","SKAT"), repos="https://cloud.r-project.org")'
```

详细依赖说明见 [docs/dependencies.md](docs/dependencies.md)。

### 4.2 测试运行

推荐测试顺序（按依赖复杂度由低到高）：

```bash
# 1. 测试基础链路
python scripts/python/main.py burden --config config/default.yaml

# 2. 测试 SKAT-O
python scripts/python/main.py skato --config config/default.yaml

# 3. 测试 GWAS Catalog 本地 TSV
python scripts/python/main.py gwas-gene-catalog --config config/gwas_gene_catalog_demo.yaml

# 4. 测试 GWAS Catalog API
python scripts/python/main.py gwas-gene-catalog --config config/gwas_gene_catalog_api_demo.yaml

# 5. 测试连续表型关联分析（需先安装 plink2）
python scripts/python/make_plink2_example.py
python scripts/python/main.py quant-assoc --config config/quant_assoc_plink2_example.yaml

# 6. 测试 regenie（需先安装 regenie）
plink2 --pedmap data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset
python scripts/python/main.py burden --config config/regenie_example.yaml
```

### 4.3 Docker 运行

```bash
# 构建镜像
docker build -t genetic-support-tool .

# 运行示例
docker run --rm -it \
  -v $(pwd)/output:/work \
  genetic-support-tool \
  burden --config config/default.yaml
```

Docker 详情见 [docs/dependencies.md](docs/dependencies.md)。

---

## 5. 功能模块总览

| 模块名称 | 方法名 | 命令行入口 | 主要外部工具 / R 包 | 当前实现状态 | 主要输出 |
|---|---|---|---|---|---|
| `burden` | Burden 集合检验 | `python scripts/python/main.py burden --config ...` | `SKAT`、`regenie` | 已实现 | `output/burden_<engine>/` |
| `skato` | SKAT-O 集合检验 | `python scripts/python/main.py skato --config ...` | `SKAT`、`regenie` | 已实现 | `output/skato_<engine>/` |
| `quant-assoc` | 连续表型关联分析 | `python scripts/python/main.py quant-assoc --config ...` | `PLINK2` | 已实现 | `output/quant_assoc_<engine>/` |
| `gwas-gene-catalog` | 基因级 GWAS Catalog 证据回填 | `python scripts/python/main.py gwas-gene-catalog --config ...` | GWAS Catalog TSV / 官方 REST API | 已实现 | `output/gwas_gene_catalog/` |

各模块的详细说明、输入输出格式见 [docs/io_files.md](docs/io_files.md)。

---

## 6. 配置文件

配置文件位于 `config/` 目录：

| 配置文件 | 用途 |
|----------|------|
| `default.yaml` | 默认配置 |
| `strong_signal_skat.yaml` | Burden 强信号测试 |
| `strong_signal_skato.yaml` | SKAT-O 强信号测试 |
| `regenie_example.yaml` | regenie 示例 |
| `quant_assoc_plink2_example.yaml` | PLINK2 示例 |
| `gwas_gene_catalog_demo.yaml` | GWAS Catalog 本地 TSV 示例 |
| `gwas_gene_catalog_api_demo.yaml` | GWAS Catalog API 示例 |
| `gwas_gene_catalog_official.yaml` | 官方完整 TSV + API 索引 |

完整参数说明见 [docs/config_reference.md](docs/config_reference.md)。

---

## 7. 运行方式

### 7.1 Burden

```bash
python scripts/python/main.py burden --config config/default.yaml
```

### 7.2 SKAT-O

```bash
python scripts/python/main.py skato --config config/default.yaml
```

### 7.3 连续表型关联分析

```bash
python scripts/python/main.py quant-assoc --config config/quant_assoc_plink2_example.yaml
```

### 7.4 GWAS Catalog 基因注释

```bash
# 本地 TSV 模式
python scripts/python/main.py gwas-gene-catalog --config config/gwas_gene_catalog_demo.yaml

# API 模式
python scripts/python/main.py gwas-gene-catalog --config config/gwas_gene_catalog_api_demo.yaml
```

---

## 8. 说明

当前版本为可运行的分析工具原型，已具备：
- 统一 CLI 入口
- 标准化输入输出
- Burden 方法
- SKAT-O 方法
- 连续表型关联分析方法
- GWAS Catalog 基因注释方法
- 示例数据
- 服务器部署说明

后续可在当前框架上继续扩展更多遗传学支撑分析模块。

---

## 9. 参考文档

- [输入输出文件说明](docs/io_files.md)
- [配置文件参数参考](docs/config_reference.md)
- [依赖安装说明](docs/dependencies.md)