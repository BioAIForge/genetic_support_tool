# 遗传学支撑分析工具说明文档

## 1. 工具概述

本工具用于构建遗传学支撑分析流程，采用“Python 总控 + R 统计分析后端”的混合架构。

工具定位：

- 统一管理遗传分析方法的输入、参数、执行和结果输出
- 以 Python 作为流程编排入口
- 以 R 包或其他外部工具作为统计分析后端
- 便于在服务器环境中部署、测试和扩展

当前已实现的方法包括：

- `burden`：Burden 集合检验
- `skato`：SKAT-O 集合检验
- `haplotype`：局部 Haplotype 分析
- `quant_assoc`：连续表型关联分析
- `gwas_overlap`：与传统 GWAS hits 的重叠和补充关系分析

后续可继续扩展的方法包括：

- 其他基于 R、Python、Java 或命令行程序的分析模块

---

## 2. 工具架构

本工具采用分层设计：

### 2.1 Python 总控层

负责：

- 命令行入口
- 配置文件解析
- 输入文件检查
- 运行参数组织
- 调用 R 脚本或其他后端程序
- 结果文件整理

### 2.2 R 统计分析层

负责：

- 执行具体统计模型
- 输出标准结果表

### 2.3 外部工具扩展层

当前已预留对以下类型后端的接入能力：

- R 包
- Python 脚本
- Java 工具
- 命令行程序，如 `PLINK2`、`bedtools`

---

## 3. 项目结构

```text
genetic_support_tool/
  config/
    default.yaml
    gwas_overlap_demo.yaml
    regenie_example.yaml
    quant_assoc_plink2_example.yaml
    strong_signal_haplotype.yaml
    strong_signal_skat.yaml
    strong_signal_skato.yaml
  data/
    example/
    example_strong_signal/
    gwas_overlap/
  output/
  scripts/
    python/
      config_utils.py
      main.py
      run_burden.py
      run_gwas_overlap.py
      run_haplotype.py
      run_regenie.py
      run_quant_assoc.py
      run_skato.py
    r/
      burden.R
      haplotype.R
      quant_assoc.R
      skato.R
  requirements.txt
  README.md
```

## 3.1 快速开始

如果你是第一次把代码拉到服务器，建议按下面顺序测试：

1. 先完成最小环境安装
2. 先跑不依赖外部命令行工具的模块，确认 Python + R 链路可用
3. 再按需安装 `plink2`、`bedtools`、`regenie`，测试扩展模块

推荐的阅读顺序是：

- 本节先看整体流程
- 详细依赖见 `7. 依赖说明`
- 服务器安装见 `8. 服务器环境安装`
- 各模块命令见 `9. 运行方式`
- 如果希望用容器，直接看 `13. Docker 使用说明`

### 3.1.1 服务器最小可运行流程

以 Ubuntu / Debian 为例，拉取代码后可先执行：

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip r-base

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite","SKAT"), repos="https://cloud.r-project.org")'
```

完成后，建议先跑一个最小 smoke test：

```bash
python scripts/python/main.py burden --config config/default.yaml
```

如果该命令可以成功生成 `output/burden_skat/` 下的结果文件，说明：

- Python 入口可用
- YAML 配置解析可用
- `Rscript` 调用链路可用
- 基础的 TSV 输入输出流程可用

### 3.1.2 推荐测试顺序

推荐按依赖复杂度由低到高测试：

1. `burden --config config/default.yaml`
2. `skato --config config/default.yaml`
3. `haplotype --config config/strong_signal_haplotype.yaml`
4. `gwas-overlap --config config/gwas_overlap_demo.yaml`
5. `quant-assoc --config config/quant_assoc_plink2_example.yaml`
6. `burden` / `skato` 的 `engine=regenie` 配置

说明：

- 第 1 步主要验证 Python + R + TSV 基础链路
- 第 2、3 步会进一步验证 `SKAT`、`haplo.stats` 等 R 包依赖
- 第 4 步会验证 `bedtools`
- 第 5 步会验证 `plink2`
- 第 6 步最后验证 `regenie`

### 3.1.3 各模块开始前要注意什么

#### Burden / SKAT-O / Haplotype

这三个模块最适合先在服务器上跑通，因为它们直接使用仓库中的 TSV 示例数据。

常见前置条件：

- Python 依赖已安装
- `Rscript` 可用
- 对应 R 包已安装

#### Quantitative Association

`quant-assoc` 依赖 `plink2`。此外，当前示例配置默认读取 `PED/MAP` 示例前缀，因此第一次测试前建议先生成示例输入：

```bash
python scripts/python/make_plink2_example.py
python scripts/python/main.py quant-assoc --config config/quant_assoc_plink2_example.yaml
```

#### Burden / SKAT-O 的 regenie 场景

`regenie` 场景建议最后测试，因为它比其他模块多依赖：

- `regenie`
- `plink2`，如果你需要先把示例 `PED/MAP` 转成 `PGEN`
- annotation / set-list / mask-def / pred 等额外输入

示例测试方式：

```bash
plink2 --pedmap data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset
python scripts/python/main.py burden --config config/regenie_example.yaml
python scripts/python/main.py skato --config config/regenie_example.yaml
```

当前实现中，`engine=regenie` 会直接保留 regenie 原始输出文件，不再额外封装统一 TSV 结果表。

---

## 4. 功能模块说明

### 4.0 模块总览

| 模块名称 | 方法名 | 命令行入口 | 主要外部工具 / R 包 | 当前实现状态 | 主要输出 |
|---|---|---|---|---|---|
| `burden` | Burden 集合检验 | `python scripts/python/main.py burden --config ...` | `SKAT`、`regenie` | 已实现 | `output/burden_<engine>/` |
| `skato` | SKAT-O 集合检验 | `python scripts/python/main.py skato --config ...` | `SKAT`、`regenie` | 已实现 | `output/skato_<engine>/` |
| `haplotype` | 局部 Haplotype 分析 | `python scripts/python/main.py haplotype --config ...` | `haplo.stats` | 已实现 | `output/haplotype/` |
| `quant-assoc` | 连续表型关联分析 | `python scripts/python/main.py quant-assoc --config ...` | `PLINK2` | 已实现 | `output/quant_assoc_<engine>/` |
| `gwas-overlap` | GWAS hits 重叠分析 | `python scripts/python/main.py gwas-overlap --config ...` | `bedtools` | 已实现 | `output/gwas_overlap/` |

### 4.1 Burden 模块

方法名称：

- `burden`

用途：

- 对一个基因或区域中的多个变异做集合层面的 Burden 检验
- 判断总体变异负荷是否与目标表型相关

当前支持引擎：

#### `engine=skat`

说明：

- 标准公共包实现
- 调用 `SKAT` R 包
- 使用 `method = "Burden"`

适用场景：

- 正式的稀有变异集合检验
- 与 `SKAT-O` 做对照分析

#### `engine=regenie`

说明：

- 基于 `regenie` 的 step 2 set-based tests
- 适用于原生 `PGEN` 或 `BED` 格式输入
- 可用于 burden 风格的 mask-based rare variant testing
- 需要先完成 `regenie` step 1，并提供 step 2 所需的 `pred` 文件

输入文件：

- `engine=skat`
  - `geno_matrix.tsv`
  - `pheno.tsv`
  - `covar.tsv`
- `engine=regenie`
  - `PGEN/PVAR/PSAM` 或 `BED/BIM/FAM`
  - phenotype file
  - covariate file
  - annotation file
  - set-list file
  - mask-def file

输出文件：

- `engine=skat`
  - `output/burden_skat/burden_result.tsv`
  - `output/burden_skat/burden_scores.tsv`
  - `output/burden_skat/run_metadata.json`
- `engine=regenie`
  - `output/burden_regenie/regenie_burden_<phenotype>.regenie`
  - `output/burden_regenie/run_metadata.json`

结果解读重点：

- `engine=skat`
  - `burden_pvalue`
- `engine=regenie`
  - 直接查看 regenie 原始结果文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段
  - `run_metadata.json` 中的 `result_file` 会指向原始 `.regenie` 文件

### 4.2 SKAT-O 模块

方法名称：

- `skato`

用途：

- 对一个基因或区域中的多个变异做 `SKAT-O` 集合检验
- 在 Burden 与 SKAT 两种统计策略之间自适应组合

实现方式：

- 调用 `SKAT` R 包
- 使用 `method = "SKATO"`
- 可选新增 `engine=regenie`
  - 调用 `regenie --step 2`
  - 使用 `--vc-tests skato`

适用场景：

- 变异集合内部效应模式可能不一致
- 希望采用比 Burden 更稳健的集合检验方法

输入文件：

- `engine=skat`
  - `geno_matrix.tsv`
  - `pheno.tsv`
  - `covar.tsv`
- `engine=regenie`
  - `PGEN/PVAR/PSAM` 或 `BED/BIM/FAM`
  - phenotype file
  - covariate file
  - annotation file
  - set-list file
  - mask-def file

输出文件：

- `engine=skat`
  - `output/skato_skat/skato_result.tsv`
  - `output/skato_skat/run_metadata.json`
- `engine=regenie`
  - `output/skato_regenie/regenie_skato_<phenotype>.regenie`
  - `output/skato_regenie/run_metadata.json`

结果解读重点：

- `engine=skat`
  - `skato_pvalue`
- `engine=regenie`
  - 直接查看 regenie 原始结果文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段
  - `run_metadata.json` 中的 `result_file` 会指向原始 `.regenie` 文件

说明：

- `SKAT-O` 重点输出集合层面的显著性
- 当前版本不输出统一的 `beta` 和 `SE`
- `engine=regenie` 时不再额外封装统一结果表，保留 regenie 原始输出

### 4.3 局部 Haplotype 模块

方法名称：

- `haplotype`

用途：

- 在局部区域内根据多个位点组合构建 haplotype
- 评估主要 haplotype 与目标表型之间的关联

实现方式：

- 调用 `haplo.stats` R 包
- 将当前 0/1/2 剂量矩阵自动转换为双等位基因输入
- 使用 `haplo.em` 估计 haplotype 频率
- 使用 `haplo.score` 输出 haplotype 全局关联检验结果

输入文件：

- `geno_matrix.tsv`
- `pheno.tsv`
- `covar.tsv`，可选

输出文件：

- `output/haplotype/haplotype_result.tsv`
- `output/haplotype/haplotype_frequency.tsv`
- `output/haplotype/run_metadata.json`

结果解读重点：

- `top_haplotype`
- `top_haplotype_frequency`
- `haplotype_pvalue`

说明：

- 当前输出的 `haplotype_pvalue` 为基于 `haplo.stats` 的全局 haplotype 关联检验 P 值
- `haplotype_effect` 在当前版本中不再作为主要解释字段，可能为 `NA`

### 4.4 连续表型关联分析模块

方法名称：

- `quant-assoc`

用途：

- 对基因矩阵中的每个位点逐一进行连续表型关联分析
- 输出位点级别的 `beta`、`SE` 和 `P 值`

实现方式：

- 由 Python 将当前 `TSV genotype matrix` 转换为临时 `PED/MAP`
- 调用 `PLINK2 --glm`
- 解析 `PLINK2` 输出为统一结果表

输入文件：

- `geno_matrix.tsv`
- `pheno.tsv`
- `covar.tsv`

输出文件：

- `output/quant_assoc_plink2/quant_assoc_result.tsv`
- `output/quant_assoc_plink2/run_metadata.json`

结果解读重点：

- `variant_id`
- `beta`
- `se`
- `p_value`

说明：

- 当前模块已升级为 `PLINK2` 版实现
- 当前模块直接使用 `PLINK2` 原生输入格式
- 运行环境中必须提前安装 `plink2`
- 默认启用 `--covar-variance-standardize`，以避免协变量量纲差异导致的数值不稳定

### 4.5 GWAS Overlap 模块

方法名称：

- `gwas-overlap`

用途：

- 将本项目结果与参考 GWAS hits 进行坐标级比较
- 判断结果属于已知信号、邻近已知信号还是新增信号

实现方式：

- 读取项目结果文件
- 支持本地 GWAS 参考 TSV，或远程下载的 GWAS Catalog 参考文件
- 将项目结果和参考结果转换为 BED
- 调用 `bedtools intersect -wao`
- 根据窗口大小和交集结果分类为 `known_exact / near_known / novel`

输入文件：

- `project_hits.tsv`
- `gwas_reference.tsv`，或远程 GWAS Catalog 下载文件
- `config.yaml`

输出文件：

- `output/gwas_overlap/gwas_overlap_result.tsv`
- `output/gwas_overlap/run_metadata.json`

结果解读重点：

- `nearest_gwas_hit`
- `distance_bp`
- `category`

---

## 5. 输入文件说明

### 5.1 基因型文件 `geno_matrix.tsv`

要求：

- 第一列必须为 `sample_id`
- 后续每一列代表一个变异位点
- 位点值通常为 `0`、`1`、`2`

示例：

```text
sample_id	var1	var2	var3	var4
S1	0	1	0	1
S2	1	0	1	0
```

### 5.2 表型文件 `pheno.tsv`

要求：

- 必须包含 `sample_id`
- 必须包含配置文件中指定的表型列

示例：

```text
sample_id	phenotype
S1	5.2
S2	6.1
```

### 5.3 协变量文件 `covar.tsv`

要求：

- 必须包含 `sample_id`
- 可包含年龄、性别、PC、批次等协变量

示例：

```text
sample_id	age	sex	PC1	PC2
S1	38	0	-0.02	0.01
S2	41	1	-0.01	0.00
```

---

## 6. 配置文件说明

默认配置文件：

- [default.yaml](D:\Genos\生信工具分析\genetic_support_tool\config\default.yaml)
- [regenie_example.yaml](D:\Genos\生信工具分析\genetic_support_tool\config\regenie_example.yaml)
- [quant_assoc_plink2_example.yaml](D:\Genos\生信工具分析\genetic_support_tool\config\quant_assoc_plink2_example.yaml)

regenie 示例模板文件：

- [example_pheno.tsv](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_pheno.tsv)
- [example_covar.tsv](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_covar.tsv)
- [example_dataset.ped](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_dataset.ped)
- [example_dataset.map](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_dataset.map)
- [example.annotations](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example.annotations)
- [example.setlist](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example.setlist)
- [example.masks](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example.masks)

示例：

```yaml
project:
  name: demo_genetic_support

analysis:
  method: burden

input:
  genotype_matrix: data/example/geno_matrix.tsv
  phenotype_file: data/example/pheno.tsv
  covariate_file: data/example/covar.tsv

phenotype:
  name: phenotype
  type: continuous

burden:
  set_id: DEMO_GENE
  engine: skat
  covariates:
    - age
    - sex
    - PC1
    - PC2

skato:
  set_id: DEMO_GENE
  covariates:
    - age
    - sex
    - PC1
    - PC2

output:
  root_dir: output
```

关键参数说明：

- `analysis.method`：方法名称
- `phenotype.name`：表型列名
- `phenotype.type`：`continuous` 或 `binary`
- `burden.set_id`：Burden 分析对象名称
- `burden.engine`：`skat` 或 `regenie`
- `burden.regenie_bin`：`regenie` 可执行文件路径，默认 `regenie`
- `burden.regenie_aaf_bins`：set-based burden 的 AAF 分箱，默认 `0.01`
- `burden.regenie_build_mask`：mask 聚合方式，默认 `max`
- `burden.regenie_ignore_pred`：是否跳过 step 1 prediction 文件，适合最小 demo，默认 `false`
- `burden.regenie_bsize`：step 2 block size，默认 `200`
- `burden.covariates`：Burden 协变量列表
- `skato.set_id`：SKAT-O 分析对象名称
- `skato.engine`：`skat` 或 `regenie`
- `skato.regenie_bin`：`regenie` 可执行文件路径，默认 `regenie`
- `skato.regenie_aaf_bins`：VC test 的 AAF 分箱，默认 `0.01`
- `skato.regenie_build_mask`：mask 聚合方式，默认 `max`
- `skato.regenie_ignore_pred`：是否跳过 step 1 prediction 文件，适合最小 demo，默认 `false`
- `skato.regenie_bsize`：step 2 block size，默认 `200`
- `skato.covariates`：SKAT-O 协变量列表
- `regenie.genotype_format`：`pfile` 或 `bfile`
- `input.pred_file`：`regenie` step 1 输出的 prediction 文件；当 `regenie_ignore_pred=false` 时必填
- `quant_assoc.covar_variance_standardize`：是否对协变量做方差标准化，建议设为 `true`
- `output.root_dir`：输出根目录

---

## 7. 依赖说明

### 7.1 Python 依赖

当前 Python 依赖：

- `PyYAML`

安装命令：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 7.2 R 依赖

最低版本要求：

- 推荐 `R >= 4.3`
- 不建议使用 `R 4.1.x` 及更低版本

说明：

- 当前工具中的 `SKAT`、`haplo.stats` 等 R 包在较旧版本 R 环境下可能存在不可安装或兼容性不足的问题
- 如果需要稳定运行 `burden(engine=skat)`、`skato`、`haplotype` 等模块，建议在较新的 R 版本中部署
- 由于 `haplo.stats` 已从 CRAN 主仓库移除，需通过 archive 包方式安装

Burden / SKAT-O / Haplotype 的基础 R 包：

- `optparse`
- `jsonlite`

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite"), repos="https://cloud.r-project.org")'
```

### 7.3 方法对应 R 包依赖

#### Burden 标准版

适用配置：

- `burden.engine: skat`

需要：

- `optparse`
- `jsonlite`
- `SKAT`

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite","SKAT"), repos="https://cloud.r-project.org")'
```

#### SKAT-O

适用方法：

- `skato`

需要：

- `optparse`
- `jsonlite`
- `SKAT`

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite","SKAT"), repos="https://cloud.r-project.org")'
```

说明：

- 如果方法底层调用 R 包，则执行环境中必须预先安装对应 R 包
- 仅安装 Python 依赖不足以运行完整工具

#### Burden / SKAT-O regenie 后端

适用配置：

- `burden.engine: regenie`
- `skato.engine: regenie`

需要：

- `regenie`

安装说明：

- `regenie` 为独立命令行工具，不通过 R 包安装
- 运行环境中需预先安装 `regenie`，并保证可在 PATH 中调用
- 如果不在 PATH 中，可在配置文件中通过 `burden.regenie_bin` 或 `skato.regenie_bin` 指定绝对路径
- 当前仓库提供了一套 toy demo 文件，可先用来验证 `regenie` 调用链路
- demo 运行前仍需先将 [example_dataset.ped](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_dataset.ped) / [example_dataset.map](D:\Genos\生信工具分析\genetic_support_tool\data\regenie\example_dataset.map) 转成 `PGEN`
- 真实分析建议使用 step 1 生成的 `pred` 文件；demo 可通过 `regenie_ignore_pred=true` 跳过
- toy demo 的 `regenie_aaf_bins` 已放宽到 `1.0`，否则示例数据会因过于常见而生成空结果
- 当前 `burden.engine=regenie` 与 `skato.engine=regenie` 会直接保留 regenie 原始输出文件，不再额外封装统一 TSV 结果表

验证命令：

```bash
regenie --help
```

#### Haplotype

适用方法：

- `haplotype`

需要：

- `optparse`
- `haplo.stats`
- `arsenal`
- `rms`

安装命令：

```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'options(timeout=300); dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","arsenal","rms"), repos="https://cloud.r-project.org")'
mkdir -p ~/R/pkgsrc
wget -c --timeout=120 --tries=5 -O ~/R/pkgsrc/haplo.stats_1.9.8.3.tar.gz https://cran.r-project.org/src/contrib/Archive/haplo.stats/haplo.stats_1.9.8.3.tar.gz
R CMD INSTALL -l ~/R/library ~/R/pkgsrc/haplo.stats_1.9.8.3.tar.gz
```

如服务器无法使用 `wget`，可改用：

```bash
curl -L --retry 5 --connect-timeout 30 --max-time 300 -o ~/R/pkgsrc/haplo.stats_1.9.8.3.tar.gz https://cran.r-project.org/src/contrib/Archive/haplo.stats/haplo.stats_1.9.8.3.tar.gz
R CMD INSTALL -l ~/R/library ~/R/pkgsrc/haplo.stats_1.9.8.3.tar.gz
```

安装完成后可用以下命令验证：

```bash
Rscript -e 'library(haplo.stats); cat("haplo.stats loaded successfully\n")'
```

#### Quantitative Association

适用方法：

- `quant-assoc`

需要：

- `plink2`

安装命令：

```bash
plink2 --version
```

说明：

- `quant-assoc` 当前不依赖 R 包
- 需要在执行环境中预先安装 `plink2`，并保证命令可在 PATH 中调用
- 如果 `plink2` 不在 PATH 中，可在配置文件中通过 `quant_assoc.plink2_bin` 指定可执行文件路径

`plink2` 安装建议：

1. 从 PLINK2 官方下载 Linux 版本压缩包
2. 解压后将 `plink2` 可执行文件放到 PATH 目录，或在配置中指定绝对路径

示例：

```bash
mkdir -p ~/tools/plink2
cd ~/tools/plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20260310.zip
unzip plink2_linux_x86_64_20260310.zip
chmod +x plink2
~/tools/plink2/plink2 --version
```

如果希望全局可用，可加入 PATH：

```bash
echo 'export PATH="$HOME/tools/plink2:$PATH"' >> ~/.bashrc
source ~/.bashrc
plink2 --version
```

说明：

- 实际下载文件名可能会随官方版本更新而变化，请以官方页面为准
- 官方下载与输入格式说明见：
  - https://www.cog-genomics.org/plink/2.0/
  - https://www.cog-genomics.org/plink/2.0/input

#### GWAS Overlap

适用方法：

- `gwas-overlap`

需要：

- `bedtools`

版本要求：

- 推荐版本：`bedtools 2.31.x`
- 最低建议版本：`bedtools >= 2.30`

安装命令：

```bash
bedtools --version
```

说明：

- `gwas-overlap` 当前不依赖 R 包
- 需要在执行环境中预先安装 `bedtools`
- 如果 `bedtools` 不在 PATH 中，可在配置文件中通过 `gwas_overlap.bedtools_bin` 指定可执行文件路径
- 当前实现核心依赖命令为 `bedtools intersect`

`bedtools` 服务器安装建议：

Ubuntu / Debian：

```bash
sudo apt-get update
sudo apt-get install -y bedtools
bedtools --version
```

如果系统仓库版本过旧，建议手动安装新版本。示例：

```bash
mkdir -p ~/tools/bedtools
cd ~/tools/bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools.static
mv bedtools.static bedtools
chmod +x bedtools
~/tools/bedtools/bedtools --version
```

如需加入 PATH：

```bash
echo 'export PATH="$HOME/tools/bedtools:$PATH"' >> ~/.bashrc
source ~/.bashrc
bedtools --version
```

如果不加入 PATH，可在配置文件中指定：

```yaml
gwas_overlap:
  bedtools_bin: /home/zj/tools/bedtools/bedtools
```

---

## 8. 服务器环境安装

### 8.1 Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip r-base
```

然后进入项目目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(c("optparse","jsonlite"), repos="https://cloud.r-project.org")'
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages("SKAT", repos="https://cloud.r-project.org")'
```

由于默认 Burden 引擎已经切换为 `SKAT`，建议在基础环境安装阶段直接完成 `SKAT` 安装。

如果还需要运行 `quant-assoc`，还需要安装 `plink2` 并确保命令可用。

如果需要运行 `gwas-overlap`，还需要安装 `bedtools` 并确保命令可用。

### 8.2 CentOS / RHEL

```bash
sudo yum install -y python3 python3-pip R
```

然后执行与上面相同的 Python 和 R 包安装步骤。

### 8.3 环境变量持久化

建议将用户级 R 包路径写入 shell 配置文件：

```bash
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
source ~/.bashrc
```

如果使用 `zsh`，则写入 `~/.zshrc`。

---

## 9. 运行方式

### 9.1 Burden 默认示例

```bash
python scripts/python/main.py burden --config config/default.yaml
```

### 9.2 Burden 强信号测试

```bash
python scripts/python/main.py burden --config config/strong_signal_skat.yaml
```

regenie 模板：

```bash
plink2 --pedmap data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset
python scripts/python/main.py burden --config config/regenie_example.yaml
```

### 9.3 SKAT-O 默认示例

```bash
python scripts/python/main.py skato --config config/default.yaml
```

### 9.4 SKAT-O 强信号测试

```bash
python scripts/python/main.py skato --config config/strong_signal_skato.yaml
```

regenie 模板：

```bash
plink2 --pedmap data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset
python scripts/python/main.py skato --config config/regenie_example.yaml
```

### 9.5 局部 Haplotype 示例

```bash
python scripts/python/main.py haplotype --config config/strong_signal_haplotype.yaml
```

### 9.6 连续表型关联分析示例

```bash
python scripts/python/main.py quant-assoc --config config/quant_assoc_plink2_example.yaml
```

### 9.7 GWAS Overlap 示例

```bash
python scripts/python/main.py gwas-overlap --config config/gwas_overlap_demo.yaml
```

---

## 10. 输出文件说明

### 10.1 Burden 输出

目录：

- `output/burden_skat/`
- `output/burden_regenie/`

文件：

- `engine=skat`
  - `burden_result.tsv`
  - `burden_scores.tsv`
  - `run_metadata.json`
- `engine=regenie`
  - `regenie_burden_<phenotype>.regenie`
  - `run_metadata.json`

`burden_result.tsv` 主要字段：

- `set_id`
- `engine`
- `phenotype_name`
- `phenotype_type`
- `n_samples`
- `n_variants`
- `burden_beta_or_logodds`
- `burden_se`
- `burden_pvalue`
- `burden_model`

解读：

- `burden_pvalue` 是最核心结果
- `engine=skat` 时，`beta` 和 `SE` 可能为 `NA`，此时重点看集合层面的 `pvalue`
- `engine=regenie` 时，不再额外生成统一封装的 `burden_result.tsv`
- `engine=regenie` 时，结果来自 regenie set-based burden 检验，直接查看原始 `.regenie` 文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段
- `engine=regenie` 时，`run_metadata.json` 中的 `result_file` 会指向原始 `.regenie` 文件

`burden_scores.tsv` 主要字段：

- `sample_id`
- `burden_score`
- `phenotype`

用途：

- 检查样本级 burden score 是否合理
- 辅助理解 burden score 与表型的关系

### 10.2 SKAT-O 输出

目录：

- `output/skato_skat/`
- `output/skato_regenie/`

文件：

- `engine=skat`
  - `skato_result.tsv`
  - `run_metadata.json`
- `engine=regenie`
  - `regenie_skato_<phenotype>.regenie`
  - `run_metadata.json`

`skato_result.tsv` 主要字段：

- `set_id`
- `engine`
- `phenotype_name`
- `phenotype_type`
- `n_samples`
- `n_variants`
- `skato_pvalue`
- `skato_model`

解读：

- `skato_pvalue` 是核心结果
- `P` 越小，说明该变异集合与表型关联越明显
- `engine=regenie` 时，不再额外生成统一封装的 `skato_result.tsv`
- `engine=regenie` 时，直接查看原始 `.regenie` 文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段
- `engine=regenie` 时，`run_metadata.json` 中的 `result_file` 会指向原始 `.regenie` 文件

### 10.3 Haplotype 输出

目录：

- `output/haplotype/`

文件：

- `haplotype_result.tsv`
- `haplotype_frequency.tsv`
- `run_metadata.json`

`haplotype_result.tsv` 主要字段：

- `set_id`
- `top_haplotype`
- `top_haplotype_frequency`
- `n_samples`
- `n_variants`
- `haplotype_effect`
- `haplotype_pvalue`
- `haplotype_model`

解读：

- `top_haplotype` 表示频率最高的 haplotype
- `top_haplotype_frequency` 表示该 haplotype 的估计频率
- `haplotype_pvalue` 表示全局 haplotype 关联检验的显著性
- 当前 `haplotype_effect` 可能为 `NA`，因为现阶段重点输出全局检验结果

### 10.4 连续表型关联分析输出

目录：

- `output/quant_assoc_plink2/`

文件：

- `quant_assoc_result.tsv`
- `run_metadata.json`

`quant_assoc_result.tsv` 主要字段：

- `variant_id`
- `chr`
- `pos`
- `ref`
- `alt`
- `n_samples`
- `beta`
- `se`
- `p_value`

### 10.5 GWAS Overlap 输出

目录：

- `output/gwas_overlap/`

文件：

- `gwas_overlap_result.tsv`
- `run_metadata.json`

`gwas_overlap_result.tsv` 主要字段：

- `project_id`
- `chr`
- `pos`
- `nearest_gwas_hit`
- `nearest_gwas_trait`
- `distance_bp`
- `category`

`category` 的常见取值：

- `known_exact`
- `near_known`
- `novel`

---

## 11. 使用说明

建议使用顺序：

1. 先安装 Python 依赖
2. 再安装所选方法需要的 R 包
3. 根据方法准备输入文件
   - `burden / skato / haplotype`：准备 `geno_matrix.tsv`、`pheno.tsv`、`covar.tsv`
   - `quant-assoc`：准备 `PLINK2` 原生输入文件和 `pheno.tsv`、`covar.tsv`
   - `gwas-overlap`：准备项目结果表和 GWAS 参考结果表
4. 编辑配置文件
5. 执行对应命令
6. 查看 `output/` 下结果文件

建议先使用强信号测试数据验证：

- 工具是否安装成功
- 方法是否能跑通
- 输出结果是否符合预期

---

## 12. 说明

当前版本为可运行的分析工具原型，已具备：

- 统一 CLI 入口
- 标准化输入输出
- Burden 方法
- SKAT-O 方法
- 局部 Haplotype 方法
- 连续表型关联分析方法
- GWAS Overlap 方法
- 示例数据
- 服务器部署说明

后续可在当前框架上继续扩展更多遗传学支撑分析模块，并在现有 5 个模块基础上增强输入适配、批量分析、显著结果汇总和可视化输出能力。

---

## 13. Docker 使用说明

本项目支持封装为 Docker 镜像运行，以统一环境依赖并简化部署。

当前 Docker 镜像中包含：

- 基于 `rocker/r2u` 的 R 运行环境
- Python 运行环境
- R 运行环境
- `SKAT`
- `haplo.stats`
- `plink2`
- `bedtools`
- 项目代码、默认配置文件和示例数据

说明：

- 当前 Docker 镜像尚未内置 `regenie`
- 如需使用 `burden.engine=regenie` 或 `skato.engine=regenie`，请在宿主环境安装 `regenie`，或基于当前镜像再构建自定义镜像

说明：


### 13.1 构建镜像

在项目根目录执行：

```bash
docker build -t genetic-support-tool .
```

### 13.2 基本运行方式

说明：

- 镜像内已包含项目代码、默认配置和示例数据
- 在其他环境中验证镜像时，不需要再额外拷贝整套项目代码
- Docker 环境中默认输出目录通过环境变量重定向到容器内 `/work/output`
- 如需保留输出结果，建议挂载 `/work`

示例：运行 `burden`

```bash
docker run --rm -it \
  -v $(pwd)/docker-output:/work \
  genetic-support-tool \
  burden --config config/default.yaml
```

示例：运行 `skato`

```bash
docker run --rm -it \
  -v $(pwd)/docker-output:/work \
  genetic-support-tool \
  skato --config config/default.yaml
```

示例：运行 `haplotype`

```bash
docker run --rm -it \
  -v $(pwd)/docker-output:/work \
  genetic-support-tool \
  haplotype --config config/strong_signal_haplotype.yaml
```

示例：运行 `quant-assoc`

```bash
docker run --rm -it \
  -v $(pwd)/docker-output:/work \
  genetic-support-tool \
  quant-assoc --config config/quant_assoc_plink2_example.yaml
```

示例：运行 `gwas-overlap`

```bash
docker run --rm -it \
  -v $(pwd)/docker-output:/work \
  genetic-support-tool \
  gwas-overlap --config config/gwas_overlap_demo.yaml
```

### 13.3 说明

- 镜像默认入口为：
  - `python3 scripts/python/main.py`
- 因此在 `docker run` 中只需要追加方法名和参数
- 如仅验证镜像自带示例，无需挂载项目代码目录
- 如需持久化输出结果，建议挂载宿主机目录到 `/work`
- 如需使用自定义配置或自定义输入数据，可额外挂载数据目录并在配置文件中使用绝对路径
