# 输入输出文件说明

本文档说明各模块的输入文件格式和输出文件结构。

## 目录

- [输入文件](#输入文件)
  - [基因型文件](#基因型文件-genomatrixtsv)
  - [表型文件](#表型文件-phenotsv)
  - [协变量文件](#协变量文件-covartsv)
  - [其他输入文件](#其他输入文件)
- [输出文件](#输出文件)
  - [Burden 模块输出](#burden-模块输出)
  - [SKAT-O 模块输出](#skat-o-模块输出)
  - [连续表型关联分析输出](#连续表型关联分析输出)
  - [GWAS Catalog 基因注释输出](#gwas-catalog-基因注释输出)

---

## 输入文件

### 基因型文件 `geno_matrix.tsv`

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

### 表型文件 `pheno.tsv`

要求：
- 必须包含 `sample_id`
- 必须包含配置文件中指定的表型列

示例：

```text
sample_id	phenotype
S1	5.2
S2	6.1
```

### 协变量文件 `covar.tsv`

要求：
- 必须包含 `sample_id`
- 可包含年龄、性别、PC、批次等协变量

示例：

```text
sample_id	age	sex	PC1	PC2
S1	38	0	-0.02	0.01
S2	41	1	-0.01	0.00
```

### 其他输入文件

#### regenie 模式

使用 `engine=regenie` 时，需要以下文件：

| 文件类型 | 说明 |
|----------|------|
| `PGEN/PVAR/PSAM` 或 `BED/BIM/FAM` | 基因型数据 |
| phenotype file | 表型文件 |
| covariate file | 协变量文件 |
| annotation file | 变异注释文件 |
| set-list file | 集合定义文件 |
| mask-def file | mask 定义文件 |

示例文件位于 `data/regenie/` 目录：

- [example_pheno.tsv](data/regenie/example_pheno.tsv)
- [example_covar.tsv](data/regenie/example_covar.tsv)
- [example_dataset.ped](data/regenie/example_dataset.ped)
- [example_dataset.map](data/regenie/example_dataset.map)
- [example.annotations](data/regenie/example.annotations)
- [example.setlist](data/regenie/example.setlist)
- [example.masks](data/regenie/example.masks)

#### GWAS Catalog 基因注释

| 文件类型 | 说明 |
|----------|------|
| `project_result_file` | 本地基因级结果表，一行一个基因 |
| `gwas_reference_file` | 本地 GWAS Catalog 参考 TSV |
| `gwas_catalog_tsv_url` | 远程 GWAS Catalog TSV 下载地址 |

---

## 输出文件

### Burden 模块输出

目录：
- `output/burden_skat/`
- `output/burden_regenie/`

#### engine=skat

| 文件 | 说明 |
|------|------|
| `burden_result.tsv` | 结果主表 |
| `burden_scores.tsv` | 样本级 Burden Score |
| `run_metadata.json` | 运行元数据 |

**burden_result.tsv 主要字段：**

| 字段 | 说明 |
|------|------|
| `set_id` | 分析对象名称 |
| `engine` | 使用的引擎 |
| `phenotype_name` | 表型名称 |
| `phenotype_type` | 表型类型 |
| `n_samples` | 样本数 |
| `n_variants` | 变异数 |
| `burden_beta_or_logodds` | 效应值 |
| `burden_se` | 标准误 |
| `burden_pvalue` | P 值（核心结果） |
| `burden_model` | 模型类型 |

**burden_scores.tsv 主要字段：**

| 字段 | 说明 |
|------|------|
| `sample_id` | 样本 ID |
| `burden_score` | Burden Score |
| `phenotype` | 表型值 |

#### engine=regenie

| 文件 | 说明 |
|------|------|
| `regenie_burden_<phenotype>.regenie` | regenie 原始结果 |
| `run_metadata.json` | 运行元数据 |

解读：直接查看 `.regenie` 文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段。

---

### SKAT-O 模块输出

目录：
- `output/skato_skat/`
- `output/skato_regenie/`

#### engine=skat

| 文件 | 说明 |
|------|------|
| `skato_result.tsv` | 结果主表 |
| `run_metadata.json` | 运行元数据 |

**skato_result.tsv 主要字段：**

| 字段 | 说明 |
|------|------|
| `set_id` | 分析对象名称 |
| `engine` | 使用的引擎 |
| `phenotype_name` | 表型名称 |
| `phenotype_type` | 表型类型 |
| `n_samples` | 样本数 |
| `n_variants` | 变异数 |
| `skato_pvalue` | P 值（核心结果） |
| `skato_model` | 模型类型 |

#### engine=regenie

| 文件 | 说明 |
|------|------|
| `regenie_skato_<phenotype>.regenie` | regenie 原始结果 |
| `run_metadata.json` | 运行元数据 |

解读：直接查看 `.regenie` 文件中的 `P` / `LOG10P`、`TEST`、`ID`、`N`、`NBURDEN` 等字段。

---

### 连续表型关联分析输出

目录：
- `output/quant_assoc_plink2/`

| 文件 | 说明 |
|------|------|
| `quant_assoc_result.tsv` | 结果主表 |
| `run_metadata.json` | 运行元数据 |

**quant_assoc_result.tsv 主要字段：**

| 字段 | 说明 |
|------|------|
| `variant_id` | 变异 ID |
| `chr` | 染色体 |
| `pos` | 位置 |
| `ref` | 参考等位基因 |
| `alt` | 替代等位基因 |
| `n_samples` | 样本数 |
| `beta` | 效应值 |
| `se` | 标准误 |
| `p_value` | P 值 |

---

### GWAS Catalog 基因注释输出

目录：
- `output/gwas_gene_catalog/`

| 文件 | 说明 |
|------|------|
| `project_results_with_gwas_catalog.tsv` | 带 GWAS Catalog 注释的主表 |
| `gwas_catalog_gene_hits.tsv` | 命中明细表 |
| `run_metadata.json` | 运行元数据 |

#### project_results_with_gwas_catalog.tsv 新增字段

| 字段 | 说明 |
|------|------|
| `gwas_catalog_found` | 是否命中 GWAS Catalog |
| `gwas_catalog_result_count` | 命中记录数 |
| `gwas_catalog_traits` | 去重后的标准化 trait |
| `gwas_catalog_reported_traits` | 去重后的 reported trait |
| `gwas_catalog_mapped_genes` | 去重后的 mapped genes |
| `gwas_catalog_top_variant` | 最显著记录的变异标识 |
| `gwas_catalog_top_pvalue` | 最显著记录的 P 值 |
| `gwas_catalog_top_study` | 最显著记录的研究信息 |
| `gwas_catalog_study_accessions` | 所有命中的 study accession |

#### gwas_catalog_gene_hits.tsv 主要字段

| 字段 | 说明 |
|------|------|
| `project_row_index` | 对应主表中的行号 |
| `gene` | 基因 symbol |
| `gwas_catalog_trait` | 标准化 trait |
| `gwas_catalog_reported_trait` | 原始 reported trait |
| `gwas_catalog_mapped_genes` | Mapped genes |
| `gwas_catalog_reported_genes` | Reported genes |
| `gwas_catalog_variant` | 变异标识 |
| `gwas_catalog_pvalue` | P 值 |
| `gwas_catalog_or_value` | OR 值 |
| `gwas_catalog_beta` | Beta 值 |
| `gwas_catalog_ci` | 置信区间 |
| `gwas_catalog_study_accession` | Study accession |
| `gwas_catalog_first_author` | 第一作者 |
| `gwas_catalog_pubmed_id` | PubMed ID |