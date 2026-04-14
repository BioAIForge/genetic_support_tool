# 配置文件参数参考

本文档提供所有配置参数的完整说明。

## 目录

- [配置文件位置](#配置文件位置)
- [通用参数](#通用参数)
- [Burden 参数](#burden-参数)
- [SKAT-O 参数](#skat-o-参数)
- [连续表型关联分析参数](#连续表型关联分析参数)
- [GWAS Catalog 基因注释参数](#gwas-catalog-基因注释参数)

---

## 配置文件位置

项目提供以下配置文件：

| 配置文件 | 用途 |
|----------|------|
| `config/default.yaml` | 默认配置 |
| `config/strong_signal_skat.yaml` | Burden 强信号测试 |
| `config/strong_signal_skato.yaml` | SKAT-O 强信号测试 |
| `config/regenie_example.yaml` | regenie 示例 |
| `config/quant_assoc_plink2_example.yaml` | PLINK2 示例 |
| `config/gwas_gene_catalog_demo.yaml` | GWAS Catalog 本地 TSV 示例 |
| `config/gwas_gene_catalog_api_demo.yaml` | GWAS Catalog API 示例 |
| `config/gwas_gene_catalog_official.yaml` | 官方完整 TSV + API 索引 |

---

## 通用参数

### project

```yaml
project:
  name: demo_genetic_support  # 项目名称
```

### analysis

```yaml
analysis:
  method: burden  # 方法名称：burden / skato / quant-assoc / gwas-gene-catalog
```

### input

```yaml
input:
  genotype_matrix: data/example/geno_matrix.tsv  # 基因型矩阵文件（TSV 格式）
  phenotype_file: data/example/pheno.tsv        # 表型文件
  covariate_file: data/example/covar.tsv       # 协变量文件
  # regenie 模式额外需要的文件：
  pred_file: data/regenie/example_pred.bin      # regenie step 1 输出的预测文件
```

### phenotype

```yaml
phenotype:
  name: phenotype      # 表型列名
  type: continuous    # 表型类型：continuous / binary
```

### output

```yaml
output:
  root_dir: output     # 输出根目录
```

---

## Burden 参数

### 基本参数

```yaml
burden:
  set_id: DEMO_GENE        # 分析对象名称（基因或区域标识）
  engine: skat             # 引擎：skat / regenie
```

### SKAT 引擎参数

```yaml
burden:
  engine: skat
  covariates:             # 协变量列表
    - age
    - sex
    - PC1
    - PC2
```

### regenie 引擎参数

```yaml
burden:
  engine: regenie
  regenie_bin: regenie                    # regenie 可执行文件路径，默认 "regenie"
  regenie_aaf_bins: 0.01                  # AAF 分箱，默认 0.01
  regenie_build_mask: max                # mask 聚合方式：max / sum
  regenie_ignore_pred: false              # 是否跳过 step 1 prediction 文件
  regenie_bsize: 200                      # step 2 block size，默认 200
  covariates:                            # 协变量列表
    - age
    - sex
```

### regenie 输入格式

```yaml
regenie:
  genotype_format: pfile  # 基因型格式：pfile (PGEN) / bfile (BED)
```

---

## SKAT-O 参数

### 基本参数

```yaml
skato:
  set_id: DEMO_GENE        # 分析对象名称
  engine: skat             # 引擎：skat / regenie
```

### SKAT 引擎参数

```yaml
skato:
  engine: skat
  covariates:             # 协变量列表
    - age
    - sex
    - PC1
    - PC2
```

### regenie 引擎参数

```yaml
skato:
  engine: regenie
  regenie_bin: regenie                    # regenie 可执行文件路径
  regenie_aaf_bins: 0.01                  # AAF 分箱
  regenie_build_mask: max                # mask 聚合方式
  regenie_ignore_pred: false              # 是否跳过 prediction 文件
  regenie_bsize: 200                      # block size
  covariates:                            # 协变量列表
    - age
    - sex
```

---

## 连续表型关联分析参数

```yaml
quant_assoc:
  plink2_bin: plink2                # plink2 可执行文件路径，默认 "plink2"
  covar_variance_standardize: true  # 是否对协变量做方差标准化，建议设为 true
```

---

## GWAS Catalog 基因注释参数

### 基本参数

```yaml
gwas_gene_catalog:
  project_result_file: data/gwas_gene_catalog/project_genes.tsv  # 本地基因级结果表
  gene_column: gene                  # gene symbol 所在列名
  source_mode: tsv                   # 数据源模式：tsv / api
```

### TSV 模式参数

```yaml
gwas_gene_catalog:
  source_mode: tsv
  gwas_reference_file: data/gwas_gene_catalog/gwas_catalog_reference.tsv  # 本地参考 TSV
  # 或使用远程下载：
  gwas_catalog_tsv_url: https://ftp.ebi.ac.uk/pub/databases/gwas/...  # 远程 TSV 地址
  match_fields:                      # 匹配字段范围
    - mapped_genes
    - reported_genes
```

### API 模式参数

```yaml
gwas_gene_catalog:
  source_mode: api
  gwas_catalog_api_base_url: https://www.ebi.ac.uk/gwas/rest/api/v2/associations  # API 地址
  api_page_size: 100                # 单页大小，默认 100
  api_max_pages: 5                  # 最大页数，不设则拉取全部
  api_extended_geneset: false       # 是否启用 extended_geneset 参数
  # 如需同时支持 reported_genes，补充索引：
  gwas_reference_file: data/gwas_gene_catalog/gwas_catalog_reference.tsv
  match_fields:
    - mapped_genes
    - reported_genes
```

---

## 参数速查表

| 参数 | 所属方法 | 说明 |
|------|----------|------|
| `analysis.method` | 通用 | 方法名称 |
| `phenotype.name` | 通用 | 表型列名 |
| `phenotype.type` | 通用 | 表型类型 |
| `input.genotype_matrix` | 通用 | 基因型矩阵文件 |
| `input.phenotype_file` | 通用 | 表型文件 |
| `input.covariate_file` | 通用 | 协变量文件 |
| `output.root_dir` | 通用 | 输出目录 |
| `burden.set_id` | burden | 分析对象名称 |
| `burden.engine` | burden | 引擎类型 |
| `burden.covariates` | burden | 协变量列表 |
| `burden.regenie_*` | burden | regenie 相关参数 |
| `skato.set_id` | skato | 分析对象名称 |
| `skato.engine` | skato | 引擎类型 |
| `skato.covariates` | skato | 协变量列表 |
| `skato.regenie_*` | skato | regenie 相关参数 |
| `quant_assoc.covar_variance_standardize` | quant-assoc | 协变量方差标准化 |
| `gwas_gene_catalog.project_result_file` | gwas-gene-catalog | 项目结果表 |
| `gwas_gene_catalog.gene_column` | gwas-gene-catalog | 基因列名 |
| `gwas_gene_catalog.source_mode` | gwas-gene-catalog | 数据源模式 |
| `gwas_gene_catalog.api_*` | gwas-gene-catalog | API 相关参数 |

---

## 参考链接

- [GWAS Catalog API 文档](https://www.ebi.ac.uk/gwas/rest/api/v2/docs)
- [GWAS Catalog 使用教程](https://www.ebi.ac.uk/gwas/docs/tut/api)
- [PLINK2 官方文档](https://www.cog-genomics.org/plink/2.0/)