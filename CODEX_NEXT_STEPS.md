# Hand this entire repository to Codex in VS Code

## Copy-paste Codex instruction (中文)

请在当前打开的 **Aging Genomics Lab v0.1.0** 仓库上实施 P1 增量开发，不要重新创建项目、替换 UI、清空目录或删除现有样例。首先阅读 README.md、AGENTS.md、docs/STATISTICAL_SCOPE.md，运行现有 pytest 和 CLI smoke test，记录基线结果。若缺依赖，先安装到 `.venv`。确保完成后 Windows `scripts/start_windows.bat` 仍能启动网页。

### P1 — 真实人类衰老转录组学案例（优先完成）

目标：增加真正可复现的公开数据分析案例，优先 NIH/NIA GEO **GSE164471**（健康人骨骼肌跨年龄 RNA-seq），备选 **GSE5086**（人类骨骼肌微阵列，适合规范化 log2 表达）。

1. 检查 GEO 官方页面和原论文的具体文件、样本来源、年龄、性别、平台、实验设计、重复测量、原始计数/normalized 数据类型和下载许可。**先形成 provenance + dataset card；不得猜字段或混淆数据类型。**
2. 建立 `src/aginggenomics/ingest/geo.py`、`scripts/fetch_public_dataset.py`，使用稳定且合法的官方下载链接；实现缓存、断点容错（如可行）、SHA-256、文件大小限制、失败重试、明确源数据引用。下载失败时不得生成伪造的“真实研究数据”。
3. 若选择 RNA-seq raw counts，增加单独的 **R/DESeq2** 分析适配器（可用 Docker/renv 保证可重复），有 >=3 独立生物样本/组时使用合理的设计矩阵，报告 size-factor normalization、dispersion、log2FC、BH-FDR、covariates 及版本。不能把 counts 直接交给当前 Welch 检验。
4. 若选择微阵列，验证数据是否已经经过 RMA/适当归一化与 log2 转换，并记录探针到基因的映射；实现谨慎的年龄分析或年轻/年长比较。优先按文献设计复现，不能以年龄组划分制造不合理的独立“年龄效应”结论。
5. 前端新增 **Public Study** 选择入口、数据集卡片（GEO accession、组织、样本量、数据类型、引用、处理方式）、下载/导入、真实/模拟明确区分，使用既有 Dashboard/CLI/报告组件。模拟数据继续默认可离线演示。
6. 新增真实数据的 QC、样本/基因过滤日志、协变量查看、批次可视化和结果导出；统计方法随输入类型明确显示，不允许静默降级到不合适模型。
7. 加单元测试、网络失败/CSV/schema 测试、fixtures、可复现的集成冒烟测试，更新 README 和 data source/provenance docs。大文件放 `data/cache/` 并加入 `.gitignore`，不要把真实数据整个提交到 GitHub。

### P1 验收

- 原有 pytest/CLI/Streamlit 功能全部保留，并通过新增测试。
- 能基于 **真实公开 GEO 原始数据**下载、清洗、形成可解释的统计结果；断网时模拟 demo 仍可工作。
- raw counts 无论是通过网页还是 CLI，都不会使用 Welch t-test。
- 输出结果清楚地标明真实/模拟、样本 provenance、数据尺度、计算模型、协变量、参数和软件版本。
- 提供可照做的 Windows VS Code 运行步骤，不要求用户提供 API key。
- Codex 完成后报告具体修改、测试数量、结果路径、实际使用的公开研究 accession 和局限性。**如未实际下载、验证真实数据，请明确未完成，不得宣称已完成。**

## After P1 (do not implement yet)

**P2: EWAS / methylation.** Review the GSE40279 accession, paper, public-use conditions, age/sex metadata, platform annotation and cell-composition inputs. Pin official assets and build a separate methylation data contract for beta values and sample metadata. Define probe and sample QC, a prespecified covariate-aware association design, multiple-testing correction and sensitivity analyses. Validate model assumptions, compare effect estimates with the original study where feasible, and export provenance and exclusion logs. Acceptance requires offline fixture tests, a verified public-data run and separate reporting in the workbench. Epigenetic-age modeling can follow after the EWAS validation.

**P3: GWAS/TWAS.** Public summary statistics import, Manhattan/QQ, valid coordinate build, significance thresholds, optional locus summaries; no private genotype data, no fabricated ancestry claims.

**P4: Stanford gene regulation.** ENCODE regulatory region BED/peak integration, reference genome tracking, interpretable sequence/feature model with clear baseline comparisons and leakage checks.

**P5: HPC/workflow.** Nextflow/Snakemake, Slurm profile, Docker, data lineage, robust chunking and performance benchmarks; optional cloud read-only demo.
