# Analysis code — BNB_small multivariable logistic regression

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22097248.svg)](https://doi.org/10.5281/zenodo.22097248)

R code for a multivariable logistic regression analysis of good functional
outcome in convalescent rehabilitation inpatients, using restricted cubic
splines (RCS) for non-linearity and substantive-model-compatible multiple
imputation (`smcfcs`) for missing covariates.

> **No patient data are included in this repository.**
> The individual patient records underlying the analysis are held under
> institutional data-governance rules and are not redistributed. Every script
> here reads from a local `data/` directory that is deliberately excluded by
> `.gitignore`. Running the code therefore reproduces the *procedure*, not the
> published numbers.

---

## Repository layout

```
.
├── 00_labels.R            English display labels, sourced by every later script
├── 01_preprocess.R … 21_welch_discharge_fim_total.R
├── LICENSE                MIT — applies to the code
├── LICENSE-docs.txt       CC BY 4.0 — applies to documentation, figures, tables
├── CITATION.cff
└── (created at run time, not tracked)
    ├── data/              inputs and intermediate .rds / .csv objects
    └── figures/           generated .png / .pdf / .tiff
```

Scripts locate the project root with the **`here`** package, so `00_labels.R`
must stay at the repository root. Open the folder as an RStudio project (or
create an empty `.here` file) before running anything.

## Requirements

R ≥ 4.2 and the following packages:

`Hmisc`, `MASS`, `ResourceSelection`, `brglm2`, `car`, `data.table`,
`detectseparation`, `flextable`, `forcats`, `furrr`, `future`, `ggcorrplot`,
`ggplot2`, `gridExtra`, `here`, `logistf`, `magick`, `mice`, `officer`,
`pROC`, `patchwork`, `renv`, `smcfcs`, `tictoc`

```r
install.packages(c(
  "Hmisc", "MASS", "ResourceSelection", "brglm2", "car", "data.table",
  "detectseparation", "flextable", "forcats", "furrr", "future", "ggcorrplot",
  "ggplot2", "gridExtra", "here", "logistf", "magick", "mice", "officer",
  "pROC", "patchwork", "renv", "smcfcs", "tictoc"
))
```

`17_session_info_license.R` writes the exact package versions used
(`sessionInfo.txt`, `package_versions.csv`) for the record.

## Running order

Run the scripts in numeric order. Sub-numbered scripts (`03-1`, `04-1`, …) are
optional re-plotting or curation steps that depend on the parent step and never
overwrite it.

| Script | Step |
|---|---|
| `00_labels.R` | Shared English labels for variables, levels and pipelines |
| `01_preprocess.R` | Data preprocessing |
| `02_mnar_systems.R` | Two sensitivity-analysis systems for MNAR |
| `03_centering_corr_3.R` | Median centring and collinearity screening |
| `03-1_corr_combined.R` | Combined correlation heat map across subsets |
| `04_linearity_rcs.R` | Linearity assessment and conversion to RCS |
| `04-1_linearity_replot_1.R` | Re-plot of the linearity figures for publication |
| `05_imputation_1.R` | Multiple imputation (MAR: JCS, MMSE, BBS) via `smcfcs` |
| `05-1_imputation_diagnostics.R` | Convergence traces, observed vs imputed distributions |
| `05-1-1_final_model_fmi_mcse.R` | Fraction of missing information and Monte-Carlo SE |
| `05-2_imputation_diagnostic_plots.R` | Imputation diagnostic figures, subsets combined |
| `06_separation_estimation.R` | Separation detection and choice of estimator |
| `06-1_separation_report.R` | Print the separation decision |
| `07_influence_diagnostics_3.R` | Influence diagnostics for the full model |
| `07-1_outlier_definition_1.R` | Fix the outlier definition from those diagnostics |
| `08_variable_selection_1.R` | Provisional model — backward selection |
| `08-1_user_model_2.R` | Manual curation of the provisional model |
| `09_interaction_screening_1.R` | Screening of two-way interaction terms |
| `10_final_model_5.R` | Final model — forward selection of interactions |
| `10-1_user_model_1.R` | Manual curation of the final model |
| `11_or_forest_2.R` | Pooled coefficient table and forest plots |
| `11-1_orcurve_provisional_1.R` | Odds-ratio curves for continuous predictors |
| `11-1-1_combine_orcurve_canvases_1.R` | Combine the OR curves onto single canvases |
| `11-2_interaction_plots.R` | Plots of the retained interactions |
| `11-2-1_merge_11-2.R` | Combine the interaction plots onto single canvases |
| `12_model_performance.R` | Calibration and discrimination |
| `12-1_gvif_collinearity.R` | GVIF decomposition for high-collinearity models |
| `12-2_variable_stability.R` | glm vs Firth estimation path, variable stability |
| `12-2-1_selection_bootstrap_optimism.R` | Optimism correction with selection inside each bootstrap |
| `12-3_implementation_performance.R` | Performance summary table for the implementation model |
| `12-4_implementation_plots.R` | Calibration and ROC panels |
| `13_table1_4.R` | Table 1 — baseline characteristics |
| `14_excluded_vs_included.R` | Analysed cohort vs patients who died or were transferred |
| `15_reporting_supplements-1.R` | Reporting quantities: minimum sample size, crude ORs, absolute risks |
| `15-1_table_s3_crude_multilevel.R` | Crude odds ratios recomputed so that a k-level categorical predictor contributes all k−1 contrasts |
| `16_table2_model_specification-1.R` | Table 2 — full model specification for the four pipelines |
| `17_session_info_license.R` | `sessionInfo`, package versions, licence files |
| `18_show_models.R` | Print the currently stored main-effect and final models |
| `19_external_validation_iecv.R` | Internal–external cross-validation |
| `20_fig2_fig3_descriptive.R` | Descriptive figures: outcome by age and disease class, and the subset forest plot |
| `21_welch_discharge_fim_total.R` | Welch two-sample *t* test of the total discharge FIM, oldest-old versus younger patients |

`15_reporting_supplements.R` and `16_table2_model_specification.R` are the
earlier revisions of the two `-1` scripts above and are kept only for
provenance; use the `-1` versions.

`15-1_table_s3_crude_multilevel.R` replaces the crude odds-ratio table produced
by section 15.2 of `15_reporting_supplements-1.R`, which reported only the first
contrast of each categorical predictor. It is a stand-alone script: it neither
modifies nor has to be run together with the parent script, and it leaves the
other outputs of section 15 untouched.

Output file names follow each script's own numbering, which is not the figure
numbering of the manuscript: `20_fig2_fig3_descriptive.R` writes
`fig2_outcome_by_age_disease` and `fig3_subset_forest`, which appear in the
manuscript as Figures 1 and 2.

## Note on comments

The source comments were removed before publication, so these files carry the
executed logic only. This README and the step table above are the intended
documentation.

## Licence

- **Code** (`*.R`): MIT — see [`LICENSE`](LICENSE)
- **Documentation, figures, tables, reporting checklists**: CC BY 4.0 — see [`LICENSE-docs.txt`](LICENSE-docs.txt)
- **Individual patient data**: not covered by either licence and not redistributed

## Citation

See [`CITATION.cff`](CITATION.cff). The repository is archived on Zenodo under
the concept DOI [10.5281/zenodo.22097248](https://doi.org/10.5281/zenodo.22097248),
which always resolves to the newest release. When citing the code, use the
**version** DOI shown on the Zenodo record for the release you actually used.

---

## 日本語

回復期リハビリテーション入院患者の良好な機能転帰を対象とした多変量ロジスティック
回帰の解析コードです。非線形性は制限付き 3 次スプライン（RCS）、欠測共変量は
substantive-model-compatible な多重代入（`smcfcs`）で扱っています。

**個票データは含まれていません。** 各スクリプトはローカルの `data/` を読みますが、
このディレクトリは `.gitignore` で除外されています。コードは手続きを再現するもので
あり、公表値そのものを再現するものではありません。

スクリプトは番号順に実行してください。枝番（`03-1`、`04-1` など）は親ステップに
依存する再描画・手動キュレーション用で、親ステップを上書きしません。`here`
パッケージでプロジェクトルートを解決するため、`00_labels.R` はリポジトリ直下に
置いたままにしてください。

`15-1_table_s3_crude_multilevel.R` は、`15_reporting_supplements-1.R` の 15.2 節が
各カテゴリ変数の第 1 対比しか出力しなかった点を修正し、k 水準の変数について k−1
対比すべてを報告する粗オッズ比表を出力します。親スクリプトを書き換えるものでは
なく、単独で実行できます。

コメントは公開前に全削除しています。上記のステップ表がドキュメントの役割を
果たします。
