# =====================================================================
# 00_setup.R — Q 系列（腫瘤快速上手）練習腳本 0：環境、資料下載、專案結構
#
# 對應影片：Q2 頁 7–8、Q3 頁 8
# 執行方式：在 RStudio 開啟專案（.Rproj），從專案根目錄逐段執行。
#           第一次執行約 10–20 分鐘（安裝套件 + 下載約 60 MB 資料）。
# 練習資料：
#   (1) 10x GBM 5k  — Human Glioblastoma Multiforme 3'v3, Cell Ranger 4.0.0, 5,604 cells
#   (2) GSE84465    — Darmanis et al. 2017, 4 GBM patients, core vs periphery, Smart-seq2
# =====================================================================

## ---- 1. 套件 -------------------------------------------------------
# CRAN
cran <- c("Seurat", "dplyr", "ggplot2", "patchwork", "clustree", "harmony", "msigdbr", "remotes",
          "mclust", "SoupX", "data.table", "survival", "survminer", "devtools", "ashr", "reshape2", "ggrepel")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c("scDblFinder", "SingleR", "celldex", "GEOquery", "DESeq2", "apeglm", "fgsea", "speckle", "infercnv",
          "EnhancedVolcano", "slingshot", "tradeSeq", "decoupleR", "TCGAbiolinks", "UCell", "AUCell",
          "TOAST",   # TOAST：MuSiC 的相依套件
          "clusterProfiler", "org.Hs.eg.db", "enrichplot",   # 06 的 ORA（GO / KEGG）與富集圖
          "MAST")                                            # 06b 的 cell-level DE（病人當共變量）
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
# GitHub 上的套件：lisi（04 §5）、CellChat（07）、liana（07 §4）、MuSiC（08 §3）
# 常見錯誤：「HTTP error 401 Bad credentials」= 環境變數 GITHUB_PAT 裡有一個過期／無效的 token。
# 下面的 helper 會先用匿名下載（不需要 token）；匿名有每小時 60 次的限制，超過再設定有效的 PAT。
install_gh <- function(repo, pkg = basename(repo)) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  old <- Sys.getenv("GITHUB_PAT"); Sys.unsetenv("GITHUB_PAT"); on.exit(if (nzchar(old)) Sys.setenv(GITHUB_PAT = old))
  ok <- tryCatch({ remotes::install_github(repo, upgrade = "never"); TRUE },
                 error = function(e) { message("  ! ", repo, " 安裝失敗：", conditionMessage(e)); FALSE })
  if (!ok) message("    → 稍後手動執行：remotes::install_github(\"", repo, "\")；仍失敗請檢查 usethis::gh_token_help()")
  invisible(ok)
}
install_gh("immunogenomics/lisi")
install_gh("jinworks/CellChat")          # 依賴多（NMF、circlize、ComplexHeatmap…），約 5–10 分鐘
install_gh("saezlab/liana")
install_gh("xuranw/MuSiC")
# —— 選配（08 軌跡的 Monocle 比較；不裝也能跑 08 的 §1，§2–3 會自動跳過）——
# install_gh("satijalab/seurat-wrappers", "SeuratWrappers")   # Seurat → cell_data_set 的橋
# install_gh("cole-trapnell-lab/monocle3")                    # Monocle3（相依較多，Windows 需 Rtools）
# BiocManager::install("monocle", update = FALSE, ask = FALSE) # Monocle2（經典 DDRTree）
# infercnv 底層的 rjags 需要「系統層級」的 JAGS 程式（不是 R 套件，R 裝不了它；裝完要重開 R）：
#   macOS  : brew install jags
#   Ubuntu : sudo apt-get install jags
#   Windows: https://sourceforge.net/projects/mcmc-jags/  安裝後重開 R

library(Seurat)
stopifnot(packageVersion("Seurat") >= "5.0.0")   # 本系列全程 Seurat v5
set.seed(1234)                                    # 全系列固定 seed
cat("Seurat", as.character(packageVersion("Seurat")), "\n")

## ---- 2. 專案結構 ---------------------------------------------------
# 一律相對路徑；不要 setwd()。
for (d in c("data", "R", "output", "output/figs")) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 3. 資料 (1)：10x GBM 5k --------------------------------------
# 官方頁面：10x Genomics Datasets → "Human Glioblastoma Multiforme: 3'v3 Whole Transcriptome Analysis"
# 若下列直連網址失效，請到該頁面「Output and supplemental files」下載同名檔案放到 data/。
gbm.url <- paste0("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/",
                  "Parent_SC3v3_Human_Glioblastoma/",
                  "Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz")
gbm.tar <- "data/gbm5k_filtered_feature_bc_matrix.tar.gz"
if (!file.exists(gbm.tar)) {
  options(timeout = 600)
  download.file(gbm.url, gbm.tar, mode = "wb")     # 約 30 MB
}
if (!dir.exists("data/gbm5k/filtered_feature_bc_matrix")) {
  untar(gbm.tar, exdir = "data/gbm5k")
}
list.files("data/gbm5k/filtered_feature_bc_matrix")
# 預期：barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

## ---- 4. 資料 (2)：GSE84465（Darmanis 2017）-------------------------
# counts（基因 × 細胞，read counts，約 20 MB）
geo.url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz"
geo.csv <- "data/GSE84465_GBM_All_data.csv.gz"
if (!file.exists(geo.csv)) {
  options(timeout = 600)
  download.file(geo.url, geo.csv, mode = "wb")
}
# metadata 由 GEOquery 於 04_multipatient.R 內下載（series matrix）。
# 若 GEO FTP 直連失敗：到 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465
# 的 Supplementary file 區手動下載同名檔案放到 data/。

## ---- 5. inferCNV 基因座標檔 -----------------------------------------
# 官方提供的 hg38 gencode v27 基因座標（gene, chr, start, end）
gpos.url <- "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt"
gpos     <- "data/hg38_gencode_v27.txt"
if (!file.exists(gpos)) try(download.file(gpos.url, gpos, mode = "wb"))
# 若下載失敗，見 05_infercnv.R §1 的替代作法（用 EnsDb 自行產生）。

## ---- 6. 檢查 --------------------------------------------------------
cat("\n== 檔案檢查 ==\n")
for (f in c(gbm.tar, "data/gbm5k/filtered_feature_bc_matrix/matrix.mtx.gz", geo.csv, gpos))
  cat(sprintf("%-60s %s\n", f, ifelse(file.exists(f), "OK", "缺少")))
sessionInfo()

# =====================================================================
# ▶ 練習 0（暖身，5 分鐘）
#  0-1 用 readLines() 讀 features.tsv.gz 的前 5 行，三個欄位各是什麼？
#  0-2 用 R.utils::countLines() 或 readLines() 數 barcodes.tsv.gz 有幾行；
#      它應該等於 Cell Ranger 判定的細胞數。跟 10x 頁面寫的 5,604 一樣嗎？
#  0-3 GSE84465 的 csv 最後 5 列是什麼？（提示：它們不是基因）
# =====================================================================
