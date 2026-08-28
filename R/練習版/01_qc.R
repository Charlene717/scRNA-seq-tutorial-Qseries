# =====================================================================
# 01_qc.R — 練習腳本 1：讀入 GBM 5k、三條 QC 線、doublet 標記
#
# 對應影片：Q2 頁 8–22（§0 SoupX 選做、§1 讀檔與先摸一遍、§2 三指標、§3 MAD 三條線、§4 scDblFinder）
# 輸入：data/gbm5k/filtered_feature_bc_matrix/（由 00_setup.R 準備）
# 輸出：output/gbm_raw.rds（過濾前）、output/gbm_qc.rds（過濾後 + doublet 標記）
# 時間：約 3–5 分鐘
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)

## ---- 0. ambient-rna (optional) ------------------------------------- Q2 頁 10–12
# 需要 Cell Ranger 的 raw + filtered 兩個矩陣同時在同一個 outs/ 資料夾。
# 10x 官網的 GBM 5k 只提供 filtered，所以這一節預設不執行（RUN_SOUPX <- FALSE）；
# 用自己的資料時改成 TRUE，路徑指到 cellranger 的 outs/。
RUN_SOUPX <- FALSE
if (RUN_SOUPX) {
  library(SoupX)
  sc <- load10X("data/my_sample/outs/")                 # 自動讀 raw + filtered
  tmp <- CreateSeuratObject(sc$toc) |> NormalizeData() |> FindVariableFeatures() |>
         ScaleData() |> RunPCA(verbose = FALSE) |> FindNeighbors(dims = 1:20) |>
         ## TODO ▶ 分群的 resolution 從掃描與 clustree 挑哪一層？（Q2 頁 33–34）
         FindClusters(resolution = ____)                 # 粗分群給 SoupX 用
  sc  <- setClusters(sc, setNames(tmp$seurat_clusters, colnames(tmp)))
  sc  <- autoEstCont(sc)                                # 估污染比例 rho
  message("estimated rho = ", round(sc$fit$rhoEst, 3))  # < 0.05 可略過；> 0.10 建議校正
  adj <- adjustCounts(sc, roundToInt = TRUE)            # 校正後整數 counts
  saveRDS(adj, "output/counts_soupx.rds")
}

## ---- 1. load ------------------------------------------------------- Q2 頁 8–9
mtx.dir  <- "data/gbm5k/filtered_feature_bc_matrix"
stopifnot("找不到資料，請先跑 00_setup.R" = dir.exists(mtx.dir))
gbm.data <- Read10X(data.dir = mtx.dir)               # 吃資料夾，不是單一檔案
dim(gbm.data)                                         # 基因 × 細胞（約 36,601 × 5,604）

gbm <- CreateSeuratObject(counts = gbm.data, project = "gbm5k",
                          min.cells = 3,              # 基因：< 3 顆細胞表現 → 刪
                          min.features = 200)         # 細胞：< 200 個基因 → 刪（粗篩）
gbm

# 先摸一遍：維度、稀疏度、粒線體基因命名（人 MT-、鼠 mt-）——Q2 頁 9
dim(gbm)
cnt <- LayerData(gbm, layer = "counts")
1 - Matrix::nnzero(cnt) / prod(dim(cnt))                # 稀疏度（稀疏矩陣不能直接 mean(cnt == 0)），> 0.9 正常
summary(Matrix::colSums(cnt))                           # 每顆細胞 UMI 總數
grep("^MT-", rownames(gbm), value = TRUE)               # 空的 → 物種或基因名有問題
cnt[c("PTPRC", "SOX2", "MBP"), 1:5]

## ---- 2. metrics ---------------------------------------------------- Q2 頁 14–16
## TODO ▶ 粒線體基因的前綴是什麼？人與鼠不同（Q2 頁 9）
gbm[["percent.mt"]] <- PercentageFeatureSet(gbm, pattern = "____")
gbm[["percent.rb"]] <- PercentageFeatureSet(gbm, pattern = "^RP[SL]") # 核糖體，供參考
saveRDS(gbm, "output/gbm_raw.rds")                    # 過濾前先存：截斷之後回不去

p <- VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, pt.size = 0.05)
ggsave("output/figs/01_qc_violin.png", p, width = 11, height = 4.5, dpi = 150, bg = "white")

p <- FeatureScatter(gbm, "nCount_RNA", "percent.mt") +
     FeatureScatter(gbm, "nCount_RNA", "nFeature_RNA")
ggsave("output/figs/01_qc_scatter.png", p, width = 11, height = 4.5, dpi = 150, bg = "white")
# 看圖：percent.mt 有沒有第二個峰？散點圖三個角落各是什麼？

## ---- 3. thresholds ------------------------------------------------- Q2 頁 17–19
# 閾值跟著「這份資料」的分布走：median ± 3×MAD。mad() 已含 1.4826 校正。
## TODO ▶ MAD 的倍數：3 是慣例，為什麼不是 2？（Q2 頁 17–18）
mad.hi <- function(x, k = ____) median(x) + k * mad(x)
## TODO ▶ MAD 的倍數：3 是慣例，為什麼不是 2？（Q2 頁 17–18）
mad.lo <- function(x, k = ____) median(x) - k * mad(x)

mt.hi <- mad.hi(gbm$percent.mt)                       # percent.mt 直接算
## TODO ▶ nFeature / nCount 要在哪個尺度上算 MAD？（Q2 頁 17）
nf.lo <- 10^mad.lo(____(gbm$nFeature_RNA))
## TODO ▶ nFeature / nCount 要在哪個尺度上算 MAD？（Q2 頁 17）
nc.hi <- 10^mad.hi(____(gbm$nCount_RNA))
thr <- c(mt.hi = mt.hi, nf.lo = nf.lo, nc.hi = nc.hi)
print(round(thr, 1))

# 把線畫在分布上，確認 mt.hi 落在「主體與第二個峰之間」
p <- ggplot(gbm@meta.data, aes(percent.mt)) + geom_histogram(bins = 80, fill = "#2563A8", alpha = .6) +
     geom_vline(xintercept = 5, linetype = 3) + geom_vline(xintercept = mt.hi, colour = "#C0392B") +
     annotate("text", x = mt.hi, y = Inf, vjust = 2, hjust = -0.1, colour = "#C0392B",
              label = sprintf("median + 3×MAD = %.1f%%", mt.hi)) + theme_classic()
ggsave("output/figs/01_mt_threshold.png", p, width = 7, height = 4, dpi = 150, bg = "white")

n.before <- ncol(gbm)
gbm <- subset(gbm, subset = percent.mt < mt.hi & nFeature_RNA > nf.lo & nCount_RNA < nc.hi)
n.after  <- ncol(gbm)
cat(sprintf("QC：%d → %d 顆（濾掉 %d，%.1f%%）\n", n.before, n.after, n.before - n.after,
            100 * (n.before - n.after) / n.before))
# 方法段模板（把數字換成你的）：
# 「細胞保留條件：percent.mt < median + 3×MAD（= X%）、nFeature 與 nCount 在 log10 尺度
#   median ± 3×MAD 內；N0 顆中保留 N1 顆。」

## ---- 4. doublets --------------------------------------------------- Q2 頁 20–22
library(scDblFinder)
sce <- scDblFinder(as.SingleCellExperiment(gbm))      # dbr 依細胞數自動估（10x ≈ 0.8%/千顆）
gbm$dbl.score <- sce$scDblFinder.score
gbm$dbl.class <- sce$scDblFinder.class
table(gbm$dbl.class)
# 先標記、不刪；02_cluster.R §4 分群後看「哪一群整群是 doublet」再決定。

## ---- 5. save -------------------------------------------------------
saveRDS(gbm, "output/gbm_qc.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 1
#  1-1 把 k 從 3 改成 2.5 與 4，各濾掉幾顆？濾掉的比例差多少？哪一個更合理？（畫圖說明）
#  1-2 用 5% 當 percent.mt 上限重跑 subset：濾掉幾顆？用 VlnPlot(group.by=…) 看被砍掉的
#      細胞 nFeature 分布——它們像垂死細胞，還是像活細胞？
#  1-3 percent.rb（核糖體）在這份資料的分布如何？要不要用它當 QC 指標？查一篇文獻支持你的決定。
#  1-4 scDblFinder 給的 doublet 率是多少？用 10x 的經驗法則（每千顆 0.8%）估 5,600 顆的預期率，
#      兩者接近嗎？若差很多，可能的原因？
#  1-5 （有 raw 矩陣者）跑 §0：rho 是多少？校正前後 MBP 在免疫細胞群的 pct 差多少？
#  進階 把 1-1 的三個 k 值各存成一個物件，跑完 02_cluster.R 後比較分群結果是否穩定。
# =====================================================================
