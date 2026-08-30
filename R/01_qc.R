# =====================================================================
# 01_qc.R — 練習腳本 1：讀入 GBM 5k、三條 QC 線、doublet 標記
#
# 對應影片：Q2 頁 8–22（§0 SoupX 選做、§1 讀檔與初步檢視、§2 三指標、§3 MAD 三條線、§4 DoubletFinder）
# 輸入：data/gbm5k/filtered_feature_bc_matrix/（由 00_setup.R 準備）
# 輸出：output/rds/01_gbm_raw.rds（過濾前）、output/rds/01_gbm_qc.rds（過濾後 + doublet 標記）
# 時間：約 8–15 分鐘（DoubletFinder 的 pK 掃描本身就要數分鐘）
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)   # 圖檔資料夾先建好，後面 ggsave 才不會失敗

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
         FindClusters(resolution = 0.5)                 # 粗分群給 SoupX 用
  sc  <- setClusters(sc, setNames(tmp$seurat_clusters, colnames(tmp)))
  sc  <- autoEstCont(sc)                                # 估污染比例 rho
  message("estimated rho = ", round(sc$fit$rhoEst, 3))  # < 0.05 可略過；> 0.10 建議校正
  adj <- adjustCounts(sc, roundToInt = TRUE)            # 校正後整數 counts
  saveRDS(adj, "output/rds/01_counts_soupx.rds")
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

# 初步檢視：維度、稀疏度、粒線體基因命名（人 MT-、鼠 mt-）——Q2 頁 9
dim(gbm)
cnt <- LayerData(gbm, layer = "counts")
1 - Matrix::nnzero(cnt) / prod(dim(cnt))                # 稀疏度（稀疏矩陣不能直接 mean(cnt == 0)）
# 判讀：多數 10x 資料 > 0.9；腫瘤（惡性細胞 RNA 量大、深度高）常見 0.85–0.90，本份約 0.86 屬正常。
# 真正該警覺的是 < 0.8（可能讀到 raw 矩陣或多重樣本合併）或 > 0.98（深度太淺）。
summary(Matrix::colSums(cnt))                           # 每顆細胞 UMI 總數
grep("^MT-", rownames(gbm), value = TRUE)               # 空的 → 物種或基因名有問題
cnt[intersect(c("PTPRC", "SOX2", "MBP"), rownames(cnt)), 1:5]   # 先 intersect：min.cells 可能已刪掉其中之一

## ---- 2. metrics ---------------------------------------------------- Q2 頁 14–16
gbm[["percent.mt"]] <- PercentageFeatureSet(gbm, pattern = "^MT-")   # 人類 MT-；小鼠 "^mt-"
gbm[["percent.rb"]] <- PercentageFeatureSet(gbm, pattern = "^RP[SL]") # 核糖體，供參考
saveRDS(gbm, "output/rds/01_gbm_raw.rds")                    # 過濾前先存：截斷之後回不去

p <- VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, pt.size = 0.05)
ggsave("output/figs/01_qc_violin.png", p, width = 11, height = 4.5, dpi = 150, bg = "white")

p <- FeatureScatter(gbm, "nCount_RNA", "percent.mt") +
     FeatureScatter(gbm, "nCount_RNA", "nFeature_RNA")
ggsave("output/figs/01_qc_scatter.png", p, width = 11, height = 4.5, dpi = 150, bg = "white")
# 看圖：percent.mt 有沒有第二個峰？散點圖三個角落各是什麼？

## ---- 3. thresholds ------------------------------------------------- Q2 頁 17–19
# 閾值跟著「這份資料」的分布走：median ± 3×MAD。mad() 已含 1.4826 校正。
mad.hi <- function(x, k = 3) median(x) + k * mad(x)
mad.lo <- function(x, k = 3) median(x) - k * mad(x)

mt.hi <- mad.hi(gbm$percent.mt)                       # percent.mt 直接算
nf.lo <- 10^mad.lo(log10(gbm$nFeature_RNA))           # nFeature / nCount 是對數常態 → 在 log10 尺度算
nc.hi <- 10^mad.hi(log10(gbm$nCount_RNA))
thr <- c(mt.hi = mt.hi, nf.lo = nf.lo, nc.hi = nc.hi)
print(round(thr, 1))

# 把線畫在分布上，確認 mt.hi 落在「主體與第二個峰之間」
p <- ggplot(gbm@meta.data, aes(percent.mt)) + geom_histogram(bins = 80, fill = "#2563A8", alpha = .6) +
     geom_vline(xintercept = 5, linetype = 3) + geom_vline(xintercept = mt.hi, colour = "#C0392B") +
     annotate("text", x = mt.hi, y = Inf, vjust = 2, hjust = -0.1, colour = "#C0392B",
              label = sprintf("median + 3×MAD = %.1f%%", mt.hi)) + theme_classic()
ggsave("output/figs/01_mt_threshold.png", p, width = 7, height = 4, dpi = 150, bg = "white")

# 三條線各自砍掉幾顆（相加會大於總移除，因為有細胞同時違反兩條）
cat(sprintf("  percent.mt   >= %8.1f%%：%4d 顆\n", mt.hi, sum(gbm$percent.mt   >= mt.hi)),
    sprintf("  nFeature_RNA <= %8.0f ：%4d 顆\n", nf.lo, sum(gbm$nFeature_RNA <= nf.lo)),
    sprintf("  nCount_RNA   >= %8.0f ：%4d 顆\n", nc.hi, sum(gbm$nCount_RNA   >= nc.hi)), sep = "")
## >>> 參考答案 ------------------------------------------------------
# 判讀：本課這份 GBM 5k 的結果是 mt 222 顆、nFeature 74 顆、nCount 0 顆（重疊 51 顆）。
#       nCount 上限砍 0 顆——log10 尺度的 3×MAD 對深度分布寬的腫瘤資料很鬆，
#       算出來的 nc.hi（144,178）落在實際最大值（117,003）之外。
#       這代表「深度高端沒有離群細胞」，不是設定錯誤。
#       想要一條必定作用的上限：nc.hi <- quantile(gbm$nCount_RNA, 0.995)；
#       不設也合理——超高深度多半是 doublet，交給 §4 的 DoubletFinder 判斷更準。
#       方法段報告時，只寫實際有作用的條件。
## <<< 參考答案

n.before <- ncol(gbm)
gbm <- subset(gbm, subset = percent.mt < mt.hi & nFeature_RNA > nf.lo & nCount_RNA < nc.hi)
n.after  <- ncol(gbm)
cat(sprintf("QC：%d → %d 顆（濾掉 %d，%.1f%%）\n", n.before, n.after, n.before - n.after,
            100 * (n.before - n.after) / n.before))
# 方法段模板（把數字換成你的）：
# 「細胞保留條件：percent.mt < median + 3×MAD（= X%）、nFeature 與 nCount 在 log10 尺度
#   median ± 3×MAD 內；N0 顆中保留 N1 顆。」

## ---- 4. doublets --------------------------------------------------- Q2 頁 20–22
# 主工具：DoubletFinder（McGinnis et al. 2019, Cell Systems）。
# 原理：合成人工 doublet 混進資料，再問每顆細胞「鄰居裡有多少是人工 doublet」（pANN 分數）。
# 因為要在 PCA 空間找鄰居，它需要一個「已經跑過標準流程」的物件——
# 這裡先做一份暫時的 tmp 只給 DoubletFinder 用；正式的 HVG / PCA / 分群在 02 才做。
library(DoubletFinder)
tmp <- NormalizeData(gbm, verbose = FALSE) |>
       FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
       ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
tmp <- FindNeighbors(tmp, dims = 1:20, verbose = FALSE) |>
       FindClusters(resolution = 0.5, verbose = FALSE)     # 粗分群，只用來估 homotypic 比例

# (a) pK：DoubletFinder 唯一需要調的參數（鄰域大小）。掃一遍，取 BCmetric 最高的那個。
sweep.stats <- summarizeSweep(paramSweep(tmp, PCs = 1:20, sct = FALSE), GT = FALSE)
bcmvn        <- find.pK(sweep.stats)
bcmvn$pK.num <- as.numeric(as.character(bcmvn$pK))
pK.sel       <- bcmvn$pK.num[which.max(bcmvn$BCmetric)]
p <- ggplot(bcmvn, aes(pK.num, BCmetric)) + geom_line() + geom_point(size = 1.5) +
     geom_vline(xintercept = pK.sel, colour = "#C0392B", linetype = 2) + theme_classic() +
     labs(x = "pK", y = "BCmetric", title = sprintf("DoubletFinder pK sweep (selected pK = %s)", pK.sel))
ggsave("output/figs/01_doubletfinder_pK.png", p, width = 7, height = 4, dpi = 150, bg = "white")
cat("選到的 pK =", pK.sel, "\n")
# 看圖：BCmetric 有一個明確的單峰才算掃得乾淨；平坦或多峰代表資料結構弱，pK 的選擇會不穩。

# (b) 期望 doublet 數。10x 經驗法則：每裝載 1,000 顆細胞約 0.8%。
dbr      <- 0.008 * ncol(tmp) / 1000                  # 例：5,261 顆 → 約 4.2%
nExp     <- round(dbr * ncol(tmp))
homo     <- modelHomotypic(tmp$seurat_clusters)       # 同型別相撞看不出來，要從期望值扣掉
nExp.adj <- round(nExp * (1 - homo))
cat(sprintf("期望 doublet 率 %.1f%%；nExp = %d → homotypic 調整後 = %d\n", 100 * dbr, nExp, nExp.adj))
# 這是 DoubletFinder 與 scDblFinder 最大的差別：期望率要「你自己給」，不是工具替你估。
# 給錯就會系統性多殺或少殺——所以這個數字必須寫進方法段。

# (c) 執行（pN = 0.25 是預設；原論文顯示結果對 pN 不敏感，對 pK 才敏感）
tmp  <- doubletFinder(tmp, PCs = 1:20, pN = 0.25, pK = pK.sel, nExp = nExp.adj, sct = FALSE)
pann <- grep("^pANN_",               colnames(tmp@meta.data), value = TRUE)[1]
dfcl <- grep("^DF.classifications_", colnames(tmp@meta.data), value = TRUE)[1]
gbm$dbl.score <- tmp@meta.data[[pann]]                # pANN：鄰居中人工 doublet 的比例
gbm$dbl.class <- tolower(tmp@meta.data[[dfcl]])       # "doublet" / "singlet"
table(gbm$dbl.class)
rm(tmp); invisible(gc())

# （選配）第二種方法交叉比對：scDblFinder（Germain et al. 2021, F1000Research）。
# 原理相近但實作獨立；兩者準確度相當（Xi & Li 2021 的 benchmark 中 DoubletFinder 準確度最高，
# 但該篇未納入 scDblFinder）。實務上兩法取交集當「高信心 doublet」最穩。
RUN_SCDBL <- TRUE
if (RUN_SCDBL && requireNamespace("scDblFinder", quietly = TRUE)) {
  set.seed(1234)            # scDblFinder 有隨機性；不固定 seed 的話結果會隨前面消耗掉的亂數而飄
  sce <- scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(gbm))
  gbm$dbl.class.scdbl <- as.character(sce$scDblFinder.class)
  print(table(DoubletFinder = gbm$dbl.class, scDblFinder = gbm$dbl.class.scdbl))
  cat(sprintf("兩法一致率 %.1f%%；兩法都判 doublet %d 顆、只有單邊判到 %d 顆\n",
              100 * mean(gbm$dbl.class == gbm$dbl.class.scdbl),
              sum(gbm$dbl.class == "doublet" & gbm$dbl.class.scdbl == "doublet"),
              sum(gbm$dbl.class != gbm$dbl.class.scdbl)))
}

## >>> 參考答案 ------------------------------------------------------
# 本課這份 GBM 5k 的參考結果（你的數字應該一樣）：
#   pK = 0.11；nExp = 221 → homotypic 調整後 198；DoubletFinder 標 198 顆（3.8%）。
#   01_doubletfinder_pK.png 是「掃得乾淨」的範例：pK = 0.11 的 BCmetric 約 6,500，
#   次高的局部峰（pK ≈ 0.16）只有約 1,400——差 4 倍以上，峰位才可信。
#   若曲線平坦、或有兩三個高度相近的峰，代表 pK 選擇不穩，「哪 198 顆」的可信度就下降。
#   scDblFinder 自行找閾值後標 415 顆（7.9%），兩法交集 159 顆、一致率 94.4%。
#
# 判讀交叉表時注意：兩法的差異主要來自「閾值畫在哪」，不是「誰排在前面」——
# DoubletFinder 標記數由你給的 nExp 決定，scDblFinder 自己找閾值，通常會抓比較多。
# 所以「只有 scDblFinder 判到」的那一批，多半是排名在交界處的細胞，不是兩法根本不同意。
## <<< 參考答案

# 先標記、不刪。分群之後（02 §4）再用群層級的證據決定要不要整群移除。
# 依據：目前的最佳實務指南建議「先把判定出的 doublet 留在資料裡，視覺化與分群後再重新評估」
#       （sc-best-practices，Heumos et al. 2023）。

## ---- 5. save -------------------------------------------------------
saveRDS(gbm, "output/rds/01_gbm_qc.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 1
#  1-1 把 k 從 3 改成 2.5 與 4，各濾掉幾顆？濾掉的比例差多少？哪一個更合理？（畫圖說明）
#  1-2 用 5% 當 percent.mt 上限重跑 subset：濾掉幾顆？用 VlnPlot(group.by=…) 看被砍掉的
#      細胞 nFeature 分布——它們像垂死細胞，還是像活細胞？
#  1-3 percent.rb（核糖體）在這份資料的分布如何？要不要用它當 QC 指標？查一篇文獻支持你的決定。
#  1-4 §4 印出的 doublet 數，正好等於你給的 nExp.adj——這不是巧合。
#      (a) 把 dbr 的 0.008 改成 0.004 與 0.016 各重跑一次 §4，標記數各變成多少？
#      (b) 由此回答：DoubletFinder 報出來的「doublet 率」是資料的性質，還是你的假設？
#      (c) 看交叉表：scDblFinder 自己找閾值後判了幾顆？「只有 scDblFinder 判到」的那批細胞，
#          它們的 dbl.score（pANN）與 nFeature 落在什麼位置？這告訴你兩法的分歧來自哪裡？
#  1-5 （有 raw 矩陣者）跑 §0：rho 是多少？校正前後 MBP 在免疫細胞群的 pct 差多少？
#  進階 把 1-1 的三個 k 值各存成一個物件，跑完 02_cluster.R 後比較分群結果是否穩定。
# =====================================================================
