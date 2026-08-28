# =====================================================================
# 05_infercnv.R — 練習腳本 5：inferCNV 惡性判定、CNV 分數與相關、三角驗證、對答案
#
# 對應影片：Q3 頁 20–26（§1 輸入與執行、§2 兩個數字、§3 三角驗證）
# 輸入：output/gbm4_unintegrated.rds（04_multipatient.R；用「未整合」那份）
# 輸出：output/infercnv/（inferCNV 輸出）、output/gbm4_malignant.rds
# 時間：inferCNV 約 10–30 分鐘（denoise、無 HMM；視機器而定）
# 注意：inferCNV 底層的 rjags 需要「系統層級」的 JAGS 程式（不是 R 套件，R 裝不了它），
#       必須先在作業系統安裝 JAGS 4.x 再重開 R；未安裝的話本腳本會在 §1 直接停下並提示。
#       還沒裝 JAGS 前，06–08 可先用 celltype_author 的 Neoplastic 當替代惡性標籤測試（見 06 §0）。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/gbm4_unintegrated.rds")

## ---- 1. run-infercnv ----------------------------------------------- Q3 頁 20–22
# 前置檢查：inferCNV 依賴 rjags，而 rjags 需要「系統層級」安裝 JAGS 4.x（不是 R 套件）
#   Windows：到 https://sourceforge.net/projects/mcmc-jags/files/ 下載 JAGS-4.x.y.exe 安裝
#   macOS  ：brew install jags      Linux：sudo apt install jags
#   安裝完「重新啟動 R / RStudio」，再跑本腳本
if (inherits(try(suppressPackageStartupMessages(library(rjags)), silent = TRUE), "try-error")) {
  stop("找不到 JAGS：請先在系統安裝 JAGS 4.x（見上方註解），重開 R 後再執行 05_infercnv.R。",
       call. = FALSE)
}
library(infercnv)
# (a) 註釋檔：參考組 = 確定正常的細胞（免疫、寡樹突）；觀察組按病人分開
refs <- c("Immune cell", "Oligodendrocyte")
stopifnot(all(refs %in% gbm4$celltype_author))
ann <- data.frame(row.names = colnames(gbm4),
                  group = ifelse(gbm4$celltype_author %in% refs,
                                 gbm4$celltype_author,
                                 paste0("obs_", gbm4$patient)))
write.table(ann, "data/infercnv_annot.txt", sep = "\t", col.names = FALSE, quote = FALSE)
table(ann$group)

# (b) 基因座標檔：00_setup.R 已下載 hg38_gencode_v27.txt（gene  chr  start  end）
gpos <- "data/hg38_gencode_v27.txt"
if (!file.exists(gpos)) {
  # 替代作法：用 EnsDb 自行產生（BiocManager::install("EnsDb.Hsapiens.v86")）
  library(EnsDb.Hsapiens.v86)
  g <- genes(EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name", "gene_seq_start", "gene_seq_end"))
  g <- as.data.frame(g); g <- g[g$seq_name %in% c(1:22, "X", "Y"), ]
  g <- g[!duplicated(g$gene_name), c("gene_name", "seq_name", "start", "end")]
  g$seq_name <- paste0("chr", g$seq_name)
  write.table(g, gpos, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# (c) 原始 counts（不是 data、不是 integrated）
cts <- LayerData(gbm4, assay = "RNA", layer = "counts")
obj <- CreateInfercnvObject(raw_counts_matrix = cts,
                            annotations_file  = "data/infercnv_annot.txt",
                            gene_order_file   = gpos,
                            ref_group_names   = refs)
obj <- infercnv::run(obj,
                     cutoff = 1,                       # Smart-seq2 用 1；10x 用 0.1
                     out_dir = "output/infercnv",
                     cluster_by_groups = TRUE,         # 每位病人各自聚類
                     denoise = TRUE,
                     HMM = FALSE,                      # 亞株分析時再開（慢很多）
                     num_threads = 4)
# 輸出的 infercnv.png 就是熱圖：上半參考（平）、下半四位病人；看 chr7 gain / chr10 loss。

## ---- 2. cnv-score-cor ---------------------------------------------- Q3 頁 23–24
obs <- read.table("output/infercnv/infercnv.observations.txt", check.names = FALSE)
ref <- read.table("output/infercnv/infercnv.references.txt",   check.names = FALSE)
cnv.all <- as.matrix(cbind(obs, ref))                 # 基因 × 細胞；已平滑，以參考為中心（≈1）

cnv.score <- colMeans((cnv.all - 1)^2)                # ① 偏離參考的程度
top       <- names(sort(cnv.score, decreasing = TRUE))[1:200]
mal.prof  <- rowMeans(cnv.all[, top])                 #   「最像惡性」的 200 顆平均 profile
cnv.cor   <- apply(cnv.all, 2, cor, y = mal.prof)     # ② 與惡性 profile 的相關（Tirosh 2016）

gbm4$cnv.score <- cnv.score[colnames(gbm4)]
gbm4$cnv.cor   <- cnv.cor[colnames(gbm4)]
p <- ggplot(gbm4@meta.data, aes(cnv.score, cnv.cor, colour = celltype_author)) +
     geom_point(size = .7, alpha = .7) + theme_classic() +
     labs(x = "CNV score (deviation from reference)", y = "CNV correlation (with malignant profile)")
ggsave("output/figs/05_cnv_scatter.png", p, width = 8, height = 6, dpi = 150, bg = "white")
# 看圖：右上（高、高）= 惡性；左下 = 正常；中間 = 不確定。閾值畫在兩群之間的谷。

## ---- 3. triangulate ------------------------------------------------ Q3 頁 25–26
# ★ 閾值請看你的散點圖後調整；以下數字只是常見起點 ★
s.hi <- quantile(gbm4$cnv.score[gbm4$celltype_author %in% refs], 0.99)   # 參考組的 99 分位當「正常上限」
c.hi <- 0.4
gbm4$malignant <- with(gbm4@meta.data, ifelse(
  cnv.score > s.hi & cnv.cor > c.hi, "malignant",
  ifelse(cnv.score <= s.hi & cnv.cor < 0.2, "normal", "unresolved")))

# 第三個證人：譜系。免疫／血管／寡樹突被標惡性 → 一律改 unresolved（多半是 doublet 或雜訊）
lineage.normal <- c("Immune cell", "Vascular", "Oligodendrocyte")
gbm4$malignant[gbm4$malignant == "malignant" & gbm4$celltype_author %in% lineage.normal] <- "unresolved"

# 第二個證人：按病人成島。惡性細胞應集中在 raw_clusters 中「單一病人為主」的群
tab <- table(gbm4$malignant, gbm4$celltype_author)
print(tab)
agree <- mean((gbm4$malignant == "malignant") == (gbm4$celltype_author == "Neoplastic"),
              na.rm = TRUE)
cat(sprintf("與作者 Neoplastic 標籤一致率（含 unresolved 視為不一致）：%.1f%%\n", 100 * agree))
# 不一致的細胞在哪？常是 Periphery 的浸潤細胞（CNV 訊號弱）
print(table(gbm4$malignant, gbm4$tissue))

p <- DimPlot(gbm4, reduction = "umap.raw", group.by = "malignant", cols = c("#C0392B", "#1A6B5A", "grey70"))
ggsave("output/figs/05_umap_malignant.png", p, width = 7, height = 6, dpi = 150, bg = "white")
saveRDS(gbm4, "output/gbm4_malignant.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 5
#  5-1 【雷區實驗】把 refs 改成 c("Astrocyte")，重跑 inferCNV（可只跑 §1）。熱圖變成什麼樣？
#      chr7 / chr10 的訊號還在嗎？用一句話解釋為什麼。
#  5-2 把 top 從 200 改成 50 與 500，cnv.cor 的分布變了多少？判定結果（malignant 數）差幾顆？
#  5-3 只看 Periphery 的細胞：作者標 Neoplastic、但你標 unresolved 的有幾顆？
#      它們的 cnv.score 分布跟 Tumor 的 Neoplastic 比如何？這告訴你浸潤細胞的什麼特性？
#  5-4 對每位病人，畫出該病人惡性細胞在 chr7 與 chr10 的平均 CNV 值（cnv.all 的列名含基因，
#      配合 gpos 找出染色體）。四位病人都有 chr7+/chr10− 嗎？
#  進階 開 HMM = TRUE 重跑其中一位病人，看 infercnv 的 subcluster 結果：這位病人有幾個亞株？
#      各亞株的私有事件是什麼？
# =====================================================================
