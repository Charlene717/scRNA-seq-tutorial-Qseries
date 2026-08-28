# =====================================================================
# 02_cluster.R — 練習腳本 2：前處理、PCA、分群、UMAP 與三個保命檢查
#
# 對應影片：Q2 頁 25–36（§1 四行 + 週期、§2 PC1 與 nPC、§3 分群、掃描與穩定性檢查、§4 每群 QC）
# 輸入：output/gbm_qc.rds（01_qc.R）
# 輸出：output/gbm_clustered.rds
# 時間：約 3–5 分鐘
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm <- readRDS("output/gbm_qc.rds")

## ---- 1. preprocess ------------------------------------------------- Q2 頁 25–27
gbm <- NormalizeData(gbm)                             # CP10K + log1p：除掉深度
## TODO ▶ 挑幾個高變異基因？（Q2 頁 25）
gbm <- FindVariableFeatures(gbm, nfeatures = ____)
head(VariableFeatures(gbm), 20)                       # 先看名單：惡性 / 免疫 / 寡樹突的基因都在嗎？
p <- LabelPoints(VariableFeaturePlot(gbm), points = head(VariableFeatures(gbm), 12), repel = TRUE)
ggsave("output/figs/02_hvg.png", p, width = 8, height = 5, dpi = 150, bg = "white")

# 週期分數：腫瘤裡增殖細胞多，先算好，分群後才能檢查它有沒有主導分群
gbm <- CellCycleScoring(gbm, s.features   = cc.genes.updated.2019$s.genes,
                             g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
table(gbm$Phase)

gbm <- ScaleData(gbm)                                 # z-score、clip ±10；預設只做 HVG
## TODO ▶ 先算幾個 PC 存著？（Q2 頁 25）
gbm <- RunPCA(gbm, npcs = ____, verbose = FALSE)

## ---- 2. pca-checks ------------------------------------------------- Q2 頁 28–31
# 保命檢查一：PC1 是生物學，不是深度
print(gbm[["pca"]], dims = 1:3, nfeatures = 8)        # GBM：PC1 幾乎一定是免疫 vs 膠質
r1 <- cor(Embeddings(gbm, "pca")[, 1], gbm$nCount_RNA)
cat(sprintf("cor(PC1, nCount) = %.2f  （|r| < 0.3 放心；0.3–0.5 常見於惡性 RNA 量大的腫瘤，PC1 載荷是生物學即可；> 0.5 要查）\n", r1))
if (abs(r1) > 0.5) warning("PC1 與深度高度相關：回去確認 NormalizeData 有跑、且跑在 counts 上")

# 保命檢查二：取幾個 PC
p <- ElbowPlot(gbm, ndims = 50)
ggsave("output/figs/02_elbow.png", p, width = 7, height = 4, dpi = 150, bg = "white")
pct  <- gbm[["pca"]]@stdev / sum(gbm[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]                 # 累積 90% 且單 PC < 5%
co2  <- sort(which(diff(pct) < -0.1), decreasing = TRUE)[1] + 1
cat("量化 cutoff：", co1, co2, "\n")
npc <- 25                                             # 腫瘤結構複雜：20–30；取多一點比取少安全

## ---- 3. cluster-umap ----------------------------------------------- Q2 頁 33–35
gbm <- FindNeighbors(gbm, dims = 1:npc)               # PC 空間建 KNN / SNN 圖
## TODO ▶ 分群的 resolution 從掃描與 clustree 挑哪一層？（Q2 頁 33–34）
gbm <- FindClusters(gbm, resolution = ____)            # 起點，不是答案
gbm <- RunUMAP(gbm, dims = 1:npc)                     # 同一組 dims

# resolution 掃描 + clustree：找穩定層
for (r in seq(0.2, 1.2, by = 0.2))
  gbm <- FindClusters(gbm, resolution = r, verbose = FALSE)
library(clustree)
p <- clustree(gbm, prefix = "RNA_snn_res.")
ggsave("output/figs/02_clustree.png", p, width = 9, height = 8, dpi = 150, bg = "white")
Idents(gbm) <- "RNA_snn_res.0.5"                      # 掃完一定指回定案欄位！
gbm$seurat_clusters <- gbm$RNA_snn_res.0.5

p <- DimPlot(gbm, label = TRUE) + DimPlot(gbm, group.by = "Phase")
ggsave("output/figs/02_umap.png", p, width = 12, height = 5, dpi = 150, bg = "white")

# 穩定性三檢查（Q2 頁 35）：換 seed、換 k.param、看群大小
## TODO ▶ 分群的 resolution 從掃描與 clustree 挑哪一層？（Q2 頁 33–34）
gbm <- FindClusters(gbm, resolution = ____, random.seed = 42, cluster.name = "seed42")
mclust::adjustedRandIndex(gbm$RNA_snn_res.0.5, gbm$seed42)     # 接近 1 = 穩
gbm <- FindNeighbors(gbm, dims = 1:npc, k.param = 30)
## TODO ▶ 分群的 resolution 從掃描與 clustree 挑哪一層？（Q2 頁 33–34）
gbm <- FindClusters(gbm, resolution = ____, cluster.name = "k30")
table(gbm$RNA_snn_res.0.5, gbm$k30)                              # 對角線清楚 = 穩
sort(table(gbm$RNA_snn_res.0.5))                                 # < 30 顆的群先存疑
gbm <- FindNeighbors(gbm, dims = 1:npc)                         # 換回預設 k = 20
Idents(gbm) <- "RNA_snn_res.0.5"

## ---- 4. per-cluster-qc --------------------------------------------- Q2 頁 36
# 保命檢查三：分群後回頭看每群的 QC、doublet、週期
qc.tab <- gbm@meta.data |>
  group_by(cluster = seurat_clusters) |>
  summarise(n   = n(),
            mt  = round(median(percent.mt), 1),
            nF  = round(median(nFeature_RNA)),
            dbl = round(mean(dbl.class == "doublet"), 2),
            cyc = round(mean(Phase != "G1"), 2)) |>
  arrange(desc(dbl))
print(qc.tab, n = 30)
write.csv(qc.tab, "output/02_per_cluster_qc.csv", row.names = FALSE)
# 判讀：
#   dbl > 0.5 且卡在兩大群之間 → 整群是 doublet 產物 → 03_annotate.R 移除並記錄
#   mt 一枝獨秀 + nF 墊底       → QC 殘渣 → 回 01_qc.R 修線
#   cyc 很高                    → 增殖狀態，可能混多種型別 → 註釋時小心

p <- VlnPlot(gbm, features = c("percent.mt", "nFeature_RNA", "dbl.score"), pt.size = 0, ncol = 3)
ggsave("output/figs/02_per_cluster_qc.png", p, width = 13, height = 4, dpi = 150, bg = "white")

## ---- 5. save -------------------------------------------------------
saveRDS(gbm, "output/gbm_clustered.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 2
#  2-1 把 npc 改成 10 與 40，各跑一次分群（res 0.5）。用 table(npc10, npc25) 交叉表比較：
#      哪些群在 npc = 10 時被合併了？它們是什麼細胞？（提示：03 的 marker）
#  2-2 在 clustree 上找出「穩定層」與「亂流層」各是哪幾個 resolution。
#      挑相鄰兩個 resolution 之間多切的一刀，用 FindMarkers(ident.1, ident.2) 看兩邊差在哪些基因，
#      判斷這一刀該不該切。
#  2-3 用 ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) 重跑 PCA 與分群。
#      增殖細胞還自成一群嗎？哪種做法更適合腫瘤資料？寫下你的理由。
#  2-4 per-cluster QC 表裡 dbl 最高的群，在 UMAP 上位於哪裡？它介於哪兩群之間？
#  進階 用 FeaturePlot(features = "nCount_RNA") 看深度在 UMAP 上的分布；
#      若某群明顯偏深或偏淺，那是生物學（惡性細胞 RNA 多）還是技術？怎麼分辨？
# =====================================================================
