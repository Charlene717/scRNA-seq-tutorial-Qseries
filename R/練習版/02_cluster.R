# =====================================================================
# 02_cluster.R — 練習腳本 2：前處理、PCA、分群、UMAP 與三個必要檢查
#
# 對應影片：Q2 頁 25–36（§1 四行 + 週期、§2 PC1 與 nPC、§3 分群、掃描與穩定性檢查、§4 每群 QC、§4b doublet 群診斷）
# 輸入：output/rds/01_gbm_qc.rds（01_qc.R）
# 輸出：output/rds/02_gbm_clustered.rds、output/tables/02_per_cluster_qc.csv、output/figs/02_*.png
# 時間：約 3–5 分鐘
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm <- readRDS("output/rds/01_gbm_qc.rds")
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)   # 圖檔資料夾先建好，後面 ggsave 才不會失敗

## ---- 1. preprocess ------------------------------------------------- Q2 頁 25–27
gbm <- NormalizeData(gbm)                             # CP10K + log1p：除掉深度
## TODO ▶ 挑幾個高變異基因？（Q2 頁 25）
gbm <- FindVariableFeatures(gbm, nfeatures = ____)
head(VariableFeatures(gbm), 20)                       # 先看名單，解讀見下
# HVG 排的是「變異量」，不是「重要性」：稀有但轉錄體極端的族群會把前段名次吃掉。
# 這份資料前 20 名幾乎全是基質／血管基因（COL3A1、COL1A1、DCN、LUM、POSTN…），
# 來自只有 140 顆的血管周／纖維母細胞群（跑完 §3 就是 cluster 9，佔 2.7%）。
# 這不代表資料裡沒有免疫或膠質——它們在 2,000 個 HVG 裡面，只是排在後面。
# 所以不要只看前 20 名就下結論，用下面這行直接確認關鍵譜系有沒有入選：
lin <- c("PTPRC", "SOX2", "MBP", "PECAM1", "COL3A1")  # 免疫／膠質／寡樹突／血管／基質
print(setNames(lin %in% VariableFeatures(gbm), lin))
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
# 必要檢查一：PC1 是生物學，不是深度
print(gbm[["pca"]], dims = 1:3, nfeatures = 8)        # GBM：PC1 幾乎一定是免疫 vs 膠質
r1 <- cor(Embeddings(gbm, "pca")[, 1], gbm$nCount_RNA)
cat(sprintf("cor(PC1, nCount) = %.2f  （|r| < 0.3 放心；0.3–0.5 常見於惡性 RNA 量大的腫瘤，PC1 載荷是生物學即可；> 0.5 要查）\n", r1))
if (abs(r1) > 0.5) warning("PC1 與深度高度相關：回去確認 NormalizeData 有跑、且跑在 counts 上")

# 必要檢查二：取幾個 PC
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
# ⚠ 這裡先「不要」寫 gbm$seurat_clusters：下面 §3 的穩定性檢查還會再呼叫 FindClusters()，
#   而 FindClusters() 不論有沒有給 cluster.name，都會順手覆寫 seurat_clusters 與 Idents。
#   定案欄位統一等所有檢查跑完之後（本節結尾）才寫回去，否則 03 會拿到錯的分群。

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

# ★ 定案：穩定性檢查全部跑完，才把 Idents 與 seurat_clusters 一起指回 res 0.5 ★
Idents(gbm) <- "RNA_snn_res.0.5"
gbm$seurat_clusters <- gbm$RNA_snn_res.0.5
cat("定案分群（npc =", npc, "、res 0.5、k = 20）共", nlevels(factor(gbm$seurat_clusters)), "群\n")

## ---- 4. per-cluster-qc --------------------------------------------- Q2 頁 36
# 必要檢查三：分群後回頭看每群的 QC、doublet、週期
qc.tab <- gbm@meta.data |>
  group_by(cluster = seurat_clusters) |>
  summarise(n   = n(),
            mt  = round(median(percent.mt), 1),
            nF  = round(median(nFeature_RNA)),
            dbl = round(mean(dbl.class == "doublet"), 2),
            cyc = round(mean(Phase != "G1"), 2)) |>
  arrange(desc(dbl))
print(qc.tab, n = 30)
write.csv(qc.tab, "output/tables/02_per_cluster_qc.csv", row.names = FALSE)
# 判讀：
#   dbl > 0.5 且卡在兩大群之間 → 整群是 doublet 產物 → 03_annotate.R 移除並記錄
#   percent.mt 明顯偏高且 nFeature 最低       → QC 殘渣 → 回 01_qc.R 調整閾值
#   cyc 很高                    → 增殖狀態，可能混多種型別 → 註釋時小心
#
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）

p <- VlnPlot(gbm, features = c("percent.mt", "nFeature_RNA", "dbl.score"), pt.size = 0, ncol = 3)
ggsave("output/figs/02_per_cluster_qc.png", p, width = 13, height = 4, dpi = 150, bg = "white")

## ---- 4b. doublet-cluster diagnostics -------------------------------- Q2 頁 36
# dbl > 0.5 的群是「整群 doublet」的候選——但整群移除是不可逆的決定，不能只憑這個數字。
# 腫瘤裡 RNA 量大的族群（惡性細胞、增殖細胞）本來就容易被誤判，所以三項證據一起看。
cand <- as.character(qc.tab$cluster[qc.tab$dbl > 0.5])
cat("整群 doublet 候選（dbl > 0.5）：", if (length(cand)) paste(cand, collapse = ", ") else "無", "\n")

if (length(cand)) {
  # ① 精確的每群 doublet 數（qc.tab 的 dbl 是四捨五入後的比例）
  print(table(cluster = gbm$seurat_clusters, gbm$dbl.class))

  # ② 深度：真 doublet 帶著兩顆細胞的 RNA，nCount 中位數應接近「兩個來源群相加」
  #    （同型 doublet 才是單一群的兩倍；異型 doublet 要拿兩個親代群的中位數相加來比）
  dep <- aggregate(cbind(nCount_RNA, nFeature_RNA) ~ seurat_clusters, gbm@meta.data, median)
  print(dep[order(-dep$nCount_RNA), ])

  # ③ 互斥譜系是否同時出現在同一顆細胞——真 doublet 的指紋。用數字看，不要只靠目測。
  lin <- c("PTPRC", "SOX2", "MBP", "PECAM1", "COL3A1")   # 免疫／膠質／寡樹突／血管／基質
  d   <- FetchData(gbm, c("seurat_clusters", lin))
  cat("\n各群表現各譜系標誌的細胞比例：\n")
  print(aggregate(d[, lin], list(cluster = d$seurat_clusters), function(x) round(mean(x > 0), 2)))
  cat("\n候選群「同一顆細胞同時表現 PTPRC 與 SOX2」的比例（全體平均可當基準）：\n")
  base <- round(mean(d$PTPRC > 0 & d$SOX2 > 0), 3)
  for (cl in cand) {
    i <- d$seurat_clusters == cl
    cat(sprintf("  cluster %-3s %.3f   （全體 %.3f）\n", cl,
                mean(d$PTPRC[i] > 0 & d$SOX2[i] > 0), base))
  }
  p <- VlnPlot(gbm, lin, group.by = "seurat_clusters", pt.size = 0, ncol = 5)
  ggsave("output/figs/02_doublet_lineage_check.png", p, width = 16, height = 4, dpi = 150, bg = "white")

  # ④ 候選群有沒有「專屬」marker，還是只是另外兩群 marker 的聯集
  for (cl in cand) {
    cat("\n--- cluster", cl, "的 top marker ---\n")
    print(head(FindMarkers(gbm, ident.1 = cl, only.pos = TRUE,
                           min.pct = 0.25, logfc.threshold = 0.5), 15))
  }
}
# 判讀：
#   ② 接近兩倍 ＋ ③ 同一群同時亮互斥譜系 ＋ ④ 沒有專屬 marker → 真 doublet，03 整群移除並記錄。
#   ④ 有明確的專屬 marker、且 ② 沒有翻倍                     → 被工具誤判的真實族群，不要刪。
# 反例就在同一張表裡：cluster 7 的 dbl = 0.20 但 cyc = 1.00，那是增殖群——
#   增殖細胞轉錄體複雜、RNA 量大，天生就容易被誤判成 doublet，那 20% 多半不是真 doublet。
#
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）

## ---- 5. save -------------------------------------------------------
saveRDS(gbm, "output/rds/02_gbm_clustered.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 2
#  2-1 把 npc 改成 10 與 40，各跑一次分群（res 0.5）。用 table(npc10, npc25) 交叉表比較：
#      哪些群在 npc = 10 時被合併了？它們是什麼細胞？（提示：03 的 marker）
#  2-2 在 clustree 上找出「穩定層」與「不穩定層」各是哪幾個 resolution。
#      挑相鄰兩個 resolution 之間多切的一刀，用 FindMarkers(ident.1, ident.2) 看兩邊差在哪些基因，
#      判斷這一刀該不該切。
#  2-3 用 ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) 重跑 PCA 與分群。
#      增殖細胞還自成一群嗎？哪種做法更適合腫瘤資料？寫下你的理由。
#  2-4 per-cluster QC 表裡 dbl 最高的群，在 UMAP 上位於哪裡？它介於哪兩群之間？
#  進階 用 FeaturePlot(features = "nCount_RNA") 看深度在 UMAP 上的分布；
#      若某群明顯偏深或偏淺，那是生物學（惡性細胞 RNA 多）還是技術？怎麼分辨？
# =====================================================================
