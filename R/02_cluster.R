# =====================================================================
# 02_cluster.R — 練習腳本 2：前處理、PCA、分群、UMAP 與三個必要檢查
#
# 對應影片：Q2 頁 25–36（§1 四行 + 週期、§2 PC1 與 nPC、§3 分群、掃描與穩定性檢查、§4 每群 QC、§4b doublet 群診斷）
# 輸入：output/rds/01_gbm_qc.rds（01_qc.R）
# 輸出：output/rds/02_gbm_clustered.rds、output/tables/02_per_cluster_qc.csv、output/figs/02_*.png
# 時間：約 3–5 分鐘
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm <- readRDS("output/rds/01_gbm_qc.rds")
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)   # 圖檔資料夾先建好，後面 ggsave 才不會失敗

## ---- 1. preprocess ------------------------------------------------- Q2 頁 25–27
gbm <- NormalizeData(gbm)                             # CP10K + log1p：除掉深度
gbm <- FindVariableFeatures(gbm, nfeatures = 2000)    # 2,000 個 HVG
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
gbm <- RunPCA(gbm, npcs = 50, verbose = FALSE)

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
gbm <- FindClusters(gbm, resolution = 0.5)            # 起點，不是答案
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
gbm <- FindClusters(gbm, resolution = 0.5, random.seed = 42, cluster.name = "seed42")
mclust::adjustedRandIndex(gbm$RNA_snn_res.0.5, gbm$seed42)     # 接近 1 = 穩
gbm <- FindNeighbors(gbm, dims = 1:npc, k.param = 30)
gbm <- FindClusters(gbm, resolution = 0.5, cluster.name = "k30")
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
## >>> 參考答案 ------------------------------------------------------
# 本課這份 GBM 5k 的結果（13 群，res 0.5 / k = 20；doublet 用 DoubletFinder，全場 198 顆）：
#   cluster 12：dbl = 0.60、nF = 6,221（全場最高）→ 唯一超過 0.5 的群；§4b 三項診斷全中，判定為真 doublet 群
#   cluster  9：dbl = 0.31、nF = 5,370 → 偏高但未過線，註釋時留意
#   cluster  7：dbl = 0.20，但 cyc = 1.00 → 這是增殖群。增殖細胞轉錄體複雜、RNA 量大，
#               本來就容易被誤判成 doublet，這 20% 多半是誤判而非真 doublet（03 標成 Cycling，不移除）
#   cluster  1：mt = 1.7（最低）、cluster 10：nF = 878（最低）→ 寡樹突與 T 細胞的典型輪廓
#   其餘九群 dbl ≤ 0.02。
# 沒有任何一群符合「percent.mt 明顯偏高且 nFeature 最低」，代表 01 的 QC 線切得乾淨。
#
# 值得注意的是判定的「集中程度」：cluster 12 只佔全部細胞 2.5%（131/5,261），
#   卻收下約四成的 doublet 判定（0.60 × 131 ≈ 79 顆 / 198 顆），約 16 倍富集。
#   判定越集中在少數群，「整群 doublet」的推論越有依據；若 198 顆平均散在 13 群，
#   那就只是分數雜訊，不該動任何一群。§4b 的 ① 會印出確切數字：
#   cluster 12 = 79 顆、9 = 43、7 = 54、0 = 16、2 = 6，其餘為 0，合計 198。
# ※ cluster 12 在 k = 30 的分群下會被併進 AC 樣膠質群而看不見（見 §3 交叉表）——
#   這就是為什麼 seurat_clusters 一定要等穩定性檢查跑完才寫回去（見 §3 結尾）。
## <<< 參考答案

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
## >>> 參考答案 ------------------------------------------------------
# 本課這份資料的判定結果（實跑驗證，seed = 1234、npc = 25、res 0.5、k = 20）：
# cluster 12 三項全中，是真 doublet 群。
#   ① 198 顆 doublet 判定裡有 79 顆（40%）集中在這一群，而它只佔 2.5% 的細胞
#      （131 / 5,261）——約 16 倍富集。判定越集中，「整群」的推論越站得住
#   ② nCount 中位數 31,638，全場最高。膠質 cluster 2（22,583）＋ 髓系 cluster 4（9,244）
#      ＝ 31,827，兩者只差 0.6%——正是兩顆細胞的 RNA 加在一起的樣子
#   ③ 同一顆細胞同時表現 PTPRC 與 SOX2 的比例 0.649，全體只有 0.038，約 17 倍。
#      免疫與膠質是互斥譜系，一顆正常細胞不會兩邊都開
#   ④ top 15 marker 全是髓系的標準基因（MS4A7、MS4A4A、PLEK、LCP2、BCL2A1、GPR183…），
#      沒有一個是它專屬的；avg_log2FC 最高只有 1.27，是髓系訊號被稀釋後的樣子
#   → 03_annotate.R 會依 output/tables/02_per_cluster_qc.csv 自動把它標成 DOUBLET 並整群移除，
#      5,261 顆剩 5,130 顆。方法段要寫：「移除 1 個整群 doublet 群集（131 顆，佔 2.5%）」。
#
# 兩個沒過線、但會被問到的群：
#   cluster 9（dbl = 0.31）：COL3A1 陽性率 0.89，是血管周／纖維母細胞。這類細胞稀少且 RNA 量大，
#     容易拿到偏高的 doublet 分數；但 nCount 沒有出現「兩群相加」、marker 專屬性明確 → 保留。
#   cluster 7（dbl = 0.20、cyc = 1.00）：增殖群，理由同上 → 保留，標成 Cycling。
## <<< 參考答案

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
