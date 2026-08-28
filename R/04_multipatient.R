# =====================================================================
# 04_multipatient.R — 練習腳本 4：四位病人（GSE84465）載入、基線、Seurat 整合（預設）與 Harmony（比較）、正負對照檢查
#
# 對應影片：Q3 頁 8–17（§1 載入、§2 建物件與 QC、§3 基線、§4 整合：Seurat CCA 預設 + Harmony 比較、§5 正負對照與 LISI）
# 輸入：data/GSE84465_GBM_All_data.csv.gz（00_setup.R）；metadata 由 GEOquery 下載
# 輸出：output/gbm4_unintegrated.rds、output/gbm4_integrated.rds（含 integrated.cca 與 harmony 兩個 reduction）
# 時間：約 5–10 分鐘（getGEO 需網路）
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)

## ---- 1. load ------------------------------------------------------- Q3 頁 8
# 注意：這個 csv.gz 其實是「空白分隔」（GEO 常見），read.csv 會讀成 0 欄；用 read.table + sep = " "
cnt <- read.table(gzfile("data/GSE84465_GBM_All_data.csv.gz"), header = TRUE, sep = " ",
                  row.names = 1, check.names = FALSE, quote = "\"")
dim(cnt)                                              # 23,460 × 3,589（基因 × 細胞）
stopifnot("欄數為 0：分隔符號不對，改試 sep = \"\\t\" 或 \",\"" = ncol(cnt) > 0)
tail(rownames(cnt))                                   # HTSeq 的統計列不是基因，要拿掉
bad <- c("no_feature", "ambiguous", "alignment_not_unique", "too_low_aQual", "not_aligned")
cnt <- cnt[!rownames(cnt) %in% bad, ]
cnt <- as(as.matrix(cnt), "dgCMatrix")                  # 轉稀疏矩陣，省記憶體

library(GEOquery)
gse  <- getGEO("GSE84465", GSEMatrix = TRUE)[[1]]     # series matrix（含每顆細胞的 characteristics）
meta <- pData(gse)
# characteristics 欄位順序不保證 → 用內容自動找欄位
find.col <- function(df, key) {
  hit <- sapply(df, function(v) any(grepl(paste0("^", key, ":"), v)))
  stopifnot(any(hit)); names(df)[which(hit)[1]]
}
get.val <- function(df, key) sub(paste0("^", key, ": *"), "", df[[find.col(df, key)]])
meta$plate   <- get.val(meta, "plate id")
meta$well    <- get.val(meta, "well")
meta$patient <- get.val(meta, "patient id")           # BT_S1 / BT_S2 / BT_S4 / BT_S6
meta$tissue  <- get.val(meta, "tissue")               # Tumor / Periphery
meta$celltype_author <- get.val(meta, "cell type")    # Neoplastic / Astrocyte / OPC / ...
meta$cell    <- paste(meta$plate, meta$well, sep = ".")   # 與 counts 欄名對齊：例 1001000173.G8
head(meta[, c("cell", "patient", "tissue", "celltype_author")])

## ---- 2. object-qc -------------------------------------------------- Q3 頁 9
common <- intersect(colnames(cnt), meta$cell)
cat("對得上的細胞：", length(common), "/", ncol(cnt), "\n")   # 應接近全部；對不上先檢查命名
cnt  <- cnt[, common]
meta <- meta[match(common, meta$cell), ]; rownames(meta) <- common

gbm4 <- CreateSeuratObject(counts = cnt,
                           meta.data = meta[, c("patient", "tissue", "celltype_author", "plate")],
                           project = "GSE84465", min.cells = 3, min.features = 500)   # Smart-seq2：下限 500
gbm4[["percent.mt"]] <- PercentageFeatureSet(gbm4, pattern = "^MT-")
table(gbm4$patient, gbm4$tissue)                      # 4 × 2 = 8 個樣本
table(gbm4$celltype_author)
p <- VlnPlot(gbm4, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "patient", pt.size = 0, ncol = 3)
ggsave("output/figs/04_qc_by_patient.png", p, width = 13, height = 4, dpi = 150, bg = "white")
# 作者已 QC；這裡若有極端離群（percent.mt > median + 3×MAD）可再濾，並記錄。

## ---- 3. baseline --------------------------------------------------- Q3 頁 11–12
# 先跑一次「不整合」的流程，留基線；之後惡性細胞的分析回來用這一份
gbm4 <- NormalizeData(gbm4) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(verbose = FALSE)
gbm4 <- FindNeighbors(gbm4, dims = 1:25) |> FindClusters(resolution = 0.5, cluster.name = "raw_clusters")
gbm4 <- RunUMAP(gbm4, dims = 1:25, reduction.name = "umap.raw")
p <- DimPlot(gbm4, reduction = "umap.raw", group.by = "patient") +
     DimPlot(gbm4, reduction = "umap.raw", group.by = "celltype_author", label = TRUE, repel = TRUE)
ggsave("output/figs/04_umap_unintegrated.png", p, width = 13, height = 5.5, dpi = 150, bg = "white")
# 看圖：Neoplastic 是否按病人分島？Immune / Oligodendrocyte 是否跨病人混合？
saveRDS(gbm4, "output/gbm4_unintegrated.rds")

# 量化「按病人分島」：每群裡病人的熵（越低越單一病人）
patient.entropy <- gbm4@meta.data |> dplyr::count(raw_clusters, patient) |>
  group_by(raw_clusters) |> mutate(p = n / sum(n)) |>
  summarise(entropy = -sum(p * log(p)), n = sum(n), top_type = NA)
print(patient.entropy)

## ---- 4. integrate ------------------------------------------------- Q3 頁 14–15
# 預設走 Seurat 自己的整合（CCA 錨點；大資料可改 RPCAIntegration）。Harmony 是「Seurat 整合不理想時」的替代方案，
# 本練習兩個都跑、放在同一個物件的兩個 reduction 裡，用同一套檢查比較。
gbm4[["RNA"]] <- split(gbm4[["RNA"]], f = gbm4$patient)   # 每位病人一個 layer（批次 = 病人）
gbm4 <- NormalizeData(gbm4) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(verbose = FALSE)

# (a) Seurat 預設：CCA 錨點整合 → reduction "integrated.cca"
gbm4 <- IntegrateLayers(gbm4, method = CCAIntegration, orig.reduction = "pca",
                        new.reduction = "integrated.cca", verbose = FALSE)
# (b) 比較用：Harmony → reduction "harmony"（需 harmony 套件）
gbm4 <- IntegrateLayers(gbm4, method = HarmonyIntegration, orig.reduction = "pca",
                        new.reduction = "harmony", verbose = FALSE)
gbm4 <- JoinLayers(gbm4)

# 兩個整合空間各自分群 + UMAP（同一組 dims / resolution 才能比較）
for (red in c("integrated.cca", "harmony")) {
  tag <- if (red == "harmony") "harmony" else "cca"
  gbm4 <- FindNeighbors(gbm4, reduction = red, dims = 1:25, graph.name = paste0(tag, "_snn")) |>
          FindClusters(graph.name = paste0(tag, "_snn"), resolution = 0.5, cluster.name = paste0(tag, "_clusters")) |>
          RunUMAP(reduction = red, dims = 1:25, reduction.name = paste0("umap.", tag))
}
p <- (DimPlot(gbm4, reduction = "umap.raw",     group.by = "patient") + ggtitle("Unintegrated")) +
     (DimPlot(gbm4, reduction = "umap.cca",     group.by = "patient") + ggtitle("Seurat CCA (default)")) +
     (DimPlot(gbm4, reduction = "umap.harmony", group.by = "patient") + ggtitle("Harmony"))
ggsave("output/figs/04_umap_integration_compare.png", p, width = 18, height = 5.5, dpi = 150, bg = "white")
p <- DimPlot(gbm4, reduction = "umap.cca", group.by = "celltype_author", label = TRUE, repel = TRUE) +
     DimPlot(gbm4, reduction = "umap.harmony", group.by = "celltype_author", label = TRUE, repel = TRUE)
ggsave("output/figs/04_umap_integrated_celltype.png", p, width = 13, height = 5.5, dpi = 150, bg = "white")

# 正負對照檢查（integration positive / negative controls）：
#   正對照 = 免疫、寡樹突細胞：應跨病人混合（熵高）→ 整合成功
#   負對照 = 惡性細胞：應保留病人差異（熵低）→ 沒有過度修正
patient_entropy <- function(obj, cluster.col) {
  obj@meta.data |>
    dplyr::filter(celltype_author %in% c("Immune cell", "Oligodendrocyte", "Neoplastic")) |>
    dplyr::count(celltype_author, cl = .data[[cluster.col]], patient) |>
    group_by(celltype_author, cl) |> mutate(p = n / sum(n)) |>
    summarise(patient_entropy = -sum(p * log(p)), n = sum(n), .groups = "drop") |>
    group_by(celltype_author) |> summarise(mean_entropy = round(weighted.mean(patient_entropy, n), 2), .groups = "drop")
}
ctrl <- bind_rows(Unintegrated = patient_entropy(gbm4, "raw_clusters"),
                  Seurat_CCA = patient_entropy(gbm4, "cca_clusters"),
                  Harmony = patient_entropy(gbm4, "harmony_clusters"), .id = "method")
print(tidyr::pivot_wider(ctrl, names_from = celltype_author, values_from = mean_entropy))
# 讀法：Immune / Oligodendrocyte 的熵 整合後 > 未整合 = 正對照通過；
#       Neoplastic 的熵若也被拉到跟免疫一樣高 = 負對照失敗（過度修正）。兩種方法都會把惡性細胞混掉一部分——這正是「兩難」。
# 決定：Seurat CCA 若正對照通過、負對照可接受，就用它；只有 CCA 明顯混不好（免疫仍按病人分島）才改用 Harmony。

# 整合原則：整合空間只用來註釋正常細胞、對齊免疫細胞；惡性細胞回未整合空間、按病人看。
saveRDS(gbm4, "output/gbm4_integrated.rds")

## ---- 5. integration-metrics ---------------------------------------- Q3 頁 16
# 把「混得好不好」變成數字：每群病人組成 + LISI，兩個整合空間各算一次
round(prop.table(table(gbm4$cca_clusters, gbm4$patient), 1), 2)          # 正常細胞的群：四位病人都該有
library(lisi)                                     # remotes::install_github("immunogenomics/lisi")
for (red in c("pca", "integrated.cca", "harmony")) {
  emb <- Embeddings(gbm4, red)[, 1:25]
  gbm4[[paste0("lisi_", sub("integrated.", "", red, fixed = TRUE))]] <- compute_lisi(emb, gbm4@meta.data, "patient")$patient
}
gbm4@meta.data |> group_by(celltype_author) |>
  summarise(across(starts_with("lisi_"), ~ round(median(.x), 2))) |> print()
# 免疫／寡樹突的 LISI 接近病人數（4）= 整合成功；Neoplastic 接近 1 = 病人差異保住了（正確，不是失敗）
# 比較 lisi_cca 與 lisi_harmony：正常細胞誰混得好、惡性細胞誰保得住——這就是選整合方法的依據。
saveRDS(gbm4, "output/gbm4_integrated.rds")

# =====================================================================
# ▶ 練習 4
#  4-1 在 umap.raw 上用 split.by = "patient" 畫 celltype_author：哪些型別每位病人都有？哪些只在某些病人？
#  4-2 用 §4 的正負對照熵：Neoplastic 群的病人熵是多少？Immune 群呢？用一句話解釋差異的生物學原因。
#  4-3 再加 RPCAIntegration（只改 method），與 CCA、Harmony 三者的正負對照熵與 LISI。
#      哪一個對惡性細胞「修得最少」？在腫瘤資料裡這是優點還是缺點？
#  4-4 只取 Immune cell 子集，比較整合前後的分群：整合後能否分出小膠質 vs 巨噬細胞？
#  4-5 §5 的 LISI：Neoplastic 與 Immune 的中位數各是多少？把 Harmony 的 theta 調到 4 重跑，Neoplastic 的 LISI 變成多少？這代表什麼？
#      （提示：P2RY12 / TMEM119 vs CD163 / LYZ）
#  進階 讀 Neftel 2019 的 Methods：他們如何處理病人效應？跟本腳本的原則一致嗎？
# =====================================================================
