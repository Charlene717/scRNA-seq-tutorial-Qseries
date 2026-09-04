# =====================================================================
# 04_multipatient.R — 練習腳本 4：四位病人（GSE84465）載入、基線、Seurat 整合（預設）與 Harmony（比較）、正負對照檢查
#
# 對應影片：Q3 頁 8–17（§1 載入、§2 建物件與 QC、§3 基線、§4 整合：Seurat CCA 預設 + Harmony 比較、§5 正負對照與 LISI）
# 輸入：data/GSE84465_GBM_All_data.csv.gz（00_setup.R）；metadata 由 GEO series matrix 下載並自行解析
# 輸出：output/rds/04_gbm4_unintegrated.rds、output/rds/04_gbm4_integrated.rds（含 integrated.cca 與 harmony 兩個 reduction）
#       output/tables/04_patient_entropy_raw.csv、04_integration_controls.csv、04_lisi_by_celltype.csv
#       output/figs/04_qc_by_patient.png、04_umap_unintegrated.png、04_umap_integration_compare.png、04_umap_integrated_celltype.png
# 時間：約 5–10 分鐘（下載 series matrix 需網路）
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. load ------------------------------------------------------- Q3 頁 8
# 注意：這個 csv.gz 其實是「空白分隔」（GEO 常見），read.csv 會讀成 0 欄；用 read.table + sep = " "
cnt <- read.table(gzfile("data/GSE84465_GBM_All_data.csv.gz"), header = TRUE, sep = " ",
                  row.names = 1, check.names = FALSE, quote = "\"")
dim(cnt)                                              # 23,465 × 3,589：含 5 列 HTSeq 統計，下一步拿掉後才是 23,460 個基因
stopifnot("欄數為 0：分隔符號不對，改試 sep = \"\\t\" 或 \",\"" = ncol(cnt) > 0)
tail(rownames(cnt))                                   # HTSeq 的統計列不是基因，要拿掉
bad <- c("no_feature", "ambiguous", "alignment_not_unique", "too_low_aQual", "not_aligned")
cnt <- cnt[!rownames(cnt) %in% bad, ]
cnt <- as(as.matrix(cnt), "dgCMatrix")                  # 轉稀疏矩陣，省記憶體

# ---- metadata：GEO series matrix（含每顆細胞的 characteristics）----
# ⚠ 不用 getGEO("GSE84465", GSEMatrix = TRUE)：它內部的 parseGSEMatrix 只要暫存檔下載不完整，
#   就會丟出 "parsing failed--expected only one '!series_data_table_begin'"（GEOquery 已知問題）。
#   改成自己下載到專案的 data/geo/（可重複使用、壞了刪掉重抓），再直接解析 !Sample_ 標頭列。
dir.create("data/geo", recursive = TRUE, showWarnings = FALSE)
sm <- "data/geo/GSE84465_series_matrix.txt.gz"
if (!file.exists(sm) || file.size(sm) < 5e4) {         # 太小 = 上次下載被截斷，重抓（正常約 142 KB、解壓 8 MB）
  options(timeout = 3600)                              # 預設 60 秒不夠，大檔會半路斷掉
  download.file(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/",
                       "matrix/GSE84465_series_matrix.txt.gz"), sm, mode = "wb")
}
hdr <- grep("^!Sample_", readLines(gzfile(sm), warn = FALSE), value = TRUE)
splitrow <- function(x) scan(text = sub("^![^\t]+\t", "", x), what = "", sep = "\t",
                             quote = "\"", quiet = TRUE)
meta <- as.data.frame(do.call(cbind, lapply(grep("^!Sample_characteristics_ch1", hdr, value = TRUE),
                                            splitrow)), stringsAsFactors = FALSE)
meta$geo_accession <- splitrow(grep("^!Sample_geo_accession", hdr, value = TRUE)[1])
stopifnot("series matrix 解析失敗：刪掉 data/geo/GSE84465_series_matrix.txt.gz 重抓一次" =
            nrow(meta) > 3000)
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
                           ## TODO ▶ 粗篩門檻：少於幾顆細胞表現的基因、少於幾個基因的細胞要丟？（Q2 頁 8）
                           project = "GSE84465", min.cells = ____, min.features = ____)
## TODO ▶ 粒線體基因的前綴是什麼？人與鼠不同（Q2 頁 9）
gbm4[["percent.mt"]] <- PercentageFeatureSet(gbm4, pattern = "____")
cat("建物件：", ncol(cnt), "顆 →", ncol(gbm4), "顆（min.features = 500 濾掉",
    ncol(cnt) - ncol(gbm4), "顆）\n")                  # 本課：3,589 → 3,539，濾掉 50 顆
table(gbm4$patient, gbm4$tissue)                      # 4 × 2 = 8 個樣本
# 注意 celltype_author 裡的 "Astocyte"：這是 GEO 原始 metadata 的拼字錯誤（少了 r）。
# 不要「順手改掉」——後面與作者標籤比對時要對得上原檔。
table(gbm4$celltype_author)
p <- VlnPlot(gbm4, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "patient", pt.size = 0, ncol = 3)
ggsave("output/figs/04_qc_by_patient.png", p, width = 13, height = 4, dpi = 150, bg = "white")
# 作者已 QC；這裡若有極端離群（percent.mt > median + 3×MAD）可再濾，並記錄。

## ---- 3. baseline --------------------------------------------------- Q3 頁 11–12
# 先跑一次「不整合」的流程，留基線；之後惡性細胞的分析回來用這一份
## TODO ▶ 挑幾個高變異基因？（Q2 頁 25）
gbm4 <- NormalizeData(gbm4) |> FindVariableFeatures(nfeatures = ____) |> ScaleData() |> RunPCA(verbose = FALSE)
## TODO ▶ 整合後用幾個維度？（Q3 頁 14）
gbm4 <- FindNeighbors(gbm4, dims = 1:____) |> FindClusters(resolution = ____, cluster.name = "raw_clusters")
## TODO ▶ 整合後用幾個維度？（Q3 頁 14）
gbm4 <- RunUMAP(gbm4, dims = 1:____, reduction.name = "umap.raw")
p <- DimPlot(gbm4, reduction = "umap.raw", group.by = "patient") +
     DimPlot(gbm4, reduction = "umap.raw", group.by = "celltype_author", label = TRUE, repel = TRUE)
ggsave("output/figs/04_umap_unintegrated.png", p, width = 13, height = 5.5, dpi = 150, bg = "white")
# 看圖：Neoplastic 是否按病人分成不同群集？Immune / Oligodendrocyte 是否跨病人混合？
saveRDS(gbm4, "output/rds/04_gbm4_unintegrated.rds")

# 量化「按病人分成不同群集」：每群裡病人的熵（越低越單一病人）
patient.entropy <- gbm4@meta.data |> dplyr::count(raw_clusters, patient) |>
  group_by(raw_clusters) |> mutate(p = n / sum(n)) |>
  summarise(entropy = -sum(p * log(p)), n = sum(n), .groups = "drop")
print(patient.entropy)
write.csv(patient.entropy, "output/tables/04_patient_entropy_raw.csv", row.names = FALSE)

## ---- 4. integrate ------------------------------------------------- Q3 頁 14–15
# 預設走 Seurat 自己的整合（CCA 錨點；大資料可改 RPCAIntegration）。Harmony 是「Seurat 整合不理想時」的替代方案，
# 本練習兩個都跑、放在同一個物件的兩個 reduction 裡，用同一套檢查比較。
## TODO ▶ Harmony 以什麼當批次？（Q3 頁 14）
gbm4[["RNA"]] <- split(gbm4[["RNA"]], f = gbm4$____)   # 每位病人一個 layer（批次 = 病人）
## TODO ▶ 挑幾個高變異基因？（Q2 頁 25）
gbm4 <- NormalizeData(gbm4) |> FindVariableFeatures(nfeatures = ____) |> ScaleData() |> RunPCA(verbose = FALSE)

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
  ## TODO ▶ 整合後用幾個維度？（Q3 頁 14）
  gbm4 <- FindNeighbors(gbm4, reduction = red, dims = 1:____, graph.name = paste0(tag, "_snn")) |>
          ## TODO ▶ 分群的 resolution 從掃描與 clustree 挑哪一層？（Q2 頁 33–34）
          FindClusters(graph.name = paste0(tag, "_snn"), resolution = ____, cluster.name = paste0(tag, "_clusters")) |>
          ## TODO ▶ 整合後用幾個維度？（Q3 頁 14）
          RunUMAP(reduction = red, dims = 1:____, reduction.name = paste0("umap.", tag))
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
ctrl.wide <- tidyr::pivot_wider(ctrl, names_from = celltype_author, values_from = mean_entropy)
print(ctrl.wide)
write.csv(ctrl.wide, "output/tables/04_integration_controls.csv", row.names = FALSE)
# 讀法：Immune / Oligodendrocyte 的熵 整合後 > 未整合 = 正對照通過；
#       Neoplastic 的熵若也被拉到跟免疫一樣高 = 負對照失敗（過度修正）。兩種方法都會把惡性細胞混掉一部分——這正是「兩難」。
# 決定的方式：正負對照與 LISI 一起看，兩個方法並排比。要有心理準備負對照多半不會「通過」——
#   你跑出來多半就是這樣。那時候比的是「誰失得少」，不是「誰過關」。
#
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）

# 整合原則：整合空間只用來註釋正常細胞、對齊免疫細胞；惡性細胞回未整合空間、按病人看。
# （存檔統一放在 §5 結尾，那時 LISI 欄位才算完；這裡先不存，免得同一個大物件寫兩次）

## ---- 5. integration-metrics ---------------------------------------- Q3 頁 16
# 把「混得好不好」變成數字：每群病人組成 + LISI，兩個整合空間各算一次
round(prop.table(table(gbm4$cca_clusters, gbm4$patient), 1), 2)          # 正常細胞的群：四位病人都該有
library(lisi)                                     # remotes::install_github("immunogenomics/lisi")
for (red in c("pca", "integrated.cca", "harmony")) {
  emb <- Embeddings(gbm4, red)[, 1:25]
  gbm4[[paste0("lisi_", sub("integrated.", "", red, fixed = TRUE))]] <- compute_lisi(emb, gbm4@meta.data, "patient")$patient
}
lisi.tab <- gbm4@meta.data |> group_by(celltype_author) |>
  summarise(across(starts_with("lisi_"), ~ round(median(.x), 2)))
print(lisi.tab)
write.csv(lisi.tab, "output/tables/04_lisi_by_celltype.csv", row.names = FALSE)
# 免疫／寡樹突的 LISI 接近病人數（4）= 整合成功；Neoplastic 接近 1 = 病人差異保住了（正確，不是失敗）
# 比較 lisi_cca 與 lisi_harmony：正常細胞誰混得好、惡性細胞誰保得住——這就是選整合方法的依據。
#
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）
saveRDS(gbm4, "output/rds/04_gbm4_integrated.rds")

# =====================================================================
# ▶ 練習 4
#  4-1 在 umap.raw 上用 split.by = "patient" 畫 celltype_author：哪些型別每位病人都有？哪些只在某些病人？
#  4-2 用 §4 的正負對照熵：Neoplastic 群的病人熵是多少？Immune 群呢？用一句話解釋差異的生物學原因。
#  4-3 再加 RPCAIntegration（只改 method），與 CCA、Harmony 三者的正負對照熵與 LISI。
#      哪一個對惡性細胞「修得最少」？在腫瘤資料裡這是優點還是缺點？
#  4-4 只取 Immune cell 子集，比較整合前後的分群：整合後能否分出小膠質 vs 巨噬細胞？
#      （提示：P2RY12 / TMEM119 vs CD163 / LYZ）
#  4-5 §5 的 LISI：Neoplastic 與 Immune 的中位數各是多少？把 Harmony 的 theta 調到 4 重跑，
#      Neoplastic 的 LISI 變成多少？這代表什麼？
#  進階 讀 Neftel 2019 的 Methods：他們如何處理病人效應？跟本腳本的原則一致嗎？
# =====================================================================
