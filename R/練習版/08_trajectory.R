# =====================================================================
# 08_trajectory.R — 練習腳本 8：軌跡分析（Slingshot + tradeSeq）——一位病人的惡性細胞
#
# 對應影片：Q3 頁 66–69（§1 Slingshot + tradeSeq、§2 Monocle3 與手動選起點、§3 Monocle2 選配）
# 輸入：output/rds/06_gbm4_final.rds
# 輸出：output/figs/08_*.pdf、output/tables/08_traj_association.csv
# 時間：§1 約 5–10 分鐘；§2 約 2 分鐘；§3（選配）約 5–10 分鐘
# 安裝（選配段）：見 00_setup.R——monocle3 + SeuratWrappers（GitHub）、monocle（Bioconductor）
# 前提：軌跡假設「連續過程」；跨病人混做會把病人差異當成軌跡，所以只在一位病人的惡性細胞內做。
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/rds/06_gbm4_final.rds")
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. trajectory -------------------------------------------------- Q3 頁 67–68
library(slingshot); library(tradeSeq)
# 只在「同一位病人的惡性細胞」內做；subset 之後整條流程重跑
mal1 <- subset(gbm4, malignant == "malignant" & patient == "BT_S2")
mal1 <- NormalizeData(mal1) |> FindVariableFeatures() |> ScaleData() |> RunPCA(verbose = FALSE) |>
        FindNeighbors(dims = 1:15) |> FindClusters(resolution = 0.4) |> RunUMAP(dims = 1:15)
# 起點要有理由：用 Neftel 分數（03 的 neftel 清單）找 OPC 樣分數最高的群
neftel <- list(OPC = c("BCAN", "PLP1", "GPR17", "FIBIN", "LHFPL3", "OLIG1"),
               MES = c("CHI3L1", "ANXA2", "ANXA1", "CD44", "VIM", "MT2A"))
mal1 <- AddModuleScore(mal1, features = neftel, name = names(neftel))
root <- names(which.max(tapply(mal1$OPC1, Idents(mal1), mean)))
sce  <- as.SingleCellExperiment(mal1)
## TODO ▶ 軌跡的起點要有生物學理由：哪一群？（Q3 頁 67–68）
sce  <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP", start.clus = ____)
mal1$pt <- slingPseudotime(sce)[, 1]
p <- FeaturePlot(mal1, features = "pt") + ggtitle(paste("pseudotime, root =", root))
ggsave("output/figs/08_pseudotime.pdf", p, width = 5, height = 4, bg = "white")
# 檢查：週期有沒有主導軌跡？（gbm4 這條多病人流程沒跑過 CellCycleScoring，先在這個子集算）
mal1 <- CellCycleScoring(mal1, s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
cor(mal1$pt, mal1$S.Score, use = "complete.obs"); cor(mal1$pt, mal1$G2M.Score, use = "complete.obs")   # |r| > 0.5 → 先回歸再做
# 沿 pseudotime 變化的基因（負二項 GAM）
keep   <- !is.na(mal1$pt)
counts <- LayerData(mal1, layer = "counts")[VariableFeatures(mal1), keep]
## TODO ▶ GAM 的節點數（Q3 頁 68）
gam    <- fitGAM(counts = as.matrix(counts), pseudotime = mal1$pt[keep], cellWeights = rep(1, sum(keep)), nknots = ____)
assoc  <- associationTest(gam); assoc <- assoc[order(assoc$pvalue), ]; head(assoc, 20)
write.csv(assoc, "output/tables/08_traj_association.csv")
# 畫 smoothers 的基因必須在 gam 模型裡（= 這個子集的 VariableFeatures）；
# 硬寫名字會踩「不在高變異清單」的雷——先跟 rownames(counts) 取交集，不足就用 associationTest 前幾名補
show.genes <- intersect(c("OLIG1", "SOX4", "CD44", "VIM"), rownames(counts))
show.genes <- head(unique(c(show.genes, rownames(assoc))), 4)
cat("smoothers 畫這些基因：", paste(show.genes, collapse = ", "), "\n")
pdf("output/figs/08_smoothers.pdf", 8, 6); for (g in show.genes) print(plotSmoothers(gam, as.matrix(counts), gene = g) + ggtitle(g)); dev.off()
# 穩健性：換 seed 重跑分群 + slingshot，pseudotime 的 Spearman 相關應 > 0.8；再換一位病人看方向

## ---- 2. monocle3（含「自己選起點」的示範）--------------------------- Q3 頁 69
# Monocle3 是最多人用的軌跡工具之一（graph-based、允許分支）。它的標準流程本來就要求你「選起點」，
# 有兩種姿勢：
#   (a) 手動互動式：order_cells(cds) 不給參數，RStudio 會跳出視窗讓你「點」起點——教學上最直觀，
#       但每次點的可能不一樣，不可重現，正式分析不建議單獨使用
#   (b) 腳本化：把「起點的理由」寫成程式碼（如 OPC 樣分數最高的細胞），可重現
PICK_ROOT_BY_HAND <- FALSE                     # 想體驗手動點選就改 TRUE（要在 RStudio 互動環境跑）
if (requireNamespace("monocle3", quietly = TRUE) && requireNamespace("SeuratWrappers", quietly = TRUE)) {
  library(monocle3)
  cds <- SeuratWrappers::as.cell_data_set(mal1)              # Seurat → cell_data_set（沿用 UMAP 與分群）
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)             # 在 UMAP 上學一張主圖（可分支）
  if (PICK_ROOT_BY_HAND && interactive()) {
    cds <- order_cells(cds)                                  # (a) 跳出視窗，點你認為的起點，按 Done
  } else {
    root.cells <- colnames(mal1)[order(mal1$OPC1, decreasing = TRUE)[1:10]]   # (b) OPC 樣分數前 10 顆
    cds <- order_cells(cds, root_cells = root.cells)
  }
  mal1$pt_m3 <- monocle3::pseudotime(cds)[colnames(mal1)]
  mal1$pt_m3[!is.finite(mal1$pt_m3)] <- NA                   # 圖上到不了的細胞是 Inf → 改 NA
  p <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = TRUE,
                  label_leaves = FALSE, label_roots = TRUE) + ggtitle("Monocle3 pseudotime")
  ggsave("output/figs/08_monocle3_pseudotime.pdf", p, width = 5.5, height = 4.5, bg = "white")
  # 跨工具檢查：兩套 pseudotime 的 Spearman 相關高（> 0.8）→ 結論不依賴工具，可信度加分
  cat("Slingshot vs Monocle3 pseudotime Spearman r =",
      round(cor(mal1$pt, mal1$pt_m3, method = "spearman", use = "complete.obs"), 3), "\n")
} else {
  message("未安裝 monocle3 / SeuratWrappers（見 00_setup.R 的選配段），跳過 §2。")
}

## ---- 3. monocle2（選配；經典 DDRTree，root_state 也是自己選）--------- Q3 頁 69
# Monocle2 是最早流行的版本（DDRTree），老論文常見。它把細胞分成幾個 State，
# orderCells(root_state = ...) 就是「使用者自己選起點」：先畫圖看哪個 State 該當起點，再指定。
if (requireNamespace("monocle", quietly = TRUE)) {
  library(monocle)
  # --- igraph ≥ 2.1 相容補丁（不加會在 orderCells() 報錯）---------------
  # igraph 2.1.0 把頂點索引的輔助函數 nei() 廢除（defunct，改名 .nei()），Monocle2 沒跟上，
  # orderCells() 走訪 MST 時會報「`nei()` was deprecated in igraph 2.1.0 and is now defunct」。
  # 注意：網路上常見的 assignInNamespace("nei", ...) 修不了——igraph 是在索引當下「函數內部」
  # 重新定義 nei 這個報錯版，蓋掉命名空間裡的任何東西。正解是反過來改 monocle：把它命名空間裡
  # 所有還在呼叫 nei( 的函數就地改寫成 .nei(。只影響本次 R session，不動安裝檔。
  if (utils::packageVersion("igraph") >= "2.1.0") {
    ns <- getNamespace("monocle")
    for (fn in ls(ns)) {
      f <- get(fn, envir = ns)
      if (!is.function(f)) next
      b <- deparse(body(f))
      if (!any(grepl("(?<![.[:alnum:]_])nei\\(", b, perl = TRUE))) next
      b <- gsub("(?<![.[:alnum:]_])nei\\(", ".nei(", b, perl = TRUE)
      body(f) <- parse(text = paste(b, collapse = "\n"))[[1]]
      assignInNamespace(fn, f, ns = "monocle")
      cat("  igraph 相容補丁：已改寫 monocle:::", fn, "\n", sep = "")
    }
  }
  cnt2 <- as.matrix(LayerData(mal1, layer = "counts")[VariableFeatures(mal1), ])
  cds2 <- newCellDataSet(cnt2,
                         phenoData   = new("AnnotatedDataFrame", data = mal1@meta.data),
                         featureData = new("AnnotatedDataFrame",
                                           data = data.frame(gene_short_name = rownames(cnt2),
                                                             row.names = rownames(cnt2))),
                         expressionFamily = negbinomial.size())
  cds2 <- estimateSizeFactors(cds2); cds2 <- estimateDispersions(cds2)
  cds2 <- reduceDimension(cds2, max_components = 2, method = "DDRTree")
  cds2 <- orderCells(cds2)                                   # 第一次先不指定，讓它自己排
  # 「自己選起點」：先看這張圖，決定哪個 State 是起點——
  p <- plot_cell_trajectory(cds2, color_by = "State") + ggtitle("Monocle2: pick your root State")
  ggsave("output/figs/08_monocle2_states.pdf", p, width = 5.5, height = 4.5, bg = "white")
  # 手動版：看圖後把數字填進去，例如 cds2 <- orderCells(cds2, root_state = 3)
  # 腳本版（可重現）：選 OPC 樣分數最高的 State
  st.score <- tapply(mal1$OPC1[colnames(cds2)], pData(cds2)$State, mean)
  cds2 <- orderCells(cds2, root_state = as.integer(names(which.max(st.score))))
  mal1$pt_m2 <- pData(cds2)$Pseudotime[match(colnames(mal1), colnames(cds2))]
  p <- plot_cell_trajectory(cds2, color_by = "Pseudotime") + ggtitle("Monocle2 pseudotime (DDRTree)")
  ggsave("output/figs/08_monocle2_pseudotime.pdf", p, width = 5.5, height = 4.5, bg = "white")
  cat("Slingshot vs Monocle2 pseudotime Spearman r =",
      round(cor(mal1$pt, mal1$pt_m2, method = "spearman", use = "complete.obs"), 3), "\n")
} else {
  message("未安裝 monocle（Monocle2，見 00_setup.R 的選配段），跳過 §3。")
}
# 三套工具的共同結論才寫進論文：起點端（OPC 樣）→ 終點端（MES 樣）的方向、以及沿路在動的基因。

sessionInfo()

# =====================================================================
# ▶ 練習 8
#  8-1 把 start.clus 改成 MES 樣分數最高的群，pseudotime 反過來了嗎？哪些基因的 associationTest 結果不變？這說明什麼？
#  8-2 換另一位病人重跑，OLIG1 → CD44 的方向一致嗎？四位病人各畫一張 smoothers 併排。
#  8-3 cor(pt, S.Score) 與 cor(pt, G2M.Score) 各是多少？若 |r| > 0.5，回歸掉週期之後軌跡還在嗎？
#  8-4 把 PICK_ROOT_BY_HAND 改成 TRUE，在 Monocle3 的視窗裡故意點「MES 樣那端」當起點：
#      pseudotime 反轉了嗎？跟 8-1 的結論合起來，寫一句「起點選擇影響什麼、不影響什麼」。
#  8-5 三套工具（Slingshot / Monocle3 / Monocle2）的 pseudotime 兩兩 Spearman 相關各是多少？
#      分支結構一致嗎？不一致時你會相信誰、為什麼？
#  進階 手動點選（order_cells 互動視窗）跟腳本化選起點各適合什麼場景？
#      為什麼正式分析建議「探索用手動、定稿用腳本」？
# =====================================================================
