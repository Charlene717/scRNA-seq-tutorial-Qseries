# =====================================================================
# 08_trajectory.R — 練習腳本 8：軌跡分析（Slingshot + tradeSeq）——一位病人的惡性細胞
#
# 對應影片：Q3 頁 72–75（§1 Slingshot + tradeSeq、§2 Monocle3 與手動選起點、§3 Monocle2 選配）
# 輸入：output/rds/06_gbm4_final.rds
# 輸出：output/figs/08_*.pdf、output/tables/08_traj_association.csv
#       （有分支時另加 08_pseudotime_lineages.pdf、08_traj_association_bylineage.csv、08_traj_diffend.csv）
# 時間：§1 約 5–10 分鐘（每多一條 lineage，fitGAM 約等比例增加）；§2 約 2 分鐘；§3（選配）約 5–10 分鐘
# 安裝（選配段）：見 00_setup.R——monocle3 + SeuratWrappers（GitHub）、monocle（Bioconductor）
# 前提：軌跡假設「連續過程」；跨病人混做會把病人差異當成軌跡，所以只在一位病人的惡性細胞內做。
# ⚠ 這份資料在 BT_S2 上會跑出分支（不只一條 lineage）。分支不是錯誤：§1 會照實際條數配權重，
#   §1c 再比較分支之間有沒有差。只有需要「單一條」的圖與早晚比較才固定用 lineage 1，並在圖名標示。
# ⚠ 軌跡對 doublet 也特別敏感：doublet 落在兩群中間，正好會被連成一條假的過渡路徑或假分支，
#   是最容易生出假故事的一種分析。本節的 GSE84465 是 Smart-seq2，沒有跑 doublet 偵測（見 02 §4c），
#   所以這裡不必找 doublet_status；換成 10x 資料做軌跡時，務必先照 02 §4c 處理過再來。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/rds/06_gbm4_final.rds")
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. trajectory -------------------------------------------------- Q3 頁 73–74
library(slingshot); library(tradeSeq)
# 只在「同一位病人的惡性細胞」內做；subset 之後整條流程重跑
mal1 <- subset(gbm4, malignant == "malignant" & patient == "BT_S2")
mal1 <- NormalizeData(mal1) |> FindVariableFeatures() |> ScaleData() |> RunPCA(verbose = FALSE) |>
        FindNeighbors(dims = 1:15) |> FindClusters(resolution = 0.4) |> RunUMAP(dims = 1:15)
# 起點要有理由：用 Neftel 分數找 OPC 樣分數最高的群。
# ⚠ 下面每組只列 6 個基因，是為了讓這一步跑得快的「示範用簡表」，不是 Neftel 的完整 signature
#   （原文 Table S2 每組約 50 個；03 用的是前 12 個）。它只用來挑起點，不要拿這個分數報狀態比例。
neftel <- list(OPC = c("BCAN", "PLP1", "GPR17", "FIBIN", "LHFPL3", "OLIG1"),
               MES = c("CHI3L1", "ANXA2", "ANXA1", "CD44", "VIM", "MT2A"))
mal1 <- AddModuleScore(mal1, features = neftel, name = names(neftel))
root <- names(which.max(tapply(mal1$OPC1, Idents(mal1), mean)))
sce  <- as.SingleCellExperiment(mal1)
sce  <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP", start.clus = root)
# Slingshot 可能找出不只一條 lineage。分支不是錯誤，是這份資料本身的結構——
# 真正危險的是「有分支卻沒發現」：直接取 [, 1] 會默默只分析第一條，其他分支專屬的細胞
# 被 NA 濾掉也沒人知道，報出來的「這條軌跡」其實只是全貌的一半。
# 所以這裡先把條數印出來，再照實際條數配權重；tradeSeq 本來就支援多條 lineage。
pt.all <- slingPseudotime(sce)                 # 含 NA 版：不屬於該 lineage 的細胞留空，描述性統計用這版
pt.mat <- slingPseudotime(sce, na = FALSE)     # 不含 NA 版：fitGAM 用這版
cw     <- slingCurveWeights(sce)               # 每顆細胞屬於各 lineage 的權重
n.lin  <- ncol(pt.mat)
cat("\n== Slingshot 找出", n.lin, "條 lineage ==\n"); print(slingLineages(sce))
if (n.lin > 1)
  cat("有分支：這是一個結果，不是錯誤。接下來這樣處理——\n",
      "  · fitGAM 吃整個 pseudotime 矩陣與 curve weights，每條分支各配一組曲線\n",
      "  · associationTest 預設是「跨所有 lineage 的整體檢定」；要分支各自的結果加 lineages = TRUE\n",
      "  · 分支之間的比較用 diffEndTest()（終點差異），見 §1c\n",
      "  · FeaturePlot 一次只能用一條上色，本節固定用 lineage 1，圖名也標清楚\n", sep = "")
mal1$pt <- pt.all[, 1]                         # 顯示與早晚比較都用 lineage 1（其他分支專屬的細胞為 NA）
p <- FeaturePlot(mal1, features = "pt") +
     ggtitle(sprintf("pseudotime | lineage 1 of %d | root = %s", n.lin, root))
ggsave("output/figs/08_pseudotime.pdf", p, width = 5, height = 4, bg = "white")
if (n.lin > 1) {                               # 有分支就每條都畫，才看得出它們分別走去哪
  pl <- lapply(seq_len(n.lin), function(i) {
    m <- mal1; m$pt_i <- pt.all[, i]
    FeaturePlot(m, features = "pt_i") + ggtitle(paste("lineage", i)) })
  ggsave("output/figs/08_pseudotime_lineages.pdf", patchwork::wrap_plots(pl, ncol = 2),
         width = 9, height = 4 * ceiling(n.lin / 2), bg = "white")
}
# 檢查：週期有沒有主導軌跡？（gbm4 這條多病人流程沒跑過 CellCycleScoring，先在這個子集算）
mal1 <- CellCycleScoring(mal1, s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
for (i in seq_len(n.lin))                      # 每條都要看：有可能只有其中一條被細胞週期帶著走
  cat(sprintf("lineage %d：cor(pt, S.Score) = %+.3f｜cor(pt, G2M.Score) = %+.3f\n", i,
              cor(pt.all[, i], mal1$S.Score,   use = "complete.obs"),
              cor(pt.all[, i], mal1$G2M.Score, use = "complete.obs")))   # |r| > 0.5 → 先回歸再做
# 沿 pseudotime 變化的基因（負二項 GAM）
# 權重矩陣決定每顆細胞對哪條曲線有貢獻：單一軌跡時它就是一整欄 1，有分支時各欄不同，
# 所以下面這行不必分兩種寫法，n.lin = 1 時也成立。
keep    <- rowSums(cw) > 0                     # 沒被任何 lineage 收下的細胞（極少）先排除
cnt.fit <- as.matrix(LayerData(mal1, layer = "counts")[VariableFeatures(mal1), keep])
# ⚠ 每多一條 lineage，fitGAM 就多配一組平滑曲線，時間大致等比例增加。
gam     <- fitGAM(counts = cnt.fit, pseudotime = pt.mat[keep, , drop = FALSE],
                  cellWeights = cw[keep, , drop = FALSE], nknots = 6)
# 排序要用 waldStat，不能用 pvalue：單一軌跡那次跑，有 38 個基因的 p 直接下溢成 0，
# order(pvalue) 在它們之間是任意順序——TAGLN（waldStat 91）會排到第一，
# 而真正最強的 GFAP（849）掉到第六，「前幾名」就變成假的。
# （這組數字來自改成分支版之前的執行，改用聯合配適後名次會變動，待重跑更新；
#   但「排序用 waldStat 不用 pvalue」這個結論不會因此改變。）
assoc  <- associationTest(gam); assoc <- assoc[order(-assoc$waldStat), ]; head(assoc, 20)
# 兩欄不要看混，它們排出來的名次常常不一樣：
#   waldStat  = 變化模式有多「明確」（效應量 ÷ 不確定性）——跟第 38 頁講 DESeq2 的 stat 同一個道理
#   meanLogFC = 變化幅度有多「大」，而且是絕對值，看不出是升還是降
# 舊版那次：GPR37L1 幅度 4.94 卻只排第 20（waldStat 132）；COL1A2 幅度僅 0.74 卻排第 9（231）。
# （同樣是分支版之前的數字，待重跑更新。）挑基因畫圖、決定「誰在動」用 waldStat；報告效應量用 meanLogFC。
write.csv(assoc, "output/tables/08_traj_association.csv")
if (n.lin > 1) {   # 整體檢定只說「這個基因在某處有變化」，沒說是哪一條分支；要分支各自的結果就加這個
  assoc.lin <- associationTest(gam, lineages = TRUE)
  write.csv(assoc.lin[order(-assoc.lin$waldStat), ], "output/tables/08_traj_association_bylineage.csv")
}
# 畫 smoothers 的基因必須在 gam 模型裡（= 這個子集的 VariableFeatures）。那要畫哪四個？直接取關聯最強的前四名，不要自己先想好名字再去湊——
# 原本寫死 c("OLIG1","SOX4","CD44","VIM")，結果 SOX4／CD44／VIM 根本不在這個子集的
# 2,000 個高變異基因裡，四格有三格是被 fallback 補進來的，投影片也就對不上。
show.genes <- head(rownames(assoc), 4)
cat("smoothers 畫這些基因：", paste(show.genes, collapse = ", "), "\n")
pdf("output/figs/08_smoothers.pdf", 8, 6)
for (g in show.genes) print(plotSmoothers(gam, cnt.fit, gene = g) + ggtitle(g))   # 有分支時每條各畫一條
dev.off()

## ---- 1c. 分支之間到底有沒有差？（只有多條 lineage 時才跑）------------ Q3 頁 74
# 有分支的時候，最值得問的不是「誰沿著軌跡在變」，而是「兩條路走到的地方一不一樣」。
# diffEndTest 比的是各 lineage 終點的表現量。如果幾乎沒有基因有差，那所謂的「分支」
# 比較可能是同一個狀態被曲線拆成兩半，而不是兩種命運——這種時候不該寫成「分化出兩群」。
if (n.lin > 1) {
  de.end <- diffEndTest(gam); de.end <- de.end[order(-de.end$waldStat), ]
  write.csv(de.end, "output/tables/08_traj_diffend.csv")
  cat("\n== 分支終點差異最大的 10 個基因 ==\n"); print(head(round(de.end, 3), 10))
  cat("終點差異 FDR < 0.05 的基因數：",
      sum(p.adjust(de.end$pvalue, "BH") < 0.05, na.rm = TRUE), "／", nrow(de.end), "\n")
}

# 每個基因是「往上走」還是「往下走」？associationTest 的 meanLogFC 是絕對值，
# 看不出方向，得自己比 pseudotime 兩端。但這裡有個陷阱：
# 一定要用 normalised 的 data 層，不能用 raw counts。Smart-seq2 每顆細胞的深度差好幾倍，
# 只要深度沿著 pseudotime 遞減，raw counts 會讓「每一個基因」都看起來在下降——
# 那不是生物學，是定序深度。（fitGAM 本身有 offset 校正深度，出問題的只有這種手動比較。）
#
# 這個檢查在本例的答案是「沒問題」，而且要把數字記下來：
# cor(pt, nCount_RNA) = -0.022（lineage 1），深度完全沒有沿軌跡走，
# 所以方向的判讀是真的表現量變化，不是技術假象。
# 一個回答「沒事」的檢查不是白跑的——它是你敢下結論的依據。
# 以下的早／晚比較全部只看 lineage 1，講的時候要說清楚是哪一條。
k1  <- which(!is.na(mal1$pt))
ptk <- mal1$pt[k1]; qq <- quantile(ptk, c(0.25, 0.75))
cat("深度沿軌跡的相關性 cor(pt, nCount_RNA) =",
    round(cor(mal1$pt, mal1$nCount_RNA, use = "complete.obs"), 3), "\n")   # 明顯偏離 0 就要小心
dat <- LayerData(mal1, layer = "data")[show.genes, k1]        # log1p(CP10K)，已除掉深度
trend <- t(sapply(show.genes, function(g) {
  x <- as.numeric(dat[g, ]); c(early = mean(x[ptk <= qq[1]]), late = mean(x[ptk >= qq[2]])) }))
cat("\n== 前四名沿 lineage 1 的方向 ==\n")
print(cbind(round(trend, 2),
            trend = ifelse(trend[, "late"] > trend[, "early"], "rises", "falls")))

# 終點端到底是什麼狀態？——這一步比上面那張表更重要，卻最常被略過。
# 我們是用 OPC 樣分數挑的「起點」，但從來沒有驗證過「終點」是不是我們以為的那個狀態。
# 拿 Neftel 分數比軌跡兩端——但兩半的證據力不一樣，要分開讀：
#   · OPC 從起點往終點下降：**這一半有循環成分**。起點本來就是「OPC 分數最高的群」選出來的，
#     所以 OPC 在起點高、往後降，有一部分是建構上就保證的，不能當成獨立驗證。
#   · MES 從起點往終點上升：**這一半才是真的資訊**。MES 完全沒有參與挑起點，
#     它跟著 pseudotime 上升，是這條軌跡確實對應到一個狀態轉變的證據。
# 所以下結論時靠的是 MES 那一半；OPC 那一半只能當「內部一致性檢查」（有沒有自相矛盾），
# 不能當「驗證」。想要更乾淨的做法：用一組沒參與挑起點的獨立基因（held-out）來檢查兩端。
# 有分支時還多一個問題：每條分支的終點狀態一樣嗎？所以下面逐條印，不要只看第一條。
cat("\n== Neftel 分數：每條 lineage 的起點端 vs 終點端 ==\n")
for (i in seq_len(n.lin)) {
  p.i <- pt.all[, i]; ki <- which(!is.na(p.i)); qi <- quantile(p.i[ki], c(0.25, 0.75))
  e <- ki[p.i[ki] <= qi[1]]; l <- ki[p.i[ki] >= qi[2]]
  cat(sprintf("lineage %d（n = %d）：OPC %+.3f → %+.3f｜MES %+.3f → %+.3f\n",
              i, length(ki), mean(mal1$OPC1[e]), mean(mal1$OPC1[l]),
              mean(mal1$MES2[e]), mean(mal1$MES2[l])))
}

# 前四名清一色往同一個方向的時候，也該看一眼「到底有沒有東西在反方向走」——沒有的話，
# 這條軌跡的主軸就不是「A 變成 B」，而只是「某些東西一路消失」，寫法要跟著改。
sig  <- head(rownames(assoc), 200)
d200 <- LayerData(mal1, layer = "data")[sig, k1]
dd   <- data.frame(early = rowMeans(d200[, ptk <= qq[1]]), late = rowMeans(d200[, ptk >= qq[2]]),
                   wald = assoc[sig, "waldStat"])
dd$delta <- dd$late - dd$early
cat("\n== 沿 lineage 1 最會「升」的 10 個 ==\n"); print(round(head(dd[order(-dd$delta), ], 10), 2))
cat("\n== 沿 lineage 1 最會「降」的 10 個 ==\n"); print(round(head(dd[order(dd$delta), ], 10), 2))

# 穩健性：換 seed 重跑分群 + slingshot，pseudotime 的 Spearman 相關應 > 0.8；再換一位病人看方向
# 有分支時還要多看一項：分支的「條數」和「從哪一群分岔」換 seed 後穩不穩。
# pseudotime 相關很高、但分支結構每次都不一樣，代表分支是分群的雜訊，不是生物學。
# （注意這裡的 0.8 和 §2/§3「跨工具」的 0.8 不是同一件事：同一套工具換 seed 本來就該很接近，
#   不同工具的圖形假設不同，標準要放寬——見 §3 結尾的分級。）

## ---- 2. monocle3（含「自己選起點」的示範）--------------------------- Q3 頁 75
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
  # 跨工具檢查：兩套 pseudotime 的 Spearman 相關（判讀標準見 §3 結尾的分級）
  # 注意比的是 Slingshot 的 lineage 1；Monocle3 若也分支，這個單一數字只涵蓋主幹那一段。
  cat("Slingshot vs Monocle3 pseudotime Spearman r =",
      round(cor(mal1$pt, mal1$pt_m3, method = "spearman", use = "complete.obs"), 3), "\n")
} else {
  message("未安裝 monocle3 / SeuratWrappers（見 00_setup.R 的選配段），跳過 §2。")
}

## ---- 3. monocle2（選配；經典 DDRTree，root_state 也是自己選）--------- Q3 頁 75
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
# ---- 跨工具一致性怎麼判讀：不要用單一門檻，用分級 -------------------- Q3 頁 75
# 「Spearman r > 0.8 才可信」是流傳很廣的說法，但它是慣例不是定律，而且對「不同工具」太嚴格：
# Slingshot 學的是一條主曲線、Monocle2 是 DDRTree 的樹、Monocle3 是 UMAP 上的圖，
# 三者的幾何假設不一樣，即使講的是同一件生物學，數值也不會貼得那麼近。
#   r > 0.8      強一致 → 可以直接寫「結論不依賴工具」
#   r 0.6–0.8    方向一致但細節有差 → 要有第三個「不是 pseudotime」的獨立佐證才寫結論
#   r < 0.6      不一致 → 先回頭查起點與細胞子集，這時不該報軌跡
# ⚠ 這裡引用的數字分兩類，改成分支版之後受影響的程度不一樣，不要一律當成過期：
#   不受影響：r = 0.757、Neftel OPC 0.437 → -0.383、MES 0.764 → 1.554。
#     它們只跟 lineage 1 的 pseudotime 有關，而 slingshot 的輸出沒有變——重跑會拿到同樣的值。
#   待重跑更新：waldStat 的排名與「前四名基因」（GFAP、BCAN、NDRG1、VEGFA）。
#     聯合配適之後每個基因的 waldStat 會變，名次也可能換人，請以你自己的輸出為準。
# 當時的讀法是：r = 0.757 落在中間帶，而第三個佐證有兩個，都在 §1 印過——
#   終點端的狀態分數（Neftel 兩端的 OPC／MES）、以及前幾名基因的方向。
# 這兩個都不是 pseudotime 的數值，所以「OPC 樣 → MES 樣」這個結論站得住。
# 反過來說：如果只有中間帶的 r 而沒有這兩個佐證，該寫的是「趨勢一致，待驗證」。
# 有分支的時候還要再加一句：這個結論講的是哪一條 lineage？§1c 的 diffEndTest 若顯示
#   兩條分支的終點差異很小，那就不該把它們寫成兩種命運。

sessionInfo()

# =====================================================================
# ▶ 練習 8
#  8-1 把 start.clus 改成 MES 樣分數最高的群，pseudotime 反過來了嗎？哪些基因的 associationTest 結果不變？這說明什麼？
#  8-2 換另一位病人重跑，OLIG1 → CD44 的方向一致嗎？四位病人各畫一張 smoothers 併排。
#  8-3 cor(pt, S.Score) 與 cor(pt, G2M.Score) 各是多少？若 |r| > 0.5，回歸掉週期之後軌跡還在嗎？
#  8-4 把 PICK_ROOT_BY_HAND 改成 TRUE，在 Monocle3 的視窗裡故意點「MES 樣那端」當起點：
#      pseudotime 反轉了嗎？跟 8-1 的結論合起來，寫一句「起點選擇影響什麼、不影響什麼」。
#  8-5 兩兩 Spearman 相關各是多少？照 §3 結尾的分級落在哪一格？落在中間帶時，你手上的
#      第三個獨立佐證是什麼（提示：不能也是 pseudotime）？
#      分支結構一致嗎？不一致時你會相信誰、為什麼？
#  8-6 這份資料跑出幾條 lineage？看 §1c 的 diffEndTest：終點差異達 FDR < 0.05 的基因多不多？
#      如果幾乎沒有，你會把這個結構寫成「分化出兩種命運」還是「同一個狀態被拆成兩半」？
#  8-7 把 nknots 從 6 改成 4 或 10，lineage 的條數會變嗎？associationTest 的前十名呢？
#      哪一個對參數比較敏感——說明為什麼「先確定結構、再談基因」比較安全。
#  進階 手動點選（order_cells 互動視窗）跟腳本化選起點各適合什麼場景？
#      為什麼正式分析建議「探索用手動、定稿用腳本」？
# =====================================================================
