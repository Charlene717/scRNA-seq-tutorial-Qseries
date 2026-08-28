# =====================================================================
# 07_cellchat.R — 練習腳本 7：細胞通訊（CellChat）——每個樣本各跑一次、六種圖、兩條件比較、LIANA 交叉驗證
#
# 對應影片：Q3 頁 54–65（§1 跑一次 CellChat、§2 路徑層級與六種圖、§3 兩條件比較、§4 LIANA）
# 輸入：output/gbm4_final.rds（06a_pseudobulk_gsea.R；含 malignant 標籤與 type 欄）
# 輸出：output/cellchat/<patient>_<tissue>.rds、output/figs/07_*.pdf
# 時間：每個樣本約 5–15 分鐘（8 個樣本，建議先跑一位病人）
# 安裝：devtools::install_github("jinworks/CellChat")；LIANA：remotes::install_github("saezlab/liana")
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(CellChat); library(patchwork)
set.seed(1234)
gbm4 <- readRDS("output/gbm4_final.rds")
gbm4[["RNA"]] <- JoinLayers(gbm4[["RNA"]])                # 04 之後 RNA 是按病人分層的；合併回單一 data 層，
if (!"data" %in% Layers(gbm4[["RNA"]])) gbm4 <- NormalizeData(gbm4)   # 否則 plotGeneExpression/VlnPlot 會拿 counts 畫
dir.create("output/cellchat", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figs", showWarnings = FALSE, recursive = TRUE)

# 通訊分析用的標籤：CNV 判定後的惡性 + 作者的正常型別；unresolved 排除
gbm4$cc_label <- ifelse(gbm4$malignant == "malignant", "Malignant",
                 ifelse(gbm4$celltype_author == "Immune cell", "Macro/MG",   # 本資料免疫細胞以髓系為主
                        gbm4$celltype_author))
gbm4 <- subset(gbm4, malignant != "unresolved")
table(gbm4$cc_label, paste(gbm4$patient, gbm4$tissue))                      # 每群 ≥ 20 顆才可靠

## ---- 1. run-per-sample --------------------------------------------- Q3 頁 55–56
run_cc <- function(obj, label = "cc_label") {
  obj$samples <- factor(paste(obj$patient, obj$tissue, sep = "_"))          # CellChat v2 要求 meta 有 samples 欄
  cc <- createCellChat(object = obj, group.by = label, assay = "RNA")       # 用 data 層（log-normalized）
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")        # 先只看分泌型；熟了可用全部
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE)     # 注意大小寫 triMean；群大小校正
  cc <- filterCommunication(cc, min.cells = 20)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cc
}
samples <- unique(paste(gbm4$patient, gbm4$tissue, sep = "_"))
cc.all <- list()
for (s in samples) {                                                        # 8 個樣本；先跑一位病人也可以
  p <- sub("_[^_]*$", "", s); t <- sub(".*_", "", s)      # 病人 ID 含底線，從最後一段切
  obj <- subset(gbm4, patient == p & tissue == t)
  if (ncol(obj) < 200) next
  cc.all[[s]] <- run_cc(obj)
  saveRDS(cc.all[[s]], paste0("output/cellchat/", s, ".rds"))
}
cc <- cc.all[["BT_S2_Tumor"]]                                               # 以下用一個樣本示範六種圖

## ---- 2. six-plots --------------------------------------------------- Q3 頁 57–62
groupSize <- as.numeric(table(cc@idents))
pdf("output/figs/07_1_circle.pdf", width = 10, height = 5); par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cc@net$count,  vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cc@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights")
dev.off()
# 圖 2：熱圖（誰送給誰）
p2 <- netVisual_heatmap(cc, measure = "weight", color.heatmap = "Reds"); pdf("output/figs/07_2_heatmap.pdf", 6, 5); print(p2); dev.off()
# 圖 4：signaling role（先做，找主角）
cc@netP$pathways
p4a <- netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing", height = 8)
p4b <- netAnalysis_signalingRole_heatmap(cc, pattern = "incoming", height = 8)
pdf("output/figs/07_4_roles.pdf", 12, 7); print(p4a + p4b); dev.off()
netAnalysis_signalingRole_scatter(cc); ggsave("output/figs/07_4_scatter.pdf", width = 5, height = 4, bg = "white")
# 圖 3：挑一條路徑拆開（以 SPP1 為例；名稱以 cc@netP$pathways 為準）
pw <- "SPP1"
pdf("output/figs/07_3_pathway.pdf", 10, 5); par(mfrow = c(1, 2))
netVisual_aggregate(cc, signaling = pw, layout = "circle")
netVisual_aggregate(cc, signaling = pw, layout = "chord")
dev.off()
netAnalysis_contribution(cc, signaling = pw); ggsave("output/figs/07_3_contribution.pdf", width = 5, height = 3, bg = "white")
netVisual_heatmap(cc, signaling = pw, color.heatmap = "Reds")
plotGeneExpression(cc, signaling = pw); ggsave("output/figs/07_3_genes.pdf", width = 8, height = 5, bg = "white")   # 驗證表現
netAnalysis_signalingRole_network(cc, signaling = pw, width = 8, height = 2.5)
# 圖 5：bubble（寫進論文的那張）
netVisual_bubble(cc, sources.use = c("Macro/MG", "Malignant"),
                 targets.use = c("Malignant", "Macro/MG", "Vascular"), remove.isolate = TRUE)
ggsave("output/figs/07_5_bubble.pdf", width = 7, height = 8, bg = "white")
# 只看幾條路徑：名稱一定要在 cc@netP$pathways 裡（每個樣本推得出的路徑不同，硬寫 MIF 這種名稱會報錯）
sig3 <- head(intersect(c("SPP1", "MIF", "VEGF", "PTN", "TGFb"), cc@netP$pathways), 3)
if (length(sig3) == 0) sig3 <- head(cc@netP$pathways, 3)
netVisual_bubble(cc, signaling = sig3, remove.isolate = TRUE)
# 回 Seurat 驗證具體的一對
VlnPlot(subset(gbm4, patient == "BT_S2" & tissue == "Tumor"), features = c("SPP1", "CD44"), group.by = "cc_label", pt.size = 0)

## ---- 3. compare-conditions ------------------------------------------ Q3 頁 63–64
cc.list <- list(Core = cc.all[["BT_S2_Tumor"]], Periphery = cc.all[["BT_S2_Periphery"]])
cc.m <- mergeCellChat(cc.list, add.names = names(cc.list))
g1 <- compareInteractions(cc.m, show.legend = FALSE, group = c(1, 2))
g2 <- compareInteractions(cc.m, show.legend = FALSE, group = c(1, 2), measure = "weight")
g1 + g2; ggsave("output/figs/07_6_compare.pdf", width = 6, height = 3, bg = "white")
pdf("output/figs/07_6_diff.pdf", 10, 5); par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cc.m, weight.scale = TRUE)                       # 紅：Periphery > Core；藍：反之
netVisual_diffInteraction(cc.m, weight.scale = TRUE, measure = "weight")
dev.off()
rankNet(cc.m, mode = "comparison", stacked = TRUE, do.stat = TRUE); ggsave("output/figs/07_6_rankNet.pdf", width = 5, height = 6, bg = "white")
# 兩條件並排的 bubble：這對「來源 → 目標」必須在兩個條件裡都有顯著互動才畫得出來；
# 只有一邊有的時候，CellChat 會丟出 seq 的 'by' 錯誤（已知 bug）——用 tryCatch 接住，退回單樣本各畫一張對照
p <- tryCatch(
  netVisual_bubble(cc.m, sources.use = "Macro/MG", targets.use = "Malignant", comparison = c(1, 2), angle.x = 45),
  error = function(e) {
    message("  Macro/MG → Malignant 在其中一個條件沒有顯著互動，改畫兩張單樣本 bubble 對照（", conditionMessage(e), "）")
    (netVisual_bubble(cc.list$Core,      sources.use = "Macro/MG", targets.use = "Malignant", remove.isolate = TRUE) + ggtitle("Core")) +
    (netVisual_bubble(cc.list$Periphery, sources.use = "Macro/MG", targets.use = "Malignant", remove.isolate = TRUE) + ggtitle("Periphery"))
  })
print(p); ggsave("output/figs/07_6_bubble_compare.pdf", p, width = 9, height = 6, bg = "white")
# 四位病人一致性：對每位病人重複 §3，收集 rankNet 的顯著路徑，取交集
# rank.list <- lapply(patients, function(p) rankNet(mergeCellChat(list(cc.all[[paste0(p,"_Tumor")]], cc.all[[paste0(p,"_Periphery")]]), add.names = c("Core","Periphery")), mode = "comparison", do.stat = TRUE, return.data = TRUE)$signaling.contribution)

## ---- 4. liana-crosscheck --------------------------------------------- Q3 頁 65
if (requireNamespace("liana", quietly = TRUE)) {
  library(liana)
  obj <- subset(gbm4, patient == "BT_S2" & tissue == "Tumor"); Idents(obj) <- "cc_label"
  li <- liana_wrap(obj, method = c("natmi", "connectome", "sca", "cellphonedb"), resource = "Consensus")
  li.agg <- liana_aggregate(li)                                              # 共識排名（aggregate_rank 越小越好）
  # dplyr:: 寫全名：org.Hs.eg.db／AnnotationDbi（06a 的 ORA 會載入）也有 select、filter，
  # 同一個 R session 先跑過 06a 再跑這裡會被蓋掉
  print(li.agg |> dplyr::filter(source == "Macro/MG", target == "Malignant") |>
        dplyr::select(source, target, ligand.complex, receptor.complex, aggregate_rank) |> head(10))
  p <- li.agg |> liana_dotplot(source_groups = "Macro/MG", target_groups = c("Malignant", "Vascular"), ntop = 15)
  print(p); ggsave("output/figs/07_liana_dotplot.pdf", p, width = 10, height = 6, bg = "white")
  write.csv(li.agg |> dplyr::select(-starts_with("natmi"), -starts_with("connectome")) |> head(500),
            "output/07_liana_top500.csv", row.names = FALSE)
  # 與 CellChat 的 bubble 對照：SPP1–CD44、MIF–CD74 兩邊都在前段，才寫進結果
}
sessionInfo()

# =====================================================================
# ▶ 練習 7
#  7-1 把 subsetDB 改成全部三類（不加 search），互動數變多少？新增的顯著路徑主要是哪一類？
#  7-2 圖 4：這個樣本裡「既送又收」的樞紐是誰？換另一位病人，樞紐一樣嗎？
#  7-3 SPP1 路徑的 netAnalysis_contribution 前兩對是什麼？用 VlnPlot 確認：配體在來源群的表現比例 > 25% 嗎？
#  7-4 四位病人各做一次 rankNet 比較，哪些路徑四位方向一致？只有一位病人顯著的路徑有幾條？
#  7-5 LIANA 的共識前 10 對與 CellChat bubble 的前 10 對重疊幾對？不重疊的原因可能是什麼？
#  進階 用 NicheNet 反推：邊緣惡性細胞相對核心上調的基因（06 的 DE），最能被哪個配體解釋？與 CellChat 的結論一致嗎？
# =====================================================================
