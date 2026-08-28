# =====================================================================
# 03_annotate.R — 練習腳本 3：門牌基因、marker、SingleR、掛名字、惡性狀態分數
#
# 對應影片：Q2 頁 38–59（§1 門牌與面板、§2 FindAllMarkers 與篩選、§3 SingleR、§4 掛名字、§4b 免疫亞群、§5 Neftel 分數與三種算法、§5b 品質檢查、§6 交付與存檔）
# 輸入：output/gbm_clustered.rds（02_cluster.R）
# 輸出：output/gbm_annotated.rds、output/03_markers.csv
# 時間：約 5–10 分鐘（SingleR 首次下載參考集約 1 GB，之後有快取）
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm <- readRDS("output/gbm_clustered.rds")
Idents(gbm) <- "seurat_clusters"

## ---- 1. gates-and-panel -------------------------------------------- Q2 頁 38–40
# 四個門牌基因：先劃四塊大陸
gates <- c("PTPRC",   # 免疫
           "SOX2",    # 膠質／惡性（正常星狀、OPC 也會亮）
           "MBP",     # 寡樹突
           "PECAM1")  # 血管內皮
p <- FeaturePlot(gbm, features = gates, ncol = 4, order = TRUE)
ggsave("output/figs/03_gates.png", p, width = 16, height = 4, dpi = 150, bg = "white")

# GBM 微環境 marker 面板
panel <- list(
  immune_gate  = "PTPRC",
  macrophage   = c("CD68", "CD163", "C1QB", "AIF1", "LYZ"),
  microglia    = c("P2RY12", "TMEM119", "CX3CR1"),
  T_cell       = c("CD3E", "CD3D", "CD2", "CD8A", "IL7R"),
  glial_malig  = c("SOX2", "OLIG2", "PTPRZ1", "EGFR", "GFAP"),
  oligo        = c("MBP", "PLP1", "MOG"),
  endothelial  = c("PECAM1", "VWF", "CLDN5"),
  pericyte     = c("PDGFRB", "RGS5", "ACTA2"),
  cycling      = c("MKI67", "TOP2A"))
p <- DotPlot(gbm, features = panel, cluster.idents = TRUE) + RotatedAxis() +
     theme(axis.text.x = element_text(size = 8))
ggsave("output/figs/03_dotplot_panel.png", p, width = 16, height = 6, dpi = 150, bg = "white")
# 看圖：沿對角線一塊塊亮起來嗎？PTPRC 亮在幾群（那是免疫大陸）？MKI67 疊在哪幾群上？

## ---- 2. markers ---------------------------------------------------- Q2 頁 41–42
## TODO ▶ FindAllMarkers 的三個門檻各是什麼意思？（Q2 頁 41–42）
markers <- FindAllMarkers(gbm, only.pos = ____, min.pct = ____, logfc.threshold = ____)
top5 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 5) |>
        select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
print(top5, n = 60)                                   # 看表順序：log2FC → pct.1/pct.2 → p
write.csv(markers, "output/03_markers.csv", row.names = FALSE)

## ---- 3. singler ---------------------------------------------------- Q2 頁 43–46
library(SingleR); library(celldex)
ref  <- celldex::HumanPrimaryCellAtlasData()          # 首次下載約 1 GB，之後快取
pred <- SingleR(test = GetAssayData(gbm, layer = "data"), ref = ref,
                ## TODO ▶ SingleR 用粗標籤還是細標籤？（Q2 頁 43）
                labels = ref$____, clusters = gbm$seurat_clusters)
xval <- data.frame(cluster = rownames(pred), SingleR = pred$labels,
                   delta   = round(pred$delta.next, 2))         # delta 小 = 第一、二名分不開
print(xval)
gbm$singler <- unname(setNames(pred$labels, rownames(pred))[as.character(gbm$seurat_clusters)])   # 群標籤寫回每顆細胞（要 unname：Seurat v5 會把名字當細胞名）
# 交叉驗證：把手動 marker 的判斷填進去（看 03_dotplot_panel.png）
# 免疫 / 寡樹突 / 血管：兩證人一致即過。
# 膠質群：SingleR 會給 Astrocyte / Neural progenitor / Neurons —— 那不是答案，惡性留給 CNV。

## ---- 4. name -------------------------------------------------------- Q2 頁 47
# ★ 依照「你自己的」DotPlot 與 SingleR 結果填寫；下面只是範例對應，每份資料的編號都不同 ★
# 範例對應（用本課 seed = 1234、Seurat 5.3 跑 GBM 5k 得到的 11 群；你的編號若不同，對照 03_dotplot_panel.png 改）
new.ids <- c("0"  = "Glial (CNV pending)",            # C1QL1 / NPSR1：膠質／惡性（NPC 樣）
             "1"  = "Oligodendrocyte",                # MAG / KLK6 / HAPLN2（SingleR 給 Astrocyte 是參考集沒有寡樹突）
             "2"  = "Glial (CNV pending)",            # SAA1 / CP / CLU：星狀樣（AC 樣）——惡性與否留給 CNV
             "3"  = "Macrophage",                     # LYZ / EREG / ANPEP：單核球來源 TAM
             "4"  = "Glial (CNV pending)",            # TRIB3 / IGFBP3 / VGF：壓力／缺氧樣惡性狀態
             "5"  = "Macrophage",                     # OLR1 / CCL4L2：TAM
             "6"  = "Cycling (CNV pending)",          # TOP2A / AURKB：增殖狀態（雷三：不是型別）
             "7"  = "Microglia",                      # CX3CR1 / P2RY12 / ADORA3
             "8"  = "Pericyte / fibroblast",          # COL3A1 / ASPN / C7
             "9"  = "T cell",                         # CD3D / TRBC2 / GZMA
             "10" = "Glial (CNV pending)")            # DLX5 / DLX6 / TAC3：神經元樣（NPC 樣）
# 沒對到的群先掛 Unassigned，不要讓腳本停下來；之後回頭補
unassigned <- setdiff(levels(gbm), names(new.ids))
if (length(unassigned)) { warning("這些群還沒掛名字，先標 Unassigned：", paste(unassigned, collapse = ", "))
  new.ids <- c(new.ids, setNames(rep("Unassigned", length(unassigned)), unassigned)) }
# doublet 群：依 02 §4 的 per-cluster 表，dbl > 0.5 的整群標 DOUBLET（本次跑最高 0.33，沒有整群 doublet）
qc.tab <- read.csv("output/02_per_cluster_qc.csv")
dbl.cl <- as.character(qc.tab$cluster[qc.tab$dbl > 0.5])
if (length(dbl.cl)) new.ids[dbl.cl] <- "DOUBLET"

gbm <- RenameIdents(gbm, new.ids)
gbm$celltype <- Idents(gbm)
if (any(gbm$celltype == "DOUBLET")) gbm <- subset(gbm, subset = celltype != "DOUBLET")   # 整群移除並記錄（方法段要寫）
p <- DimPlot(gbm, label = TRUE, repel = TRUE) + NoLegend()
ggsave("output/figs/03_umap_annotated.png", p, width = 7, height = 6, dpi = 150, bg = "white")
table(gbm$celltype)

## ---- 4b. immune-subsets --------------------------------------------- Q2 頁 48–50
# 層級式註釋：免疫大陸單獨拿出來，整條流程重跑（subset 之後一定重算 HVG 與 PCA）
imm <- subset(gbm, celltype %in% c("Macrophage", "Microglia", "T cell"))
## TODO ▶ 挑幾個高變異基因？（Q2 頁 25）
imm <- NormalizeData(imm) |> FindVariableFeatures(nfeatures = ____) |> ScaleData() |> RunPCA(npcs = 30)
ElbowPlot(imm)
imm <- FindNeighbors(imm, dims = 1:20) |> FindClusters(resolution = 0.6) |> RunUMAP(dims = 1:20)
imm.panel <- c("P2RY12", "TMEM119", "CX3CR1",                       # 小膠質
               "CD163", "LYZ", "TGFBI", "S100A8", "VCAN",           # 血液來源 TAM / 單核球
               "CD3E", "CD8A", "GZMK", "CD4", "IL7R", "FOXP3",     # T
               "PDCD1", "HAVCR2", "NKG7", "GNLY", "CD1C", "MKI67")
DotPlot(imm, features = imm.panel) + RotatedAxis()
## TODO ▶ FindAllMarkers 的三個門檻各是什麼意思？（Q2 頁 41–42）
imm.markers <- FindAllMarkers(imm, only.pos = ____, min.pct = ____, logfc.threshold = ____)
# 對照 DotPlot 掛第二層名字（編號依你的資料）：
# imm.ids <- c("0" = "Microglia", "1" = "Blood-derived TAM", "2" = "Monocyte", "3" = "CD8 T",
#              "4" = "CD4 T", "5" = "Treg", "6" = "NK", "7" = "DC", "8" = "Cycling TAM")
# imm$celltype_l2 <- imm.ids[as.character(Idents(imm))]
# gbm$celltype_l2 <- gbm$celltype; gbm$celltype_l2[colnames(imm)] <- imm$celltype_l2   # 寫回全體
# 命名規範（Q2 頁 57）：celltype_l1（大陸）/ celltype_l2（型別）/ celltype_l3（狀態）/ celltype_conf
gbm$celltype_l1 <- dplyr::case_when(
  grepl("CNV pending", gbm$celltype) ~ "Malignant", gbm$celltype %in% c("Macrophage", "Microglia", "T cell") ~ "Immune",
  gbm$celltype == "Oligodendrocyte" ~ "Oligo", gbm$celltype == "Endothelial" ~ "Vascular", TRUE ~ "Stromal")

## ---- 5. malignant-states ------------------------------------------- Q2 頁 51–53
# Neftel et al. 2019 (Cell) 四種狀態的 meta-module。★ 完整基因集請用論文 Table S2（每組 50 個）；
# 這裡列出每組前 12 個作示範，練習 3-3 請你補完整。
neftel <- list(
  MES1 = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2","IFITM3"),
  MES2 = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3"),
  AC   = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1","GFAP","FABP7","BCAN","PON2"),
  OPC  = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR"),
  NPC1 = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1"),
  NPC2 = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11"))
glial <- subset(gbm, subset = celltype == "Glial (CNV pending)")
glial <- AddModuleScore(glial, features = neftel, name = "state_")
colnames(glial@meta.data)[grep("^state_", colnames(glial@meta.data))] <- names(neftel)
# 兩軸表示（Neftel Fig. 2 的做法，簡化版）：
sc <- glial@meta.data[, names(neftel)]
glial$MES <- pmax(sc$MES1, sc$MES2); glial$NPC <- pmax(sc$NPC1, sc$NPC2)
glial$axis_x <- (pmax(glial$OPC, glial$NPC) - pmax(glial$AC, glial$MES))   # 左 OPC/NPC，右 AC/MES
glial$axis_y <- ifelse(glial$axis_x > 0, glial$NPC - glial$OPC, glial$MES - glial$AC)
glial$state  <- apply(glial@meta.data[, c("MES", "AC", "OPC", "NPC")], 1, function(v) names(v)[which.max(v)])
p <- ggplot(glial@meta.data, aes(axis_x, axis_y, colour = state)) + geom_point(size = .6, alpha = .6) +
     geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_classic() +
     labs(x = "← OPC/NPC-like   AC/MES-like →", y = "state strength")
ggsave("output/figs/03_neftel_states.png", p, width = 7, height = 6, dpi = 150, bg = "white")
table(glial$state, glial$Phase)                       # 週期一起看：增殖細胞落在哪個狀態？
st <- setNames(rep(NA_character_, ncol(gbm)), colnames(gbm)); st[colnames(glial)] <- glial$state; gbm$state <- unname(st)

# 三種打分數算法（Q2 頁 53）：AddModuleScore（上）、UCell、AUCell
if (requireNamespace("UCell", quietly = TRUE)) {
  glial <- UCell::AddModuleScore_UCell(glial, features = neftel, name = "_UCell")
  print(cor(glial$MES1, glial$MES1_UCell))               # 三種算法排序通常高度一致
}

## ---- 5b. annotation-qc ---------------------------------------------- Q2 頁 54–56
round(prop.table(table(gbm$celltype, gbm$singler), 1), 2)          # 交叉表（需先把 SingleR 群標籤寫回 gbm$singler）
top5 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 5)
gbm <- ScaleData(gbm, features = unique(c(VariableFeatures(gbm), top5$gene)))   # 熱圖用的基因要先 scale 過
DoHeatmap(subset(gbm, downsample = 100), features = top5$gene, group.by = "celltype")
gbm <- BuildClusterTree(gbm, dims = 1:20); PlotClusterTree(gbm)     # 分群樹
sil <- cluster::silhouette(as.integer(Idents(gbm)), dist(Embeddings(gbm, "pca")[, 1:20]))
tapply(sil[, "sil_width"], Idents(gbm), mean)                       # < 0.1 的群存疑

## ---- 6. deliverables-and-save --------------------------------------- Q2 頁 58–59
library(patchwork)
p1 <- DimPlot(gbm, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DotPlot(gbm, features = panel, group.by = "celltype") + RotatedAxis()
p3 <- VlnPlot(gbm, features = c("nFeature_RNA", "percent.mt"), group.by = "celltype", pt.size = 0, ncol = 2)
dir.create("output/figs", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figs/03_fig1_annotation.pdf", (p1 | p2) / p3, width = 14, height = 9, bg = "white")
comp <- gbm@meta.data |> dplyr::count(celltype) |> mutate(pct = round(100 * n / sum(n), 1))
write.csv(comp, "output/03_composition.csv", row.names = FALSE)
comp
saveRDS(gbm, "output/gbm_annotated.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 3
#  3-1 把 SingleR 的 labels 換成 ref$label.fine 重跑。免疫亞群分得更細了嗎？膠質群變成什麼？
#  3-2 找出 T 細胞群，用 DotPlot 看 CD8A / CD4 / IL7R / PDCD1 / HAVCR2 / TIGIT：
#      GBM 的 T 細胞有耗竭訊號嗎？它們佔全部細胞的比例是多少？（回想 Q1：GBM 是免疫冷）
#  3-3 從 Neftel 2019 的 Table S2 取得完整的六組 meta-module（每組 50 個基因），
#      替換 neftel 清單重跑 §5。四種狀態的比例變了多少？
#  3-4 小膠質（P2RY12+）與血液來源巨噬細胞（CD163+/LYZ+）在你的資料裡分得開嗎？
#      用 FindMarkers 比較兩群，前十名基因跟文獻（Müller 2017, Genome Biol）一致嗎？
#  3-5 §4b：小膠質與血液來源 TAM 各佔免疫細胞多少？耗竭 CD8 T（PDCD1+HAVCR2+）佔 CD8 T 多少？
#  3-6 §5b：哪一群的輪廓係數最低？它跟哪一群最像？該併還是該留？用 marker 說明。
#  進階 用 AddModuleScore 給「所有細胞」算 cycling 分數（MKI67, TOP2A, CENPF, …）；
#      比較各型別的增殖比例。哪一型最高？合理嗎？
# =====================================================================
