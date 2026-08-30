# =====================================================================
# 03_annotate.R — 練習腳本 3：譜系標誌基因、marker、SingleR、命名、惡性狀態分數
#
# 對應影片：Q2 頁 38–59（§1 譜系標誌與 marker 面板、§2 FindAllMarkers 與篩選、§3 SingleR、§4 命名、§4b 免疫亞群、§5 Neftel 分數與三種算法、§5b 品質檢查、§6 交付與存檔）
# 輸入：output/rds/02_gbm_clustered.rds（02_cluster.R）
# 輸出：output/rds/03_gbm_annotated.rds、output/tables/03_markers.csv（另有 03_immune_markers、03_composition）
# 時間：約 5–10 分鐘（SingleR 首次下載參考集約 1 GB，之後有快取）
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)   # 圖檔資料夾先建好，後面 ggsave 才不會失敗
gbm <- readRDS("output/rds/02_gbm_clustered.rds")
Idents(gbm) <- "seurat_clusters"

# 保險：02 的定案欄位應該是 res 0.5。對不上代表 02 的穩定性檢查（FindClusters）覆寫了 seurat_clusters
if (!is.null(gbm$RNA_snn_res.0.5) &&
    !identical(as.character(gbm$seurat_clusters), as.character(gbm$RNA_snn_res.0.5)))
  stop("seurat_clusters 不等於 RNA_snn_res.0.5：回 02_cluster.R 確認 seurat_clusters 是在\n",
       "  所有 FindClusters()（含 seed42 / k30 穩定性檢查）之後才寫回的。", call. = FALSE)
cat("分群數：", nlevels(Idents(gbm)), "群\n")

## ---- 1. gates-and-panel -------------------------------------------- Q2 頁 38–40
# 四個譜系標誌基因：先劃四個主要類群
gates <- c("PTPRC",   # 免疫
           "SOX2",    # 膠質／惡性（正常星狀、OPC 也會亮）
           "MBP",     # 寡樹突
           "PECAM1")  # 血管內皮
gates <- intersect(gates, rownames(gbm))              # 資料裡沒有的基因先剔除，畫圖才不會報錯
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
panel <- lapply(panel, intersect, rownames(gbm)); panel <- panel[lengths(panel) > 0]   # 同上
p <- DotPlot(gbm, features = panel, cluster.idents = TRUE) + RotatedAxis() +
     theme(axis.text.x = element_text(size = 8))
ggsave("output/figs/03_dotplot_panel.png", p, width = 16, height = 6, dpi = 150, bg = "white")
# 看圖：沿對角線一塊塊亮起來嗎？PTPRC 亮在幾群（那是免疫類群）？MKI67 疊在哪幾群上？

## ---- 2. markers ---------------------------------------------------- Q2 頁 41–42
markers <- FindAllMarkers(gbm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top5 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 5) |>
        select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
print(top5, n = 60)                                   # 看表順序：log2FC → pct.1/pct.2 → p
write.csv(markers, "output/tables/03_markers.csv", row.names = FALSE)

## ---- 3. singler ---------------------------------------------------- Q2 頁 43–46
library(SingleR); library(celldex)
ref  <- celldex::HumanPrimaryCellAtlasData()          # 首次下載約 1 GB，之後快取
pred <- SingleR(test = GetAssayData(gbm, layer = "data"), ref = ref,
                labels = ref$label.main, clusters = gbm$seurat_clusters)
xval <- data.frame(cluster = rownames(pred), SingleR = pred$labels,
                   delta   = round(pred$delta.next, 2))         # delta 小 = 第一、二名分不開
print(xval)
gbm$singler <- unname(setNames(pred$labels, rownames(pred))[as.character(gbm$seurat_clusters)])   # 群標籤寫回每顆細胞（要 unname：Seurat v5 會把名字當細胞名）
# 交叉驗證：把手動 marker 的判斷填進去（看 03_dotplot_panel.png）
# 免疫 / 寡樹突 / 血管：兩項獨立證據一致即可接受。
# 膠質群：SingleR 會給 Astrocyte / Neural progenitor / Neurons —— 那不是答案，是提醒（參考集裡沒有 GBM 惡性細胞）。
# 惡性與否需要基因體層級的證據（CNV 推斷、突變基因型）加上跨病人專屬性；單一病人的資料湊不齊，
# 所以本腳本的膠質群一律停在 (undetermined)，這是這份資料能給的最後結論。

## ---- 4. name -------------------------------------------------------- Q2 頁 47
# ★ 依照「你自己的」DotPlot 與 SingleR 結果填寫；下面只是範例對應，每份資料的編號都不同 ★
# 範例對應（seed = 1234、Seurat 5.3、npc = 25、res 0.5、k = 20 跑 GBM 5k 得到的 13 群）
# ※ 這是「範例」，不是答案：務必先看你自己的 03_dotplot_panel.png 與上面的 xval 表再定案。
new.ids <- c("0"  = "Glial (undetermined)",           # C1QL1 / NPSR1：膠質／惡性（NPC 樣）
             "1"  = "Oligodendrocyte",                # MAG / KLK6 / HAPLN2（SingleR 給 Astrocyte 是參考集沒有寡樹突）
             "2"  = "Glial (undetermined)",           # SAA1 / CP / CLU：星狀樣（AC 樣）——惡性與否待多重證據判定
             "3"  = "Glial (undetermined)",           # TRIB3 / IGFBP3 / VGF：壓力／缺氧樣惡性狀態
             "4"  = "Macrophage",                     # OLR1 / CCL4L2：TAM
             "5"  = "Macrophage",                     # MRC1 / LYVE1 / SDS：血管周巨噬細胞
             "6"  = "Macrophage",                     # FCN1 / EREG / ANPEP：單核球剛分化來的 TAM
             "7"  = "Cycling (undetermined)",         # TOP2A / AURKB：增殖狀態（錯誤三：不是型別）
             "8"  = "Microglia",                      # CX3CR1 / P2RY12 / ADORA3
             "9"  = "Pericyte / fibroblast",          # COL3A1 / ASPN / C7
             "10" = "T cell",                         # CD3D / TRBC2 / GZMA
             "11" = "Glial (undetermined)",           # DLX5 / DLX6 / TAC3：神經元樣（NPC 樣）
             "12" = "Glial (undetermined)")           # ← 02 §4b 三項診斷全中（見該節判讀），是真 doublet 群；下面會自動標成 DOUBLET 並移除
# 沒對到的群先掛 Unassigned，不要讓腳本停下來；之後回頭補
unassigned <- setdiff(levels(gbm), names(new.ids))
if (length(unassigned)) { warning("這些群還沒命名，先標 Unassigned：", paste(unassigned, collapse = ", "))
  new.ids <- c(new.ids, setNames(rep("Unassigned", length(unassigned)), unassigned)) }
# doublet 群：依 02 §4 的 per-cluster 表，dbl > 0.5 的群列為「整群 doublet」候選。
# ⚠ 這條閾值只是「篩出候選」，不是判定。整群移除是不可逆的決定，要三項證據一起看：
#   ① nCount 中位數是否接近「兩個來源群相加」（真 doublet 帶著兩顆細胞的 RNA；
#      同型 doublet 才是單一群的兩倍，異型 doublet 要拿兩個親代群的中位數相加來比）
#   ② 同一顆細胞內是否同時表現互斥的譜系標記（PTPRC + SOX2、MBP + CD3E …）
#   ③ 該群有沒有「專屬」marker，還是只是另外兩群 marker 的聯集
#   三項都成立 → 整群移除並在方法段記錄；若有專屬 marker 且 nCount 沒翻倍，
#   那是被工具誤判的真實族群（腫瘤裡 RNA 量大的惡性細胞最常被誤判）——不要刪。
# 這三項的診斷程式與判讀在 02_cluster.R §4b（對應該節輸出的 ②③④），跑 03 之前先看過那一段的輸出，
# 確認候選群真的該刪，再往下執行。
qc.tab <- read.csv("output/tables/02_per_cluster_qc.csv")
dbl.cl <- as.character(qc.tab$cluster[qc.tab$dbl > 0.5])
if (length(dbl.cl)) new.ids[dbl.cl] <- "DOUBLET"

gbm <- RenameIdents(gbm, new.ids)
gbm$celltype <- Idents(gbm)
if (any(gbm$celltype == "DOUBLET")) gbm <- subset(gbm, subset = celltype != "DOUBLET")   # 整群移除並記錄（方法段要寫）
gbm$celltype <- droplevels(gbm$celltype); Idents(gbm) <- "celltype"   # 清掉空的 level，圖例才不會多出空類別
p <- DimPlot(gbm, label = TRUE, repel = TRUE) + NoLegend()
ggsave("output/figs/03_umap_annotated.png", p, width = 7, height = 6, dpi = 150, bg = "white")
table(gbm$celltype)

## ---- 4b. immune-subsets --------------------------------------------- Q2 頁 48–50
# 層級式註釋：免疫類群單獨拿出來，整條流程重跑（subset 之後一定重算 HVG 與 PCA）
imm <- subset(gbm, celltype %in% c("Macrophage", "Microglia", "T cell"))
imm <- NormalizeData(imm) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 30)
# 看一下印出來的載荷：PC1 是髓系 vs T（對的），但 PC2 的一端是 MAG / PLP1 / KLK6 / CLDN11（寡樹突）、
#   PC5 的一端是 FABP7 / PTPRZ1 / GFAP / BCAN（膠質）。免疫子集裡出現這些基因有兩種可能：
#   環境 RNA（Q2 講的 SoupX），或 TAM 真的吞了髓鞘與腫瘤碎片（myelin-laden macrophage，文獻有記載）。
#   分辨方法：SoupX 校正後再看一次——訊號整片消失是環境 RNA，只集中在特定亞群才是生物學。
ggsave("output/figs/03_imm_elbow.png", ElbowPlot(imm), width = 6, height = 4, dpi = 150, bg = "white")
imm <- FindNeighbors(imm, dims = 1:20) |> FindClusters(resolution = 0.6) |> RunUMAP(dims = 1:20)
imm.panel <- c("P2RY12", "TMEM119", "CX3CR1",                       # 小膠質
               "CD163", "LYZ", "TGFBI", "S100A8", "VCAN",           # 血液來源 TAM / 單核球
               "CD3E", "CD8A", "GZMK", "CD4", "IL7R", "FOXP3",     # T
               "PDCD1", "HAVCR2", "NKG7", "GNLY", "CD1C", "MKI67")
imm.panel <- intersect(imm.panel, rownames(imm))
p <- DotPlot(imm, features = imm.panel) + RotatedAxis()
ggsave("output/figs/03_imm_dotplot.png", p, width = 11, height = 5, dpi = 150, bg = "white")
imm.markers <- FindAllMarkers(imm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(imm.markers, "output/tables/03_immune_markers.csv", row.names = FALSE)
# 對照 DotPlot 標上第二層名稱（編號依你的資料）：
# imm.ids <- c("0" = "Microglia", "1" = "Blood-derived TAM", "2" = "Monocyte", "3" = "CD8 T",
#              "4" = "CD4 T", "5" = "Treg", "6" = "NK", "7" = "DC", "8" = "Cycling TAM")
# imm$celltype_l2 <- imm.ids[as.character(Idents(imm))]
# gbm$celltype_l2 <- gbm$celltype; gbm$celltype_l2[colnames(imm)] <- imm$celltype_l2   # 寫回全體
# 命名規範（Q2 頁 57）：celltype_l1（主要類群）/ celltype_l2（型別）/ celltype_l3（狀態）/ celltype_conf
# 注意：這裡「不」直接寫 Malignant。marker 只能定譜系：膠質瘤惡性細胞的正常對應細胞（星狀、OPC）
#       就在同一塊組織裡，表現量高度重疊，所以沒有任何一組 marker 能區分惡性與正常膠質。
#       在證據到位前就叫 Malignant，就是錯誤二（用型別標籤預設了結論）。
gbm$celltype_l1 <- dplyr::case_when(
  grepl("undetermined", gbm$celltype)                       ~ "Glial (undetermined)",
  gbm$celltype %in% c("Macrophage", "Microglia", "T cell")  ~ "Immune",
  gbm$celltype == "Oligodendrocyte"                         ~ "Oligo",
  grepl("Endothelial|Pericyte|fibroblast", gbm$celltype)    ~ "Vascular / stromal",
  TRUE                                                      ~ "Other")
table(gbm$celltype_l1)

## ---- 5. malignant-states ------------------------------------------- Q2 頁 51–53
# Neftel et al. 2019 (Cell) 四種狀態的 meta-module。★ 完整基因集請用論文 Table S2（每組 50 個）；
# 這裡列出每組前 12 個作示範，練習 3-3 請你補完整。
# 注意方向性：這四種狀態是 Neftel 在「已經確認是惡性」的細胞裡定義出來的。分數只說明這群細胞
#             偏向哪一種狀態，不能反過來當成惡性的證據——正常的星狀細胞與 OPC 一樣會在
#             AC / OPC 上拿到高分。跑完這一節，膠質群仍然是 (undetermined)。
neftel <- list(
  MES1 = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2","IFITM3"),
  MES2 = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3"),
  AC   = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1","GFAP","FABP7","BCAN","PON2"),
  OPC  = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR"),
  NPC1 = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1"),
  NPC2 = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11"))
glial <- subset(gbm, subset = celltype == "Glial (undetermined)")
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
## >>> 參考答案 ------------------------------------------------------
# 本課的結果（2,379 顆膠質細胞，實跑驗證）：
#   MES 1,562（66%）／AC 537（23%）／NPC 258（11%）／OPC 22（0.9%）
#   MES 這麼高要留意：這裡的 MES1/MES2 只取了前 12 個基因，其中 VIM、CD44、ANXA1、ANXA2、
#     MT2A、SOD2 都是到處都表現的基因，本來就容易拿高分；解離壓力也會推高 MES 訊號。
#     是這顆腫瘤真的 MES 主導，還是截斷基因集造成的？練習 3-3 換成完整 50 個基因就能分辨。
#   週期一起看：NPC 有 71% 落在 S/G2M、AC 45%、MES 只有 26%（OPC 只有 22 顆，比例僅供參考）。
#     增殖集中在 NPC 樣狀態，跟 Neftel 的觀察一致。
## <<< 參考答案
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
p <- DoHeatmap(subset(gbm, downsample = 100), features = unique(top5$gene), group.by = "celltype")
ggsave("output/figs/03_heatmap_top5.png", p, width = 12, height = 10, dpi = 150, bg = "white")
gbm <- BuildClusterTree(gbm, dims = 1:20)                           # 分群樹（base 繪圖，用 png() 存）
png("output/figs/03_cluster_tree.png", width = 1200, height = 900, res = 150); PlotClusterTree(gbm); dev.off()

# 輪廓係數：dist() 是 O(n²)，細胞多要先抽樣，否則記憶體會爆
set.seed(1234)
cells.sil <- if (ncol(gbm) > 20000) sample(colnames(gbm), 20000) else colnames(gbm)
sil <- cluster::silhouette(as.integer(Idents(gbm)[cells.sil]),
                           dist(Embeddings(gbm, "pca")[cells.sil, 1:20]))
tapply(sil[, "sil_width"], Idents(gbm)[cells.sil], mean)            # < 0.1 的群存疑

## ---- 6. deliverables-and-save --------------------------------------- Q2 頁 58–59
library(patchwork)
p1 <- DimPlot(gbm, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DotPlot(gbm, features = panel, group.by = "celltype") + RotatedAxis()
p3 <- VlnPlot(gbm, features = c("nFeature_RNA", "percent.mt"), group.by = "celltype", pt.size = 0, ncol = 2)
ggsave("output/figs/03_fig1_annotation.pdf", (p1 | p2) / p3, width = 14, height = 9, bg = "white")
## >>> 參考答案 ------------------------------------------------------
# 本課這份 GBM 5k 的結果（實跑驗證；13 群、移除整群 doublet 的 cluster 12 後剩 5,130 顆）：
#   Glial (undetermined) 2,379（46.4%）／Macrophage 1,178（23.0%）／Oligodendrocyte 850（16.6%）
#   Cycling (undetermined) 276（5.4%）／Microglia 171（3.3%）／Pericyte / fibroblast 140（2.7%）
#   T cell 136（2.7%）—— T 細胞只有 2.7%，正是 Q1 說的「GBM 是免疫冷腫瘤」
# §5b 的輪廓係數：T cell 0.69、Oligodendrocyte 0.61、Microglia 0.49、Glial 0.32、
#   Cycling 0.24、Macrophage 0.17、Pericyte / fibroblast 0.16。後兩個偏低，但都有理由：
#   Macrophage 是 cluster 4/5/6 合併的（本來就不會是一顆緊緻的球），
#   Pericyte / fibroblast 只有 140 顆又貼著血管——低分不等於該併，要有 marker 依據才動。
## <<< 參考答案
comp <- gbm@meta.data |> dplyr::count(celltype) |> mutate(pct = round(100 * n / sum(n), 1))
write.csv(comp, "output/tables/03_composition.csv", row.names = FALSE)
comp
saveRDS(gbm, "output/rds/03_gbm_annotated.rds")
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
#  3-7 §4b 的 imm PCA 印出來，PC2 的一端是 MAG / PLP1 / KLK6 / CLDN11（寡樹突）、
#      PC5 的一端是 FABP7 / PTPRZ1 / GFAP / BCAN（膠質）。免疫細胞裡為什麼會有這些基因？
#      兩個假設：(a) 環境 RNA；(b) TAM 真的吞了髓鞘與腫瘤碎片（myelin-laden macrophage）。
#      用資料分辨，不要用直覺：
#        ① 算每個免疫亞群裡 MBP / PLP1 陽性細胞的比例與平均表現（FetchData + aggregate）；
#        ② 看 output/tables/03_immune_markers.csv 裡這兩個基因的 pct.1 與 pct.2。
#      判準：均勻散在所有亞群、pct.2 也高 → 環境 RNA；只集中在一兩個亞群，
#            且該亞群有自己的專屬 marker 撐著 → 生物學。
#      寫下你的結論，以及「什麼結果會推翻它」——後面這句才是這題的重點。
#      加分：跑 SoupX 校正後重做一次。訊號是整片消失，還是只降一點？
#  進階 用 AddModuleScore 給「所有細胞」算 cycling 分數（MKI67, TOP2A, CENPF, …）；
#      比較各型別的增殖比例。哪一型最高？合理嗎？
# =====================================================================
