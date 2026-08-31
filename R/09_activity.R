# =====================================================================
# 09_activity.R — 練習腳本 9：路徑與轉錄因子活性（decoupleR / PROGENy、SCENIC）
#
# 對應影片：Q3 頁 70–72（§1 PROGENy 路徑活性與條件比較、§2 SCENIC regulon）
# 輸入：output/rds/06_gbm4_final.rds
# 輸出：output/figs/09_*.pdf、output/tables/09_progeny_*.csv
# 時間：decoupleR 約 2 分鐘；SCENIC（pySCENIC，Python）數小時，為選配
# PROGENy 問「哪條訊號路徑活著」（footprint 基因）；SCENIC 問「哪個轉錄因子在驅動」（regulon）。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/rds/06_gbm4_final.rds")
# 04 的 IntegrateLayers 之後，RNA 是「按病人分層」的（data.BT_S1、data.BT_S2…），
# 沒有單一的 "data" 層；少了下面這一行，LayerData(gbm4, layer = "data") 會直接失敗。
# 凡是要一次拿到全部細胞的表現量矩陣（06b、07、09、10）都需要先 JoinLayers。
gbm4[["RNA"]] <- JoinLayers(gbm4[["RNA"]])
if (!"data" %in% Layers(gbm4[["RNA"]])) gbm4 <- NormalizeData(gbm4)
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. pathway-activity（decoupleR / PROGENy）---------------------- Q3 頁 71–72
library(decoupleR)
net <- get_progeny(organism = "human", top = 500)
mat <- as.matrix(LayerData(gbm4, layer = "data"))
act <- run_mlm(mat = mat, net = net, .source = "source", .target = "target", .mor = "weight", minsize = 5)
act.w <- act |> tidyr::pivot_wider(id_cols = source, names_from = condition, values_from = score) |>
         tibble::column_to_rownames("source") |> as.matrix()
gbm4[["progeny"]] <- CreateAssayObject(act.w[, colnames(gbm4)])
DefaultAssay(gbm4) <- "progeny"
FeaturePlot(gbm4, features = c("Hypoxia", "JAK-STAT", "EGFR", "TGFb"), ncol = 4); ggsave("output/figs/09_progeny_umap.pdf", width = 16, height = 4, bg = "white")
# 條件比較：惡性細胞的路徑活性 → 病人 × 部位 平均 → 配對檢定（單位 = 病人）
# 「n 先於 p」在這裡要問兩次：幾位病人配得成對，以及每一格是幾顆細胞平均出來的。
# 這份資料的邊緣樣本幾乎都是正常腦組織，惡性細胞很少——只看平均值看不出這件事，
# 一格若只有幾顆細胞，那個平均不穩，方向翻過來完全正常，不是生物學。
md <- gbm4@meta.data |> mutate(hyp = act.w["Hypoxia", colnames(gbm4)]) |>
      dplyr::filter(malignant == "malignant")
pa <- md |> group_by(patient, tissue) |> summarise(hyp = mean(hyp), n = dplyr::n(), .groups = "drop") |>
      tidyr::pivot_wider(names_from = tissue, values_from = c(hyp, n))
for (v in c("n_Tumor", "n_Periphery")) pa[[v]][is.na(pa[[v]])] <- 0
cat("\n== 每位病人各部位：Hypoxia 活性平均與惡性細胞數 ==\n"); print(as.data.frame(pa))
ok <- pa$n_Tumor > 0 & pa$n_Periphery > 0        # 兩側都有惡性細胞才配得成對
MINCELL <- 20                                    # 少於這個數的那一格，平均值不穩，方向翻過來很正常
thin <- ok & (pa$n_Tumor < MINCELL | pa$n_Periphery < MINCELL)
cat("進入配對檢定的病人數 n =", sum(ok), "／ 共", nrow(pa), "位\n")
if (any(thin)) cat("注意：", paste(pa$patient[thin], collapse = "、"),
                   "有一側不到", MINCELL, "顆惡性細胞，那一格的平均只是幾顆細胞的平均\n")
print(t.test(pa$hyp_Tumor[ok], pa$hyp_Periphery[ok], paired = TRUE))
# 本例的邊緣側惡性細胞數：BT_S1 = 1、BT_S2 = 13、BT_S4 = 17、BT_S6 = 0。
# 三個配得成對的病人，邊緣那一格全都不到 20 顆；BT_S1 那個 3.80 是「一顆細胞」的值。
# 所以這裡有兩層問題，而且第一層就足以判出局：
#   ① 進到檢定裡的東西不可信——一顆細胞的平均不是那個病人邊緣的缺氧活性。
#      BT_S1 之所以方向相反，最合理的解釋就是這個，不是生物學。
#   ② 就算數字可信，n = 3 也測不到東西：p = 0.834、平均差 0.30、95% CI [-5.19, 5.80]。
#      區間寬到從 -5 跨到 +5，意思是「這個檢定看不出來」，不是「兩個部位沒有差別」——
#      不顯著與沒差異是兩件事，第 47 頁講的功效就是在講這個。
# 對照組：第 42 頁的 GSEA 看得到「核心缺氧」，那是免疫細胞、8 個樣本、幾千個基因的排名；
#         第 67 頁的軌跡看得到缺氧基因上升，那是惡性細胞內部沿 pseudotime 的連續變化。
#         同一套生物學，換個統計單位就從看得到變成看不到——差別不在生物學。
write.csv(pa, "output/tables/09_progeny_hypoxia_by_sample.csv", row.names = FALSE)
DefaultAssay(gbm4) <- "RNA"
## ---- 2. scenic（選配）----------------------------------------------- Q3 頁 72
# SCENIC（pySCENIC，Python）：R 端匯出 loom，跑完讀回
#   library(SeuratDisk); SaveLoom(gbm4, "output/rds/09_gbm4.loom")   # 或 loomR / anndata
#   pyscenic grn output/09_gbm4.loom hs_hgnc_tfs.txt -o output/09_adj.csv --num_workers 8
#   pyscenic ctx output/09_adj.csv hg38_*.feather --annotations_fname motifs-v10.tbl --expression_mtx_fname output/09_gbm4.loom -o output/09_reg.csv
#   pyscenic aucell output/09_gbm4.loom output/09_reg.csv -o output/09_auc.loom
if (file.exists("output/tables/09_scenic_auc.csv")) {                      # regulon × cell（自 auc.loom 匯出）
  auc <- read.csv("output/tables/09_scenic_auc.csv", row.names = 1, check.names = FALSE)
  gbm4[["scenic"]] <- CreateAssayObject(as.matrix(auc)[, colnames(gbm4)])
  DoHeatmap(subset(gbm4, downsample = 100), features = c("SOX2(+)", "OLIG2(+)", "SOX10(+)", "SPI1(+)", "CEBPB(+)", "TCF7(+)", "ERG(+)"),
            assay = "scenic", group.by = "type")   # cc_label 是 07 建的，06 的物件裡沒有
  DefaultAssay(gbm4) <- "RNA"
}

sessionInfo()

# =====================================================================
# ▶ 練習 9
#  9-1 先看每一格的惡性細胞數，再決定那個 p 值值不值得讀。跟 06a 的 GSEA（HALLMARK_HYPOXIA）比，
#      為什麼同一份資料、同一個生物學，一邊看得到、一邊看不到？（提示：統計單位與每格的細胞數）
#  9-2 EGFR 與 JAK-STAT 的活性在哪種細胞最高？用 FeaturePlot 對照 celltype_author 說明。
#  9-3 把 §1 的「病人 × 部位平均 → 配對 t 檢定」流程套到 TGFb：結論是什麼？單位為什麼是病人不是細胞？
#  進階 跑 pySCENIC（§2 註解的三步），比較惡性細胞的 SOX2(+)/OLIG2(+) 與 TAM 的 SPI1(+)/CEBPB(+) 活性分布。
# =====================================================================
