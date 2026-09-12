# =====================================================================
# 09_activity.R — 練習腳本 9：路徑與轉錄因子活性（decoupleR / PROGENy、SCENIC）
#
# 對應影片：Q3 頁 76–78（§1 PROGENy 路徑活性與條件比較、§2 SCENIC regulon）
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

## ---- 1. pathway-activity（decoupleR / PROGENy）---------------------- Q3 頁 77–78
library(decoupleR)
net <- get_progeny(organism = "human", top = 500)
mat <- as.matrix(LayerData(gbm4, layer = "data"))
act <- run_mlm(mat = mat, net = net, .source = "source", .target = "target", .mor = "weight", minsize = 5)
act.w <- act |> tidyr::pivot_wider(id_cols = source, names_from = condition, values_from = score) |>
         tibble::column_to_rownames("source") |> as.matrix()
# 用 data= 而不是位置參數。CreateAssayObject() 的第一個參數是 counts，
# 但 PROGENy 活性是推出來的連續分數（會有負值），不是計數；放進 counts 之後
# 任何「以為拿到的是 counts」的後續函式都可能默默算錯。
gbm4[["progeny"]] <- CreateAssayObject(data = act.w[, colnames(gbm4)])
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
MINCELL <- 20                                    # 少於這個數的那一格，平均值不穩，方向翻過來很正常
paired <- pa$n_Tumor > 0 & pa$n_Periphery > 0    # 兩側都有惡性細胞才配得成對
ok     <- paired & pa$n_Tumor >= MINCELL & pa$n_Periphery >= MINCELL   # 還要每一格都夠厚
thin   <- paired & !ok
cat("配得成對的病人數 =", sum(paired), "／ 其中每格都 ≥", MINCELL, "顆的 =", sum(ok),
    "／ 共", nrow(pa), "位\n")
if (any(thin)) cat("排除：", paste(pa$patient[thin], collapse = "、"),
                   "有一側不到", MINCELL, "顆惡性細胞，那一格的平均只是幾顆細胞的平均\n")
# ★ 門檻要真的擋得住。前一版只把不足的病人印成警告，t.test 照跑照印 p 值——
#   那等於嘴上說不可信、手上還是產出了一個正式的統計結果。這裡改成：不夠就不做推論，
#   只畫描述性的配對圖。這正是這門課要教的：軟體跑得動，不代表這個分析該跑。
if (sum(ok) >= 3) {
  print(t.test(pa$hyp_Tumor[ok], pa$hyp_Periphery[ok], paired = TRUE))
} else {
  cat("\n>> 每格 ≥", MINCELL, "顆的配對病人只有", sum(ok), "位，不做配對檢定。\n",
      "   下面只畫描述性的配對變化圖；報告裡要寫的是「資料條件不允許做這個比較」，\n",
      "   而不是一個沒有意義的 p 值。\n")
  draw_pairs <- function() {                       # 螢幕與檔案各畫一次，不然這張圖只活在 RStudio 裡
    matplot(t(as.matrix(pa[paired, c("hyp_Tumor", "hyp_Periphery")])), type = "b", pch = 16,
            xaxt = "n", ylab = "PROGENy Hypoxia", xlab = "",
            main = sprintf("Descriptive only (n = %d pairs, all thin)", sum(paired)))
    axis(1, at = 1:2, labels = c("Tumor", "Periphery"))
    legend("topright", legend = sprintf("%s (peri n=%d)", pa$patient[paired], pa$n_Periphery[paired]),
           col = seq_len(sum(paired)), lty = 1, pch = 16, bty = "n", cex = 0.8)
  }
  draw_pairs()
  pdf("output/figs/09_progeny_hypoxia_paired.pdf", 5, 4); draw_pairs(); dev.off()
}
# 本例的邊緣側惡性細胞數：BT_S1 = 1、BT_S2 = 13、BT_S4 = 17、BT_S6 = 0。
# 三個配得成對的病人，邊緣那一格全都不到 20 顆；BT_S1 那個 3.80 是「一顆細胞」的值。
# 所以這裡有兩層問題，而且第一層就足以判出局：
#   ① 進到檢定裡的東西不可信——一顆細胞的平均不是那個病人邊緣的缺氧活性。
#      BT_S1 之所以方向相反，最合理的解釋就是這個，不是生物學。
#   ② 就算數字可信，n = 3 也測不到東西：p = 0.834、平均差 0.30、95% CI [-5.19, 5.80]。
#      （這三個數字是把上面那三位硬跑一次配對 t 檢定得到的——腳本刻意不跑，列在這裡
#        是要讓你看到「就算跑了也沒用」。想自己驗證：t.test(pa$hyp_Tumor[paired], pa$hyp_Periphery[paired], paired = TRUE)）
#      區間寬到從 -5 跨到 +5，意思是「這個檢定看不出來」，不是「兩個部位沒有差別」——
#      不顯著與沒差異是兩件事，第 52 頁講的功效就是在講這個。
# 對照組：第 47 頁的 GSEA 看得到「核心缺氧」，那是免疫細胞、8 個樣本、幾千個基因的排名；
#         第 72 頁的軌跡看得到缺氧基因上升，那是惡性細胞內部沿 pseudotime 的連續變化。
#         同一套生物學，換個統計單位就從看得到變成看不到——差別不在生物學。
write.csv(pa, "output/tables/09_progeny_hypoxia_by_sample.csv", row.names = FALSE)
DefaultAssay(gbm4) <- "RNA"
## ---- 2. scenic（選配）----------------------------------------------- Q3 頁 78
# SCENIC（pySCENIC，Python）：R 端匯出 loom，跑完讀回
#   library(SeuratDisk); SaveLoom(gbm4, "output/rds/09_gbm4.loom")   # 或 loomR / anndata
#   pyscenic grn output/09_gbm4.loom hs_hgnc_tfs.txt -o output/09_adj.csv --num_workers 8
#   pyscenic ctx output/09_adj.csv hg38_*.feather --annotations_fname motifs-v10.tbl --expression_mtx_fname output/09_gbm4.loom -o output/09_reg.csv
#   pyscenic aucell output/09_gbm4.loom output/09_reg.csv -o output/09_auc.loom
if (file.exists("output/tables/09_scenic_auc.csv")) {                      # regulon × cell（自 auc.loom 匯出）
  auc <- read.csv("output/tables/09_scenic_auc.csv", row.names = 1, check.names = FALSE)
  gbm4[["scenic"]] <- CreateAssayObject(data = as.matrix(auc)[, colnames(gbm4)])   # AUC 也是分數，同上
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
#  9-3 把 §1 的「病人 × 部位平均 → 門檻檢查 → 配對檢定」流程套到 TGFb。
#      注意：細胞數跟 Hypoxia 那一輪完全一樣，所以門檻一樣會擋下來——這就是答案的一半。
#      另一半：如果你把 MINCELL 調低讓它跑得出 p 值，那個 p 值可以寫進論文嗎？單位為什麼是病人不是細胞？
#  進階 跑 pySCENIC（§2 註解的三步），比較惡性細胞的 SOX2(+)/OLIG2(+) 與 TAM 的 SPI1(+)/CEBPB(+) 活性分布。
# =====================================================================
