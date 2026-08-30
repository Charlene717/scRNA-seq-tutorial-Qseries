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
pa <- gbm4@meta.data |> mutate(hyp = act.w["Hypoxia", colnames(gbm4)]) |>
      dplyr::filter(malignant == "malignant") |> group_by(patient, tissue) |> summarise(hyp = mean(hyp), .groups = "drop") |>
      tidyr::pivot_wider(names_from = tissue, values_from = hyp)
t.test(pa$Tumor, pa$Periphery, paired = TRUE); write.csv(pa, "output/tables/09_progeny_hypoxia_by_sample.csv", row.names = FALSE)
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
            assay = "scenic", group.by = "cc_label")
  DefaultAssay(gbm4) <- "RNA"
}

sessionInfo()

# =====================================================================
# ▶ 練習 9
#  9-1 Hypoxia 活性在核心 vs 邊緣的配對檢定 p 值多少？跟 06a 的 GSEA（HALLMARK_HYPOXIA）方向一致嗎？
#  9-2 EGFR 與 JAK-STAT 的活性在哪種細胞最高？用 FeaturePlot 對照 celltype_author 說明。
#  9-3 把 §1 的「病人 × 部位平均 → 配對 t 檢定」流程套到 TGFb：結論是什麼？單位為什麼是病人不是細胞？
#  進階 跑 pySCENIC（§2 註解的三步），比較惡性細胞的 SOX2(+)/OLIG2(+) 與 TAM 的 SPI1(+)/CEBPB(+) 活性分布。
# =====================================================================
