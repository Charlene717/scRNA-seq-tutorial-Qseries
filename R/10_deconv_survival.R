# =====================================================================
# 10_deconv_survival.R — 練習腳本 10：反卷積與存活分析（MuSiC + TCGA-GBM）——把 n = 4 帶到 n = 170
#
# 對應影片：Q3 頁 73–75（§1 反卷積、§2 KM 與 Cox）
# 輸入：output/gbm4_final.rds；TCGA-GBM bulk（TCGAbiolinks 自動下載，需連網）
# 輸出：output/10_deconv_tcga_gbm.csv、output/figs/10_km_*.pdf
# 時間：TCGA 下載約 10 分鐘（僅第一次），其餘約 5 分鐘
# 這是「單細胞產生假說 → 公開世代驗證」那條路：用 4 位病人的型別比例假說，到 TCGA 的 170 份 bulk 驗證。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/gbm4_final.rds")
dir.create("output/figs", showWarnings = FALSE, recursive = TRUE)

## ---- 1. deconvolution ----------------------------------------------- Q3 頁 75
library(MuSiC); library(TCGAbiolinks); library(SummarizedExperiment); library(survival); library(survminer)
ref <- as.SingleCellExperiment(JoinLayers(gbm4))
ref$celltype_l1 <- ifelse(gbm4$malignant == "malignant", "Malignant",
                   ifelse(gbm4$celltype_author == "Immune cell", "Immune",
                   ifelse(gbm4$celltype_author == "Oligodendrocyte", "Oligo",
                   ifelse(gbm4$celltype_author == "Vascular", "Vascular", "Other"))))
q <- GDCquery(project = "TCGA-GBM", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
# ★ Windows 有 260 字元的路徑長度上限：TCGA 的檔案是「兩層 UUID + 超長檔名」，
#   放在很深的專案資料夾底下解壓會整批報「無法建立檔案」。所以下載目錄用磁碟根附近的短路徑。
#   （已下載過的檔案會自動跳過；重跑只補缺的。macOS / Linux 沒這個限制，仍建議獨立資料夾。）
GDC_DIR <- if (.Platform$OS.type == "windows") "C:/GDCdata" else "~/GDCdata"
dir.create(GDC_DIR, showWarnings = FALSE, recursive = TRUE)
GDCdownload(q, directory = GDC_DIR, files.per.chunk = 50)
bulk <- GDCprepare(q, directory = GDC_DIR)                                   # 約 170 份，含臨床欄位
bulk.mtx <- assay(bulk, "unstranded"); rownames(bulk.mtx) <- rowData(bulk)$gene_name
bulk.mtx <- bulk.mtx[!duplicated(rownames(bulk.mtx)), ]
common <- intersect(rownames(bulk.mtx), rownames(ref))
est <- music_prop(bulk.mtx = bulk.mtx[common, ], sc.sce = ref[common, ], clusters = "celltype_l1", samples = "patient")
prop <- as.data.frame(est$Est.prop.weighted); write.csv(prop, "output/10_deconv_tcga_gbm.csv")
## ---- 2. survival ---------------------------------------------------- Q3 頁 75
# 存活：免疫（巨噬）比例上下半
clin <- as.data.frame(colData(bulk))[rownames(prop), ]
clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up) / 30.4
clin$event <- as.integer(clin$vital_status == "Dead")
clin$mac_hi <- prop$Immune > median(prop$Immune)
fit <- survfit(Surv(time, event) ~ mac_hi, data = clin)
p <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE, xlab = "Months"); pdf("output/figs/10_km_macrophage.pdf", 7, 6); print(p); dev.off()
coxph(Surv(time, event) ~ mac_hi + age_at_index, data = clin)               # 實際分析還要調整 MGMT、IDH
sessionInfo()

# =====================================================================
# ▶ 練習 10
#  10-1 用 Malignant 比例（純度）分組做 KM，跟巨噬比例的結果方向相同嗎？兩者相關係數多少？
#  10-2 在 Cox 模型加入 age 之後，mac_hi 的 HR 變化多少？這代表什麼？
#  10-3 參考組（sc 端）把 Other 拆成 Astrocyte / OPC / Neuron 重跑：Immune 的估計比例變多少？
#       反卷積對參考的敏感度告訴你什麼？
#  進階 用 BayesPrism 重做 §1，比較兩種反卷積估的巨噬比例（相關係數、Bland–Altman 圖）。
# =====================================================================
