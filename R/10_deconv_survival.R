# =====================================================================
# 10_deconv_survival.R — 練習腳本 10：反卷積與存活分析（MuSiC + TCGA-GBM）——把 n = 4 帶到數百位病人
#
# 對應影片：Q3 頁 79–82（§1 反卷積、§2 KM 與 Cox）
# 輸入：output/rds/06_gbm4_final.rds；TCGA-GBM Bulk（TCGAbiolinks 自動下載，需連網）
# 輸出：output/tables/10_deconv_tcga_gbm.csv、output/figs/10_km_*.pdf
# 時間：TCGA 下載約 10 分鐘（僅第一次），其餘約 5 分鐘
# 這是「單細胞產生假說 → 公開世代驗證」那條路：用 4 位病人的型別比例假說，到 TCGA 的 Bulk 世代驗證。
# 注意：TCGA 下載回來的是「檔案」不是「病人」——同一位病人常有兩三份 aliquot，還混著復發與正常組織。
# 進統計之前要先整理成「一位病人一筆、只留原發腫瘤」，否則等於把同一個死亡事件重複計算。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/rds/06_gbm4_final.rds")
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- 1. deconvolution ----------------------------------------------- Q3 頁 82
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
# GDCdownload 的 api 方法是「抓一包 tar → 在『當前工作目錄』解開 → 把資料檔搬進 directory」。
# tar 裡附的 MANIFEST.txt 不在搬移名單內，會留在專案根目錄；而且它只保留最後一批的內容，
# 當成下載紀錄看是不完整的。順手收進資料目錄，專案根目錄就不會多出一個沒人用的孤兒檔。
if (file.exists("MANIFEST.txt")) {                       # 跨磁碟（專案在 E:、GDC_DIR 在 C:）不能用 file.rename
  if (file.copy("MANIFEST.txt", file.path(GDC_DIR, "MANIFEST.txt"), overwrite = TRUE)) file.remove("MANIFEST.txt")
}
bulk <- GDCprepare(q, directory = GDC_DIR)                                   # 本例回來 391 個檔案，含臨床欄位
bulk.mtx <- assay(bulk, "unstranded"); rownames(bulk.mtx) <- rowData(bulk)$gene_name
bulk.mtx <- bulk.mtx[!duplicated(rownames(bulk.mtx)), ]
# ★ 統計單位的問題，在這裡換到 Bulk 這一層。TCGA barcode 的第 4 段是樣本型別：
#   01 = 原發腫瘤、02 = 復發、11 = 癌旁正常組織。三種混在一起做存活分析沒有意義。
#   而且同一位病人常有兩三份 aliquot（例如 TCGA-06-0743 的 -1849-01 與 -A96S-41 是同一位），
#   不去重就等於把同一個死亡事件算兩次——跟第 34 頁「cell-level DE 把細胞當樣本」是同一個錯，
#   只是這裡重複的不是細胞而是定序檔案。
bc <- colnames(bulk.mtx); part <- substr(bc, 1, 12); styp <- substr(bc, 14, 15)
cat("\n== TCGA 樣本型別（01 原發／02 復發／11 正常）==\n"); print(table(styp))
sel <- which(styp == "01"); sel <- sel[!duplicated(part[sel])]          # 只留原發，每位病人一份
cat("檔案數", length(bc), "→ 原發腫瘤", sum(styp == "01"), "→ 去重後的病人數", length(sel), "\n")
bulk.mtx <- bulk.mtx[, sel]
common <- intersect(rownames(bulk.mtx), rownames(ref))
est <- music_prop(bulk.mtx = bulk.mtx[common, ], sc.sce = ref[common, ], clusters = "celltype_l1", samples = "patient")
prop <- as.data.frame(est$Est.prop.weighted); write.csv(prop, "output/tables/10_deconv_tcga_gbm.csv")
## ---- 2. survival ---------------------------------------------------- Q3 頁 82
# 存活：免疫（巨噬）比例上下半
clin <- as.data.frame(colData(bulk))[rownames(prop), ]
clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up) / 30.4
clin$event <- as.integer(clin$vital_status == "Dead")
# 時間是 NA 的人會被 coxph 默默刪掉。刪掉誰要先看一眼：
# 如果 NA 集中在還活著的人（days_to_last_follow_up 沒填），刪掉之後就只剩死亡的人，
# KM 曲線會被系統性拉低——那是選擇性刪除，不是隨機遺漏。
cat("\n== vital_status × 存活時間是否為 NA ==\n"); print(table(clin$vital_status, is.na(clin$time)))
cat("可分析人數", sum(!is.na(clin$time)), "／ 死亡事件", sum(clin$event[!is.na(clin$time)]), "\n")
# 本例的答案很難看，但正因為難看才要印出來：54 位存活者裡有 53 位的
# days_to_last_follow_up 是 NA，於是「還活著的人」幾乎整批被刪掉，
# 剩下 229 位裡有 227 位是死亡事件（事件率 99%）。KM 曲線因此被系統性拉低。
# 這是公開資料的常態，不是這支腳本的 bug——但它必須寫進報告的限制，不能默默略過。
# 要補救就得另外抓臨床追蹤表（GDCquery_clinic 或 clinical supplement）把追蹤時間補回來；
# 補不回來時，能說的只有組間比較（兩組同樣被削），不能報絕對中位存活。
clin$mac_hi <- prop$Immune > median(prop$Immune)
fit <- survfit(Surv(time, event) ~ mac_hi, data = clin)
p <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE, xlab = "Months"); pdf("output/figs/10_km_macrophage.pdf", 7, 6); print(p); dev.off()
print(survdiff(Surv(time, event) ~ mac_hi, data = clin))                    # log-rank 的 p 要印出來，不能只看圖上那個字
# 年齡是這裡的陽性對照：GBM 的年齡效應是已知的，它若沒出來，代表臨床欄位或時間軸接錯了。
print(summary(coxph(Surv(time, event) ~ mac_hi + age_at_index, data = clin)))   # 實際分析還要調整 MGMT、IDH
# 本例實測：391 個檔案 → 372 份原發 → 去重後 284 位病人 → 229 位可分析（227 個死亡事件）。
#   巨噬比例   HR 1.10（95% CI 0.85–1.44），log-rank p = 0.5，Cox p = 0.46 → 看不出關聯
#   年齡       HR 1.03/歲（1.02–1.04），p = 1.8e-06                       → 陽性對照有出來
# 陽性對照有出來這件事很重要：它證明臨床欄位、時間軸、模型都接對了，
# 所以「巨噬比例看不出關聯」是一個可以報的結果，不是「程式壞了」。
# 為什麼會是 null？最可能的原因是解析度：這裡的 Immune 把所有免疫細胞併成一格，
# 而文獻連結到預後的是特定的 TAM 狀態，不是巨噬細胞總量（見練習 10-3）。
sessionInfo()

# =====================================================================
# ▶ 練習 10
#  10-1 用 Malignant 比例（純度）分組做 KM，跟巨噬比例的結果方向相同嗎？兩者相關係數多少？
#  10-2 在 Cox 模型加入 age 之後，mac_hi 的 HR 變化多少？這代表什麼？
#  10-3 參考組（sc 端）把 Other 拆成 Astrocyte / OPC / Neuron 重跑：Immune 的估計比例變多少？
#       反卷積對參考的敏感度告訴你什麼？
#  進階 用 BayesPrism 重做 §1，比較兩種反卷積估的巨噬比例（相關係數、Bland–Altman 圖）。
# =====================================================================
