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
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
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
# 重複的基因 symbol：多個 Ensembl 基因（ENSG）對到同一個名字。
# 最省事的做法是 !duplicated() 留第一個，但列的順序是照 Ensembl ID 排的、跟生物學無關，
# 等於隨機挑一個留下，也寫不出一條像樣的方法段。下面 §1b 改成三步處理。
# ⚠ 不要反射性地用 rowsum() 加總。加總對「轉錄本 → 基因」是對的（那些本來就是同一個基因的片段），
#   但這裡的重複是「不同的基因座剛好共用一個名字」，性質完全不同：
#     · _PAR_Y：擬體染色體區的基因在 X 與 Y 各註記一次。這確實是同一個基因，加總無妨
#       （Y 那一份通常是 0，因為讀序都被指派到 X）。
#     · 其餘的撞名：兩個不同的基因座剛好共用一個 HGNC 符號。這一類**不是**同一個基因，
#       加起來就是把兩個不同的東西混成一個，而且事後從數字上看不出來。
#   注意這裡的矩陣已經是**基因層級**（一列 = 一個 ENSG），轉錄本早在上游就併進母基因了，
#   所以撞名跟 isoform 無關——這也是為什麼「加總」在這一層不適用。
#   至於這份檔案裡除了 _PAR_Y 之外還有哪些撞名、各是什麼類型，下面 §1b 的第①步會印出來給你看。
#   （不要憑印象假設，先看過再決定怎麼處理——練習 10-4 就是在練這件事。）
#
# 那要怎麼辦？先講一件常被誤導的事：**「全程用 Ensembl ID 就沒事了」並不成立。**
#   轉換只是被延後，不會消失——只要你要跟「用符號當索引的外部資源」對接，就得轉一次：
#     · 本節的 MuSiC：參考組 ref 是單細胞物件，rownames 就是基因符號（見下面的 intersect），
#       所以 bulk 這一側非得在符號空間不可，否則對不起來。
#     · 06a 的富集：msigdbr 給的基因集是 gene_symbol；KEGG 那一路還要 Entrez。
#     · marker 清單、CellChatDB 也都是符號。
#   真正能控制的不是「要不要轉」，而是**在哪裡轉、轉的時候有沒有看見**。
#
# 所以本節留在符號空間，但把折疊規則寫明，分三步做（下面 §1b）：
#   ① 先看重複的組成——不知道重複是什麼，就沒有立場決定怎麼處理
#   ② 用 gene_type 拆撞名：判斷依據是「註解」不是「這批資料的數值」，換一批病人結果一樣
#      ⚠ 別預期它能解決大部分。本例 1,233 個撞名列裡，這一步只拆掉 7 列。
#         看第①步印出來的組成就知道為什麼：撞名絕大多數發生在 misc_RNA（937）、snoRNA（176）、
#         snRNA（69）彼此之間，同一個符號底下根本沒有蛋白編碼版本可以留，所以這一步不動它們。
#         這正是「先看組成再決定怎麼處理」的價值——先看，你才知道第②步在這份資料上是小角色，
#         真正在做事的是第③步，而第③步是依賴資料的選擇（疑慮見下）。
#   ③ 殘餘的留表現量最高的那一列，並把這是一個「選擇」講清楚
#
# ⚠ 第③步要保留的疑慮，教的時候要一起講：
#   · 這是依賴資料的選擇。留哪一列取決於這個 cohort 的表現量，換一批病人可能留到不同的基因座；
#     要把兩個 cohort 併起來時，同一個符號兩邊可能不是同一個東西。
#   · 如果殘餘的重複是「兩個真的不同的基因座剛好共用一個名字」，留一個等於把另一個
#     默默從分析裡刪掉。挑豐度高的是一個決定，不是一個事實。
#   · 更保守的做法是把這些語意不明的符號整個排除（反卷積用幾千個基因，少一小撮不痛不癢，
#     而默默留錯一個事後看不出來）。本課選擇保留、但要求你把數量與規則寫進方法段——
#     練習 10-4 會請你兩種都跑一次，看反卷積比例差多少。
#
## ---- 1b. 重複基因符號的三步處理 --------------------------------------
rd <- rowData(bulk)      # 與上面第 39 行同一個來源，列順序仍然一一對應
cat("\n== 重複的基因符號 ==\n")
cat("重複列數：", sum(duplicated(rownames(bulk.mtx))), " / 全部 ", nrow(bulk.mtx), "\n", sep = "")

# ① 看組成：這些重複是什麼？（欄位不一定存在，先確認再用，不要假設）
dup.name <- rownames(bulk.mtx) %in% rownames(bulk.mtx)[duplicated(rownames(bulk.mtx))]
if ("gene_id" %in% names(rd)) {
  par.y <- grepl("_PAR_Y$", rd$gene_id)
  cat("其中 _PAR_Y（擬體染色體區，X 與 Y 各記一次，是同一個基因）：", sum(par.y & dup.name), "\n", sep = "")
}
if ("gene_type" %in% names(rd)) {
  cat("重複列的 gene_type 組成：\n"); print(sort(table(rd$gene_type[dup.name]), decreasing = TRUE))
} else {
  cat("（這份 rowData 沒有 gene_type 欄位，跳過第②步的註解篩選）\n")
}

# ② 用註解拆撞名——只動撞名的那幾列，不要整批篩。
#    同一個符號底下如果同時有蛋白編碼與非蛋白編碼，非蛋白編碼的那幾列拿掉即可。
#    這一步不看數值，換一批病人結果一樣——這正是它比「挑表現量高」可靠的地方。
#    ★ 不要圖方便寫成 bulk.mtx[rd$gene_type == "protein_coding", ]：那會把所有沒撞名的
#      lncRNA 與假基因也一起刪掉，等於偷偷換掉整個基因池，不再只是「解決撞名」。
#      （若你確實想只用蛋白編碼做反卷積，那是另一個要獨立說明的決定，不要藏在這一步裡。）
if ("gene_type" %in% names(rd)) {
  pc        <- rd$gene_type == "protein_coding"
  has.pc    <- rownames(bulk.mtx) %in% rownames(bulk.mtx)[dup.name & pc]   # 這個符號有蛋白編碼版本
  drop.amb  <- dup.name & !pc & has.pc                                     # 撞名、非蛋白編碼、且有蛋白編碼可留
  cat("依註解拆掉的撞名列（非蛋白編碼、且同名有蛋白編碼版本）：", sum(drop.amb), "\n", sep = "")
  bulk.mtx <- bulk.mtx[!drop.amb, ]
  rd       <- rd[!drop.amb, ]
  cat("剩下的重複：", sum(duplicated(rownames(bulk.mtx))), " 列\n", sep = "")
}

# ③ 殘餘的：留表現量最高的那一列（不是留第一個——列的順序是照 Ensembl ID 排的，跟生物學無關）
n.before <- nrow(bulk.mtx)
ord     <- order(rowSums(bulk.mtx), decreasing = TRUE)
sym.ord <- rownames(bulk.mtx)[ord]
collapsed.sym <- unique(sym.ord[duplicated(sym.ord)])   # 有列被折疊掉的那些符號，等一下要量它們的影響
bulk.mtx <- bulk.mtx[ord, ][!duplicated(sym.ord), ]
cat("殘餘重複依表現量折疊：", n.before, " → ", nrow(bulk.mtx), " 列",
    "（丟掉 ", n.before - nrow(bulk.mtx), " 列，方法段要寫這個數字與規則）\n", sep = "")
# ★ 統計單位的問題，在這裡換到 Bulk 這一層。TCGA barcode 的第 4 段是樣本型別：
#   01 = 原發腫瘤、02 = 復發、11 = 癌旁正常組織。三種混在一起做存活分析沒有意義。
#   而且同一位病人常有兩三份 aliquot（例如 TCGA-06-0743 的 -1849-01 與 -A96S-41 是同一位），
#   不去重就等於把同一個死亡事件算兩次——跟第 32 頁「cell-level DE 把細胞當樣本」是同一個錯，
#   只是這裡重複的不是細胞而是定序檔案。
bc <- colnames(bulk.mtx); part <- substr(bc, 1, 12); styp <- substr(bc, 14, 15)
cat("\n== TCGA 樣本型別（01 原發／02 復發／11 正常）==\n"); print(table(styp))
sel <- which(styp == "01"); sel <- sel[!duplicated(part[sel])]          # 只留原發，每位病人一份
cat("檔案數", length(bc), "→ 原發腫瘤", sum(styp == "01"), "→ 去重後的病人數", length(sel), "\n")
bulk.mtx <- bulk.mtx[, sel]
# 折疊規則到底影響了多少？不要用猜的——量出來。
# 如果被折疊過的符號根本不在共同基因裡，那第③步那個「依賴資料的選擇」對本節結論沒有影響，
# 你可以放心地在方法段寫一句話帶過；反過來如果佔比不小，就該認真考慮
# 「整批排除語意不明的符號」那條更保守的路（練習 10-4）。這個檢查決定了你該用哪一種寫法。
common <- intersect(rownames(bulk.mtx), rownames(ref))
cat("\n== 折疊規則的影響範圍 ==\n")
cat("Bulk 與單細胞參考的共同基因：", length(common), "\n", sep = "")
cat("被折疊過的符號 ", length(collapsed.sym), " 個，其中進到共同基因的：",
    sum(collapsed.sym %in% common), " 個\n", sep = "")
## TODO ▶ 反卷積的參考用哪一層型別、哪一欄當樣本？（Q3 頁 82）
est <- music_prop(bulk.mtx = bulk.mtx[common, ], sc.sce = ref[common, ], clusters = "____", samples = "____")
prop <- as.data.frame(est$Est.prop.weighted); write.csv(prop, "output/tables/10_deconv_tcga_gbm.csv")
## ---- 2. survival ---------------------------------------------------- Q3 頁 82
# 存活：免疫細胞「總」比例的上下半
# ⚠ 命名要跟算出來的東西一致。這裡的 prop$Immune 是所有免疫細胞合起來的比例，
#   不是巨噬細胞比例——參考組在 §1 只分到 Immune 這一層，拆不出 TAM。
#   把變數叫 mac_、圖檔叫 macrophage，讀的人會以為看到的是巨噬細胞，那是過度解讀。
clin <- as.data.frame(colData(bulk))[rownames(prop), ]
clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up) / 30.4
clin$event <- as.integer(clin$vital_status == "Dead")
# vital_status 有第三類 Not Reported（本例 1 位）。上面這兩行等於把它當成「存活、設限」處理，
# 這是可接受的預設，但它是一個決定不是事實——人數多的時候要單獨列出來或做敏感度分析。
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
# ⚠ 補不回來時，**不能因為「兩組都被削」就推論組間比較不受影響**。
#   刪除的比例相同，不代表 HR 無偏——關鍵在於刪除的機制取決於「有沒有死」：
#   那 53 位其實還活著的病人，從來沒有進入他們本該參與的 risk set，
#   而 Cox 的每一份資訊都來自「某個時點死掉的人 vs 當時仍在 risk set 裡的人」。
#   KM 與 Cox 都依賴 censoring 與結果無關（non-informative censoring）這個假設，
#   現在缺失幾乎全部集中在存活者身上，這個假設站不住。
#   所以絕對存活率、中位存活時間、以及組間的 HR 與 p 值，
#   在這一版資料上都只能當**探索性示範**，不能當成 population-level 的正式推論。
#   要變成正式結論，得先從 GDC 的 clinical / follow-up supplement 把 days_to_last_follow_up 補齊。
clin$imm_hi  <- prop$Immune > median(prop$Immune)      # 二分：畫 KM 用
clin$imm_pct <- prop$Immune * 100                      # 連續：Cox 用（不用先切成兩組）
fit <- survfit(Surv(time, event) ~ imm_hi, data = clin)
p <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE, xlab = "Months"); pdf("output/figs/10_km_immune.pdf", 7, 6); print(p); dev.off()
print(survdiff(Surv(time, event) ~ imm_hi, data = clin))                    # log-rank 的 p 要印出來，不能只看圖上那個字
# 年齡是這裡的陽性對照：GBM 的年齡效應是已知的，它若沒出來，代表臨床欄位或時間軸接錯了。
# 中位數切兩半會丟掉組內的變異，只適合畫圖；Cox 以連續變項為主要分析，
# HR 讀成「免疫比例每多 1 個百分點」的風險比。兩種都印出來，結論要一致才站得住。
print(summary(coxph(Surv(time, event) ~ imm_pct + age_at_index, data = clin)))  # 主要分析：連續
print(summary(coxph(Surv(time, event) ~ imm_hi  + age_at_index, data = clin)))  # 對照：二分
# 實際分析還要調整 MGMT 甲基化、IDH 狀態；免疫比例也受腫瘤純度影響，正式報告要做敏感度分析。
# 本例實測（2026-09-11，三步折疊版；R 4.4.1 / MuSiC 1.0.0 / TCGAbiolinks 2.32.0）：
#   重複符號   1,233 列 / 60,660；_PAR_Y 44；依註解拆掉 7 列；殘餘 1,226 列依表現量折疊
#              → 60,653 → 59,427 列。共同基因 MuSiC 實際用了 17,304 個、5 個型別
#   病人流程   391 個檔案 → 372 份原發 → 去重後 284 位病人 → 229 位可分析（227 個死亡事件）
#   免疫總比例（連續，主要分析）HR 1.006／百分點（95% CI 0.991–1.022），p = 0.45
#   免疫總比例（二分，對照）    HR 1.105（95% CI 0.850–1.437），p = 0.46；log-rank p = 0.5
#     → 兩種切法結論一致：在 complete-case 子集裡看不出關聯
#       （探索性，理由見上面的 censoring 說明；二分版只用來畫 KM，不當主要分析）
#   年齡       HR 1.028/歲（1.016–1.039），p = 1.6e-06                    → 陽性對照有出來
#   ⚠ 換一版折疊規則（例如練習 10-4 的「整批排除」）這些數字會微幅改變，
#     但只要「折疊規則的影響範圍」那個檢查印出來的數字小，結論方向不會因此翻轉。
# 陽性對照有出來仍然很重要，但它能說的有限度：
#   能說 → 已知的年齡效應被偵測到，支持臨床欄位、時間單位、病人對接與模型設定大致合理，
#          免疫比例的 null 也不是單純因為程式整個失效。
#   不能說 → 它**無法**排除「存活者追蹤時間大量缺失」造成的 selection bias。
#          年齡與免疫比例是不同的 predictor，同一個 selection 機制對兩者的影響不必相同；
#          年齡效應夠強所以還看得到，不代表比較弱的關聯沒有被這個機制吃掉。
# 所以正確的寫法是：在「取得到存活時間」的 complete-case 子集裡，沒有觀察到免疫比例與存活的明顯關聯；
#   由於存活者的追蹤時間幾乎全面缺失，此結果僅供探索，不能當成
#   「TCGA-GBM 中免疫比例與存活無關」的正式證據。
# 為什麼會是 null？最可能的原因是解析度：這裡的 Immune 把所有免疫細胞併成一格，
# 而文獻連結到預後的是特定的 TAM 狀態，不是巨噬細胞總量（見練習 10-3）。
sessionInfo()

# =====================================================================
# ▶ 練習 10
#  10-1 用 Malignant 比例（純度）分組做 KM，跟免疫比例的結果方向相同嗎？兩者相關係數多少？
#  10-2 在 Cox 模型加入 age 之後，imm_pct 的 HR 變化多少？這代表什麼？
#  10-3 參考組（sc 端）把 Other 拆成 Astrocyte / OPC / Neuron 重跑：Immune 的估計比例變多少？
#       反卷積對參考的敏感度告訴你什麼？
#  10-4 §1b 的三步折疊，第③步是「依賴資料的選擇」。把它換成兩種替代做法各跑一次，
#       比較免疫比例（相關係數＋最大差值）與 Cox 的 HR：
#       (a) 最偷懶的版本：!duplicated() 留第一個（列序照 Ensembl ID，等於隨機挑）。
#       (b) 最保守的版本：把殘餘撞名的符號整批排除，不留任何一列。
#       三者差多少？差很小的話，原因是什麼（提示：看 §1 印出的「折疊規則的影響範圍」）？
#       如果差很小，方法段還需要寫這一段嗎——為什麼「影響小」本身也是要報告的結果？
#       兩個延伸問題：(c) 如果改成 rowsum() 加總，哪一類重複會被加錯、為什麼？
#       (d) 有人主張「全程用 Ensembl ID 就不會有這個問題」——看一下 §1 的 intersect 那一行，
#       這個主張在本節成立嗎？要成立的話，還得多做什麼？
#  進階 用 BayesPrism 重做 §1，比較兩種反卷積估的免疫比例（相關係數、Bland–Altman 圖）。
# =====================================================================
