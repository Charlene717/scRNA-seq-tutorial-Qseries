# =====================================================================
# 05_infercnv.R — 練習腳本 5：inferCNV 惡性判定、CNV 分數與相關、三角驗證、與作者標籤比對
#
# 對應影片：Q3 頁 20–26（§1 輸入與執行、§2 兩個數字、§3 三角驗證）
# 輸入：output/rds/04_gbm4_unintegrated.rds（04_multipatient.R；用「未整合」那份）
# 輸出：output/05_infercnv/（inferCNV 原生輸出）、output/rds/05_gbm4_malignant.rds
# 時間：inferCNV 約 10–30 分鐘（denoise、無 HMM；視機器而定）
# 注意：inferCNV 底層的 rjags 需要「系統層級」的 JAGS 程式（不是 R 套件，R 裝不了它），
#       必須先在作業系統安裝 JAGS 4.x 再重開 R；未安裝的話本腳本會在 §1 直接停下並提示。
#       還沒裝 JAGS 前，06–08 可先用 celltype_author 的 Neoplastic 當替代惡性標籤測試（見 06 §0）。
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2)
set.seed(1234)
gbm4 <- readRDS("output/rds/04_gbm4_unintegrated.rds")

## ---- 1. run-infercnv ----------------------------------------------- Q3 頁 20–22
# 前置檢查：inferCNV 依賴 rjags，而 rjags 需要「系統層級」安裝 JAGS 4.x（不是 R 套件）
#   Windows：到 https://sourceforge.net/projects/mcmc-jags/files/ 下載 JAGS-4.x.y.exe 安裝
#   macOS  ：brew install jags      Linux：sudo apt install jags
#   安裝完「重新啟動 R / RStudio」，再跑本腳本
if (inherits(try(suppressPackageStartupMessages(library(rjags)), silent = TRUE), "try-error")) {
  stop("找不到 JAGS：請先在系統安裝 JAGS 4.x（見上方註解），重開 R 後再執行 05_infercnv.R。",
       call. = FALSE)
}
library(infercnv)
# (a) 註釋檔：參考組 = 確定正常的細胞（免疫、寡樹突）；觀察組按病人分開
refs <- c("Immune cell", "Oligodendrocyte")
stopifnot(all(refs %in% gbm4$celltype_author))
ann <- data.frame(row.names = colnames(gbm4),
                  group = ifelse(gbm4$celltype_author %in% refs,
                                 gbm4$celltype_author,
                                 paste0("obs_", gbm4$patient)))
write.table(ann, "data/infercnv_annot.txt", sep = "\t", col.names = FALSE, quote = FALSE)
table(ann$group)

# (b) 基因座標檔：00_setup.R 已下載 hg38_gencode_v27.txt（gene  chr  start  end）
gpos <- "data/hg38_gencode_v27.txt"
if (!file.exists(gpos)) {
  # 替代作法：用 EnsDb 自行產生（BiocManager::install("EnsDb.Hsapiens.v86")）
  library(EnsDb.Hsapiens.v86)
  g <- genes(EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name", "gene_seq_start", "gene_seq_end"))
  g <- as.data.frame(g); g <- g[g$seq_name %in% c(1:22, "X", "Y"), ]
  g <- g[!duplicated(g$gene_name), c("gene_name", "seq_name", "start", "end")]
  g$seq_name <- paste0("chr", g$seq_name)
  write.table(g, gpos, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# (c) 原始 counts（不是 data、不是 integrated）
cts <- LayerData(gbm4, assay = "RNA", layer = "counts")
obj <- CreateInfercnvObject(raw_counts_matrix = cts,
                            annotations_file  = "data/infercnv_annot.txt",
                            gene_order_file   = gpos,
                            ref_group_names   = refs)
# ⚠ inferCNV 會把每個步驟的中間結果存進 out_dir，下次跑同一個 out_dir 時「直接從最後一步接續」，
#   log 會出現 "Checking for saved results / Using backup from step 22"。中途失敗要續跑時很方便，
#   但只要你改了 cutoff、refs、denoise 這類參數，它「不會」重算——會安靜地把舊結果再給你一次。
#   下面這道檢查專門擋這件事：拿輸入的 RDS 跟快取裡「最早」的步驟檔比時間。
#   要比最早的那個：run.final 每跑一次都會重寫，看它會被騙過去。
out.dir <- "output/05_infercnv"
in.rds  <- "output/rds/04_gbm4_unintegrated.rds"
# STALE.OK 是「臨時通行證」，不是設定：只在你確定輸入內容其實沒變（例如上游只改了註解、
# 只多存了幾個 csv）時暫時改成 TRUE，跑完請立刻改回 FALSE。
# 一直放著 TRUE 等於這道檢查不存在——那正是它要防的事。
STALE.OK <- FALSE
steps <- list.files(out.dir, pattern = "^[0-9]{2}_.*\\.infercnv_obj$", full.names = TRUE)
stale <- length(steps) && min(file.mtime(steps)) < file.mtime(in.rds)
if (stale && !STALE.OK) {
  stop("inferCNV 快取比輸入舊：", out.dir, " 裡的中間結果是用更早版本的 ", basename(in.rds), " 算出來的。\n",
       "  真的換了輸入或參數 → 刪掉整個 ", out.dir, "（或換一個 out_dir）再跑。\n",
       "  確定輸入內容沒變 → 把 STALE.OK 暫時設成 TRUE 跳過這道檢查，跑完改回 FALSE。", call. = FALSE)
}
if (STALE.OK) warning("STALE.OK = TRUE：這一輪跳過了 inferCNV 快取過期檢查。\n",
                      "  跑完記得把它改回 FALSE，否則下次真的換了參數，你會拿到舊結果而不自知。",
                      call. = FALSE, immediate. = TRUE)
obj <- infercnv::run(obj,
                     cutoff = 1,                       # Smart-seq2 用 1；10x 用 0.1
                     out_dir = out.dir,
                     cluster_by_groups = TRUE,         # 每位病人各自聚類
                     denoise = TRUE,
                     HMM = FALSE,                      # 亞株分析時再開（慢很多）
                     num_threads = 4)
# 輸出的 infercnv.png 就是熱圖：上半參考（平）、下半四位病人；看 chr7 gain / chr10 loss。

## ---- 2. cnv-score-cor ---------------------------------------------- Q3 頁 23–24
# CNV 矩陣直接從 run() 回傳的物件拿，不要去讀 infercnv.observations.txt：
#   那兩個文字檔只有在 plot_cnv(write_expr_matrix = TRUE) 時才會寫出來，run() 預設不寫，
#   讀了會得到「無法開啟連接」。物件裡的 expr.data 就是同一份資料，而且省掉幾百 MB 的文字讀寫。
# 若 R 重開過、obj 已經不在 session 裡，從 run() 存下的最終物件讀回來即可（不必重跑 inferCNV）：
if (!exists("obj")) obj <- readRDS("output/05_infercnv/run.final.infercnv_obj")
cnv.all <- obj@expr.data                              # 基因 × 細胞；已平滑、已 denoise，以參考為中心（≈1）
stopifnot("inferCNV 的細胞與 gbm4 對不上：確認讀的是同一份 04 物件" =
            setequal(colnames(cnv.all), colnames(gbm4)))

cnv.score <- colMeans((cnv.all - 1)^2)                # ① 偏離參考的程度
top       <- names(sort(cnv.score, decreasing = TRUE))[1:200]
mal.prof  <- rowMeans(cnv.all[, top])                 #   「最像惡性」的 200 顆平均 profile
cnv.cor   <- apply(cnv.all, 2, cor, y = mal.prof)     # ② 與惡性 profile 的相關（Tirosh 2016）

gbm4$cnv.score <- cnv.score[colnames(gbm4)]
gbm4$cnv.cor   <- cnv.cor[colnames(gbm4)]
# 閾值由「參考組自己的分布」長出來，不要寫死數字：參考組是確定正常的細胞，
# 它們的分位數就是「正常細胞最多長成什麼樣」，超過才有資格談惡性。
# 換一份資料、換一組參考組，兩條線會自己跟著動，不用回來改數字。
ref.cell <- gbm4$celltype_author %in% refs
s.hi <- quantile(gbm4$cnv.score[ref.cell], 0.99)   # 正常的 CNV 分數上限
c.hi <- quantile(gbm4$cnv.cor[ref.cell],   0.99)   # 正常的相關上限 → 超過才談惡性
c.lo <- quantile(gbm4$cnv.cor[ref.cell],   0.90)   # 低於這條才算「明確正常」
cat(sprintf("參考組 %d 顆決定的閾值：cnv.score > %.4f、cnv.cor > %.3f（明確正常：cnv.cor < %.3f）\n",
            sum(ref.cell), s.hi, c.hi, c.lo))

x.max <- quantile(gbm4$cnv.score, 0.995)           # 少數離群點會把 x 軸壓扁；只縮放視野，不刪點
p <- ggplot(gbm4@meta.data, aes(cnv.score, cnv.cor, colour = celltype_author)) +
     geom_point(size = .7, alpha = .7) +
     geom_vline(xintercept = s.hi, linetype = 2, colour = "grey35") +
     geom_hline(yintercept = c.hi, linetype = 2, colour = "grey35") +
     geom_hline(yintercept = c.lo, linetype = 3, colour = "grey60") +
     annotate("text", x = x.max, y = c.hi, hjust = 1, vjust = -0.6, size = 3.2, colour = "grey35",
              label = sprintf("cnv.cor > %.3f", c.hi)) +
     annotate("text", x = x.max, y = c.lo, hjust = 1, vjust = 1.5, size = 3.2, colour = "grey55",
              label = sprintf("cnv.cor < %.3f = 明確正常", c.lo)) +
     annotate("text", x = s.hi, y = max(gbm4$cnv.cor), hjust = -0.06, vjust = 1, size = 3.2, colour = "grey35",
              label = sprintf("cnv.score > %.4f", s.hi)) +
     coord_cartesian(xlim = c(0, x.max)) + theme_classic() +
     labs(x = "CNV score (deviation from reference)", y = "CNV correlation (with malignant profile)")
ggsave("output/figs/05_cnv_scatter.png", p, width = 8, height = 6, dpi = 150, bg = "white")
# 看圖：右上（兩條虛線之外）= 惡性；點線以下 = 正常；中間那條帶 = 不確定。
# 分位數只是起點：如果你的圖上兩群之間有明顯的谷，把線移到谷底會比分位數更好。

## ---- 3. triangulate ------------------------------------------------ Q3 頁 25–26
# 三條閾值在 §2 已由參考組算出（s.hi / c.hi / c.lo）。覺得不合就回 §2 調，這裡只做判定。
gbm4$malignant <- with(gbm4@meta.data, ifelse(
  cnv.score > s.hi & cnv.cor > c.hi, "malignant",
  ifelse(cnv.score <= s.hi & cnv.cor < c.lo, "normal", "unresolved")))

# 第三項證據：譜系。這幾種細胞在 GBM 裡不可能是腫瘤細胞——免疫與血管來自別的胚層，
# 寡樹突與神經元是終末分化。它們被標成惡性一定是 doublet 或雜訊，一律改 unresolved。
# ★ OPC 與 Astocyte 刻意不放進來：OPC-like / AC-like 正是惡性狀態的名字，
#   把它們擋掉等於先射箭再畫靶。那兩群本來就該留在 unresolved 讓證據說話。
lineage.normal <- c("Immune cell", "Vascular", "Oligodendrocyte", "Neuron")
gbm4$malignant[gbm4$malignant == "malignant" & gbm4$celltype_author %in% lineage.normal] <- "unresolved"

# 第二項證據：按病人成島。惡性細胞應集中在 raw_clusters 中「單一病人為主」的群
tab <- table(gbm4$malignant, gbm4$celltype_author)
print(tab)
# 與作者標籤比對。整體一致率會被大量的「正常細胞判對」灌高，所以三個數字要一起看：
mal  <- gbm4$malignant == "malignant"; neo <- gbm4$celltype_author == "Neoplastic"
cat(sprintf("整體一致率 %.1f%%｜精確率 %.1f%%（判為惡性的裡面有幾成真的是）｜召回率 %.1f%%（作者的惡性細胞抓回幾成）\n",
            100 * mean(mal == neo), 100 * sum(mal & neo) / sum(mal), 100 * sum(mal & neo) / sum(neo)))
print(table(gbm4$malignant, gbm4$tissue))
# 自我檢查：參考組是確定正常的細胞，落進 unresolved 的比例應該接近 c.lo 的分位數（這裡是 10%）。
#   遠超過的話，代表參考組裡混了不該混的東西（例如把腫瘤旁的反應性細胞也當成正常）。
cat(sprintf("參考組落進 unresolved 的比例：%.1f%%（c.lo 取 90 分位 → 預期約 10%%）\n",
            100 * mean(gbm4$malignant[ref.cell] == "unresolved")))
## >>> 參考答案 ------------------------------------------------------
# 本課的結果（實跑驗證；refs = 免疫 + 寡樹突，閾值由參考組分位數自動決定）：
#   參考組 1,915 顆決定出：cnv.score > 0.0137、cnv.cor > 0.210（明確正常：cnv.cor < 0.115）
#   malignant 633 顆、normal 2,051 顆、unresolved 855 顆（合計 3,539）
#   整體一致率 87.6%｜精確率 99.4%（629/633）｜召回率 59.1%（629/1,064）
#   參考組落進 unresolved 的比例 10.2%（預期約 10%，見下面 ①）
#
# 整體一致率會騙人：絕大部分是「正常細胞判成不是惡性」灌上去的。要看的是後面兩個數字。
#   精確率 99.4% → 判為惡性的 633 顆裡，只有 4 顆不是作者標的 Neoplastic
#   召回率 59.1% → 作者標為 Neoplastic 的細胞，抓回約六成
# 這個取捨是刻意的，不是失敗：這批細胞接下來要拿去做差異表達，
#   混進一顆正常細胞的代價，遠大於少算一顆惡性細胞。寧可漏，不可錯。
#
# 對照組：把閾值寫死成 c.hi = 0.4、c.lo = 0.2（本課早期版本）會得到
#   malignant 508、normal 2,418、unresolved 613；一致率 84.2%、精確率 99.8%、召回率 47.7%。
#   換成由參考組決定之後，召回率 47.7% → 59.1%（多抓回 122 顆），精確率 99.8% → 99.4%。
#   這就是「讓資料決定閾值」比「抄一個數字」好的地方——而且換一份資料不用回來改。
# 譜系否決的效果：把 Neuron 加進 lineage.normal 之後，4 顆被誤判的神經元從 malignant 移到
#   unresolved，精確率 98.7% → 99.4%，召回率完全不動（神經元本來就不是 Neoplastic）。
#   這說明第三項證據的角色是「否決」，不是「加分」——它只會把錯的拿掉，不會多抓對的回來。
#
# 兩個要看懂的副作用：
#   ① 有 195 顆參考組細胞（免疫 187 + 寡樹突 8）落進 unresolved，佔參考組的 10.2%。
#      這是 c.lo 取 90 分位的必然結果，不是錯誤——定義上就會有一成的參考細胞在線上面。
#      它同時是個好用的自我檢查：如果遠超過一成，代表你的參考組裡混了不該混的東西。
#   ② OPC 有 213 顆 unresolved、187 顆 normal。OPC 正是 OPC-like 惡性狀態的正常對應細胞，
#      本來就最難分——這裡不硬判，正是誠實的做法。
#
# 看 05_cnv_scatter.png：這兩個數字不是平分秋色的。
#   縱軸 cnv.cor 把 Neoplastic 和正常細胞分得很開；橫軸 cnv.score 兩者重疊很多，
#   還被一顆 0.105 的離群點拉扁（所以圖上用 coord_cartesian 縮到 99.5 分位）。
#   在這份資料裡真正在做事的是「與惡性 profile 的相關」，不是「偏離參考的程度」。
#   換一份資料不一定是同一個數字在做事——所以兩個都算、都畫出來，再決定閾值。
## <<< 參考答案

p <- DimPlot(gbm4, reduction = "umap.raw", group.by = "malignant", cols = c("#C0392B", "#1A6B5A", "grey70"))
ggsave("output/figs/05_umap_malignant.png", p, width = 7, height = 6, dpi = 150, bg = "white")
saveRDS(gbm4, "output/rds/05_gbm4_malignant.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 5
#  5-1 【錯誤示範】把 refs 改成 c("Astocyte")（注意 GEO 原檔就是少一個 r），重跑 inferCNV（可只跑 §1）。
#      ★ 一定要同時把 out_dir 換成另一個資料夾，例如 "output/05_infercnv_badref"——
#        沿用原資料夾的話 inferCNV 會直接載入舊結果，你會以為「換了參考組結果沒變」。熱圖變成什麼樣？
#      chr7 / chr10 的訊號還在嗎？用一句話解釋為什麼。
#  5-2 把 top 從 200 改成 50 與 500，cnv.cor 的分布變了多少？判定結果（malignant 數）差幾顆？
#  5-3 只看 Periphery 的細胞：作者標 Neoplastic、但你標 unresolved 的有幾顆？
#      它們的 cnv.score 分布跟 Tumor 的 Neoplastic 比如何？這告訴你浸潤細胞的什麼特性？
#  5-4 對每位病人，畫出該病人惡性細胞在 chr7 與 chr10 的平均 CNV 值（cnv.all 的列名含基因，
#      配合 gpos 找出染色體）。四位病人都有 chr7+/chr10− 嗎？
#  進階 開 HMM = TRUE 重跑其中一位病人，看 infercnv 的 subcluster 結果：這位病人有幾個亞株？
#      各亞株的私有事件是什麼？
# =====================================================================
