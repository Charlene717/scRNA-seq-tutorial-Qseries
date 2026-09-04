# =====================================================================
# 05_infercnv.R — 練習腳本 5：inferCNV 惡性判定、CNV 分數與相關、三角驗證、與作者標籤比對
#
# 對應影片：Q3 頁 20–26（§1 輸入與執行、§2 兩個數字、§3 三角驗證）
# 輸入：output/rds/04_gbm4_unintegrated.rds（04_multipatient.R；用「未整合」那份）
# 輸出：output/rds/05_infercnv/（inferCNV 原生輸出）、output/rds/05_gbm4_malignant.rds、output/figs/05_infercnv*.png
# 時間：inferCNV 約 10–30 分鐘（denoise、無 HMM；視機器而定）
# 注意：inferCNV 底層的 rjags 需要「系統層級」的 JAGS 程式（不是 R 套件，R 裝不了它），
#       必須先在作業系統安裝 JAGS 4.x 再重開 R；未安裝的話本腳本會在 §1 直接停下並提示。
#       還沒裝 JAGS 前，06–08 可先用 celltype_author 的 Neoplastic 當替代惡性標籤測試（見 06 §0）。
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
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
## TODO ▶ 誰是「確定正常」的參考組？為什麼不是星狀細胞？（Q3 頁 20、22）
refs <- c("____", "____")
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
#   下面這道檢查專門擋這件事。
#   比的是「輸入內容 + 參數」的指紋，不是檔案時間：把 04 重跑一次，RDS 的時間就變新了，
#   內容卻一模一樣；用時間比會在這種情況誤判，逼你白白重跑半小時。指紋只有真的換了東西才會不同。
out.dir <- "output/rds/05_infercnv"       # inferCNV 的原生輸出整包放在 rds/ 底下，跟其他分析物件同一層

# 參數先寫成變數，run() 再引用同一組——指紋跟真正跑的參數綁在一起，改了一定會被抓到。
# 寫成兩份（一份給 run()、一份給指紋）遲早會不同步，那時檢查就形同虛設。
## TODO ▶ Smart-seq2 與 10x 的 cutoff 差 10 倍，這份資料該用哪個？（Q3 頁 22）
CUTOFF        <- ____
CLUSTER.GROUP <- TRUE     # 每位病人各自聚類
DENOISE       <- TRUE
USE.HMM       <- FALSE    # 亞株分析時再開（慢很多）
sig <- function(v) {                                     # 內容指紋：有 digest 就用它，沒有退成便宜的摘要
  v <- as.character(v)
  if (requireNamespace("digest", quietly = TRUE)) digest::digest(v)
  else paste(length(v), sum(nchar(v)), v[1], v[length(v)], sep = "|")
}
run.key <- list(dim     = dim(cts), total = sum(cts),
                cells   = sig(colnames(cts)), genes = sig(rownames(cts)),
                groups  = sig(ann$group),     refs  = sort(refs),
                genepos = sig(readLines(gpos, warn = FALSE)),
                params  = list(cutoff = CUTOFF, cluster_by_groups = CLUSTER.GROUP,
                               denoise = DENOISE, HMM = USE.HMM))
stamp    <- file.path(out.dir, "_run_stamp.rds")         # §1b 跑完後寫：這份快取是用哪一組指紋算出來的
cached   <- list.files(out.dir, pattern = "\\.infercnv_obj$")
prev.key <- if (file.exists(stamp)) readRDS(stamp)$key else NULL

if (length(cached) && is.null(prev.key)) {
  # 舊版腳本留下的快取沒有指紋。不要用猜的：01_incoming_data 裡就存著當初送進去的矩陣，直接比對。
  inc  <- file.path(out.dir, "01_incoming_data.infercnv_obj")
  old  <- if (file.exists(inc)) try(readRDS(inc)@count.data, silent = TRUE) else NULL
  same <- !is.null(old) && !inherits(old, "try-error") &&
          identical(colnames(old), colnames(cts)) && all(rownames(old) %in% rownames(cts)) &&
          isTRUE(all.equal(sum(old), sum(cts[rownames(old), ])))
  if (!same)
    stop("out_dir 裡有舊版腳本留下的快取，而且它的輸入跟現在這一份對不上。\n",
         "  要重算：unlink(\"", out.dir, "\", recursive = TRUE) 或換一個 out_dir，再跑一次。",
         call. = FALSE)
  message("舊快取的輸入與現在這一份一致，沿用並補寫指紋。\n",
          "  注意：舊快取沒有留下當初的 cutoff / refs / denoise，這三項無法回頭確認；\n",
          "  若你在那次之後改過它們，請刪掉 out_dir 重跑。")
} else if (length(cached) && !identical(prev.key, run.key)) {
  chg <- names(run.key)[!vapply(names(run.key),
                                function(k) identical(prev.key[[k]], run.key[[k]]), logical(1))]
  stop("inferCNV 快取跟現在的設定對不上，變動的是：", paste(chg, collapse = "、"), "。\n",
       "  inferCNV 不會因為參數變了就重算，所以先在這裡停下來。\n",
       "  要重算：unlink(\"", out.dir, "\", recursive = TRUE) 或換一個 out_dir，再跑一次。",
       call. = FALSE)
}
obj <- infercnv::run(obj,
                     cutoff = CUTOFF,                  # 與上面指紋用的是同一組參數
                     out_dir = out.dir,
                     cluster_by_groups = CLUSTER.GROUP,
                     denoise = DENOISE,
                     HMM = USE.HMM,
                     num_threads = 4)
# 輸出的 infercnv.png 就是熱圖：上半參考（平）、下半四位病人；看 chr7 gain / chr10 loss。
# 熱圖是 inferCNV 自己寫在 out_dir 裡的，檔名沒有腳本編號。複製一份到 figs/ 並補上 05_ 前綴，
# 交付時所有的圖就都在同一個資料夾，不用再去翻原生輸出。
for (f in c("infercnv.png", "infercnv.preliminary.png", "infercnv_subclusters.png")) {
  if (file.exists(file.path(out.dir, f)))
    file.copy(file.path(out.dir, f), file.path("output/figs", sub("^infercnv", "05_infercnv", f)), overwrite = TRUE)
}

## ---- 1b. 中間結果要不要留（磁碟 vs 續跑）----------------------------
# 這一步不影響任何分析結果，只影響磁碟。先把事實列出來，再自己決定。
#
# inferCNV 每做完一步就把整個物件存一份到 out_dir。本例（4 位病人、3,589 顆細胞）：
#   13 個中間步驟檔（01_incoming … 22_denoise、preliminary）   約 2.9 GB
#   run.final.infercnv_obj（§2 之後唯一會讀的）                 約 95 MB
# 細胞數越多差距越大——換成 10x 的幾萬顆，中間檔會是好幾十 GB。
#
# 留著的好處，只有一個：**中途失敗可以續跑**。inferCNV 下次會從最後一個成功的步驟接上去
# （log 會出現 "Using backup from step 22"），不必從頭再跑十幾分鐘。
# 但這個好處只在「同一組參數、跑到一半掛掉」時有用。順利跑完之後就用不到了：
# 改了 cutoff / refs / denoise 本來就該重算，不該用舊快取。
#
# 清掉的代價，要看清楚：
#   ① 之後任何一次中斷都要從第一步重跑（本例約 10–30 分鐘）。
#   ② 上面那道「快取跟設定對不上」的檢查靠的是下面寫的 _run_stamp.rds（幾百 bytes），
#      不是步驟檔本身。中間檔清掉，檢查一樣有效——不是把安全網拿掉。
#   ③ run.final 一定會保留，§2 之後完全不受影響；熱圖也已經複製到 figs/ 了。
#
# 建議：還在調參數、或機器容易中斷 → 保持 FALSE。教學跑完、確定不再改參數 → 改 TRUE 收回 2.9 GB。
CLEAN.INTERMEDIATE <- FALSE
saveRDS(list(key = run.key, when = Sys.time()), stamp)   # 先寫指紋，再清才安全
if (CLEAN.INTERMEDIATE) {
  junk <- setdiff(list.files(out.dir, pattern = "\\.infercnv_obj$", full.names = TRUE),
                  file.path(out.dir, "run.final.infercnv_obj"))
  if (length(junk)) {
    cat("清掉 inferCNV 中間結果：", length(junk), "個檔案、",
        round(sum(file.size(junk)) / 1024^3, 2), "GB；run.final.infercnv_obj 保留\n")
    file.remove(junk)
  }
}

## ---- 2. cnv-score-cor ---------------------------------------------- Q3 頁 23–24
# CNV 矩陣直接從 run() 回傳的物件拿，不要去讀 infercnv.observations.txt：
#   那兩個文字檔只有在 plot_cnv(write_expr_matrix = TRUE) 時才會寫出來，run() 預設不寫，
#   讀了會得到「無法開啟連接」。物件裡的 expr.data 就是同一份資料，而且省掉幾百 MB 的文字讀寫。
# 若 R 重開過、obj 已經不在 session 裡，從 run() 存下的最終物件讀回來即可（不必重跑 inferCNV）：
if (!exists("obj")) obj <- readRDS(file.path(out.dir, "run.final.infercnv_obj"))
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

# 第二項證據：按病人各成一個群集。惡性細胞應集中在 raw_clusters 中「單一病人為主」的群
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
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）

p <- DimPlot(gbm4, reduction = "umap.raw", group.by = "malignant", cols = c("#C0392B", "#1A6B5A", "grey70"))
ggsave("output/figs/05_umap_malignant.png", p, width = 7, height = 6, dpi = 150, bg = "white")
saveRDS(gbm4, "output/rds/05_gbm4_malignant.rds")
sessionInfo()

# =====================================================================
# ▶ 練習 5
#  5-1 【錯誤示範】把 refs 改成 c("Astocyte")（注意 GEO 原檔就是少一個 r），重跑 inferCNV（可只跑 §1）。
#      ★ 一定要同時把 out_dir 換成另一個資料夾，例如 "output/rds/05_infercnv_badref"——
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
