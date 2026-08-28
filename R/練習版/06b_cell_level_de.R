# =====================================================================
# 06b_cell_level_de.R — 練習腳本 6b：不能 pseudobulk 時的差異表達（cell-level DE 的正確姿勢）
#
# 對應影片：Q3 頁 48–49（三道防線：病人共變量、逐病人一致性、標籤置換）
# 什麼時候用這支：pseudobulk 需要「每個條件 ≥ 2–3 個生物重複」。當病人太少（06a 的 de.summary
#   顯示某型別不到 MIN_PAIRS）、或整份資料只有一位病人時，退而在細胞層級做——但要知道
#   p 值是被灌大的（假重複），所以這支的重點不是 p 值，而是三道防線 + 「假說產生」的定位。
# 輸入：output/gbm4_final.rds（06a 跑過即有；此處以 OPC 為例——它在 06a 因兩部位都夠的病人只有
#       1 位而被跳過，正是「不能 pseudobulk」的實例）
# 輸出：output/06b_de_<型別>_cell_level.csv、output/figs/06b_*.png
# 時間：約 3–5 分鐘
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm4 <- readRDS("output/gbm4_final.rds")
gbm4[["RNA"]] <- JoinLayers(gbm4[["RNA"]])
if (!"data" %in% Layers(gbm4[["RNA"]])) gbm4 <- NormalizeData(gbm4)
gbm4$tissue <- factor(gbm4$tissue, levels = c("Periphery", "Tumor"))
gbm4$type <- ifelse(gbm4$malignant == "malignant", "Malignant", gbm4$celltype_author)

## ---- 1. when-you-land-here ------------------------------------------ Q3 頁 48
# 先誠實面對：這型別的樣本結構長什麼樣？
TYPE <- "OPC"                                        # 06a 的 de.summary 裡被跳過的那型
obj  <- subset(gbm4, cells = colnames(gbm4)[gbm4$type == TYPE])
print(table(obj$patient, obj$tissue))
# 讀法：兩部位都 ≥ 20 顆的病人只有 BT_S2——配對 pseudobulk 不成立。
# 從這裡開始，所有結果的定位都是「假說產生」，不是可發表的差異基因表。

## ---- 2. cell-level-de-with-covariate -------------------------------- Q3 頁 48–49
Idents(obj) <- "tissue"
# (a) 先跑 Seurat 預設：Wilcoxon（大多數教學和論文的起手式）——當基準
de.wilcox <- FindMarkers(obj, ident.1 = "Tumor", ident.2 = "Periphery",
                         logfc.threshold = 0.25, min.pct = 0.1)   # test.use 預設 = "wilcox"
de.wilcox <- de.wilcox |> tibble::rownames_to_column("gene") |> arrange(p_val_adj)

# (b) 防線一：檢定時把病人放進模型。Wilcoxon 做不到，改用 MAST 的 latent.vars
#（只有一位病人時 latent.vars 沒有意義，拿掉即可——但那也代表結論只屬於這位病人）。
de <- FindMarkers(obj, ident.1 = "Tumor", ident.2 = "Periphery",
                  ## TODO ▶ 不能 pseudobulk 時用哪個檢定、把什麼放進共變量？（Q3 頁 48–49）
                  test.use = "____", latent.vars = "____",   # ★ 病人當共變量
                  logfc.threshold = 0.25, min.pct = 0.1)
de <- de |> tibble::rownames_to_column("gene") |> arrange(p_val_adj)
head(de, 15)

# (c) 兩者對照：Wilcoxon 沒扣病人效應，通常「顯著」更多；被 MAST 拉掉的那些多半是病人差異假扮的
cat(sprintf("padj < 0.05 —— Wilcoxon（Seurat 預設）：%d；MAST + patient：%d；前 15 名重疊：%d\n",
            sum(de.wilcox$p_val_adj < 0.05, na.rm = TRUE), sum(de$p_val_adj < 0.05, na.rm = TRUE),
            length(intersect(head(de.wilcox$gene, 15), head(de$gene, 15)))))
only.wilcox <- setdiff(head(de.wilcox$gene, 30), head(de$gene, 30))
cat("只在 Wilcoxon 前 30 名的基因（很可能是病人效應）：", paste(head(only.wilcox, 10), collapse = ", "), "\n")
# p 值仍然偏小（細胞不是獨立觀察值），MAST 的排序只是「比較可信」，不是可信。

## ---- 3. per-patient-consistency ------------------------------------- Q3 頁 49
# 防線二：逐病人算 logFC，方向一致的才留。這是 cell-level DE 最重要的一張表。
per.pat <- function(g) {
  sapply(unique(obj$patient), function(p) {
    o <- subset(obj, patient == p)
    if (min(table(factor(o$tissue, levels = levels(obj$tissue)))) < 5) return(NA_real_)  # 單邊 < 5 顆不算
    m <- LayerData(o, layer = "data")[g, ]
    mean(m[o$tissue == "Tumor"]) - mean(m[o$tissue == "Periphery"])
  })
}
top <- head(de$gene, 30)
cons <- t(sapply(top, per.pat))
cons.df <- data.frame(gene = top, cons,
                      n_evaluable = rowSums(!is.na(cons)),
                      n_same_dir  = rowSums(sign(cons) == sign(rowMeans(cons, na.rm = TRUE)), na.rm = TRUE))
print(cons.df)
# NA 不是錯誤：NA = 該病人某一側 < 5 顆、無法評估。整欄 NA 正是「為什麼不能 pseudobulk」的視覺化——
# 這個型別只有 BT_S2 可評估（n_evaluable = 1），一致性防線根本使不上力。
if (max(cons.df$n_evaluable) <= 1)
  message("  ⚠ 只有 ", sum(!is.na(cons[1, ])), " 位病人可評估——一致性無從檢查，結論只屬於這（幾）位病人，報告限制段落要明寫。")
de$consistent <- de$gene %in% cons.df$gene[cons.df$n_evaluable >= 2 & cons.df$n_same_dir == cons.df$n_evaluable]
# 只有一位病人可評估時（n_evaluable = 1），「一致性」無從談起——這正是要誠實寫進報告的限制。

## ---- 4. permutation-check ------------------------------------------- Q3 頁 49
# 防線三：標籤置換。把 tissue 標籤在「病人內」隨機打散重跑，看能撈到幾個「顯著」基因。
# 真實訊號應該遠多於置換後的數目；差不多多，就代表你看到的多半是假重複灌出來的。
n.perm <- 20                                          # 教學用 20 次；正式分析 ≥ 100
perm.sig <- replicate(n.perm, {
  o <- obj
  o$tissue <- ave(as.character(o$tissue), o$patient, FUN = sample)  # 病人內打散
  Idents(o) <- factor(o$tissue, levels = c("Periphery", "Tumor"))
  d <- suppressWarnings(FindMarkers(o, ident.1 = "Tumor", ident.2 = "Periphery",
                                    test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1))
  sum(d$p_val_adj < 0.05, na.rm = TRUE)
})
cat(sprintf("真實標籤 padj<0.05：%d；置換後中位數：%.0f（範圍 %d–%d）\n",
            sum(de$p_val_adj < 0.05, na.rm = TRUE), median(perm.sig), min(perm.sig), max(perm.sig)))
p <- ggplot(data.frame(n = perm.sig), aes(n)) + geom_histogram(bins = 15, fill = "grey70") +
     geom_vline(xintercept = sum(de$p_val_adj < 0.05, na.rm = TRUE), colour = "#D62728", linewidth = 1) +
     theme_classic() + labs(x = "significant genes under permuted labels", y = "count",
                            title = sprintf("%s: real (red) vs permuted", TYPE))
ggsave("output/figs/06b_permutation.png", p, width = 6, height = 4, dpi = 150, bg = "white")

## ---- 5. deliver-with-caveats ---------------------------------------- Q3 頁 49
# 交付：表上明寫方法與限制。這份表的用途是「挑基因去驗證」，不是結論。
de$method <- "cell-level MAST, patient as latent var; hypothesis-generating (insufficient replicates for pseudobulk)"
write.csv(de, sprintf("output/06b_de_%s_cell_level.csv", gsub("[^A-Za-z0-9]+", "_", TYPE)), row.names = FALSE)
sessionInfo()

# =====================================================================
# ▶ 練習 6b
#  6b-1 §2(c) 印出的「只在 Wilcoxon 前 30 名」基因，逐一畫 VlnPlot(split.by = "patient")：
#       它們的差異是不是主要來自某一位病人？這就是共變量在扣的東西。
#  6b-2 對 Malignant（06a 裡只有 2 對病人）也跑一次這支：§3 的一致性表跟 06a 的 pseudobulk 結果
#       （06_de_Malignant.csv，06a 產出）前 20 名重疊多少？兩種方法各抓到什麼對方沒有的？
#  6b-3 把 §4 的 n.perm 提高到 100：置換分布的右尾碰得到真實值嗎？寫一句「這型別的 DE 可信度」結論。
#  6b-4 只有一位病人兩個部位都有時（例如把 obj 換成 BT_S2 的 OPC），三道防線各剩哪些還能做？
#       報告裡的限制段落該怎麼寫？
#  進階 讀 muscat 的 vignette：它的 mixed model（dream / MM-dream）與本腳本的做法差在哪？
#       什麼情況值得升級到 mixed model？
# =====================================================================
