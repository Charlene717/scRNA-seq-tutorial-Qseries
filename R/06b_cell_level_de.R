# =====================================================================
# 06b_cell_level_de.R — 練習腳本 6b：不能 pseudobulk 時的差異表達（cell-level DE 的正確姿勢）
#
# 對應影片：Q3 頁 48–49（三道防線：病人共變量、逐病人一致性、標籤置換）
# 什麼時候用這支：pseudobulk 需要「每個條件 ≥ 2–3 個生物重複」。當病人太少（06a 的 de.summary
#   顯示某型別不到 MIN_PAIRS）、或整份資料只有一位病人時，退而在細胞層級做——但要知道
#   p 值是被灌大的（假重複），所以這支的重點不是 p 值，而是三道防線 + 「假說產生」的定位。
# 輸入：output/rds/06_gbm4_final.rds（06a 跑過即有；此處以 OPC 為例——它在 06a 因兩部位都夠的病人只有
#       1 位而被跳過，正是「不能 pseudobulk」的實例）
# 輸出：output/tables/06b_de_<型別>_cell_level.csv、output/figs/06b_*.png
# 時間：約 3–5 分鐘
# =====================================================================
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
gbm4 <- readRDS("output/rds/06_gbm4_final.rds")
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
                  test.use = "MAST", latent.vars = "patient",   # ★ 病人當共變量
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
cons.df$verdict <- ifelse(cons.df$n_evaluable < 2, "無法判斷（可評估病人 < 2）",
                   ifelse(cons.df$n_same_dir == cons.df$n_evaluable, "方向一致", "方向不一致"))
de$consistent <- de$gene %in% cons.df$gene[cons.df$n_evaluable >= 2 & cons.df$n_same_dir == cons.df$n_evaluable]
cat(sprintf("通過一致性防線的基因：%d / %d（n_evaluable < 2 一律不算通過）\n",
            sum(de$consistent), nrow(de)))
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
cat(sprintf("真實標籤 padj<0.05：%d；置換後中位數：%.0f（範圍 %d–%d，%d 次）\n",
            sum(de$p_val_adj < 0.05, na.rm = TRUE), median(perm.sig), min(perm.sig), max(perm.sig), n.perm))
# ★ 先確認這個置換有沒有意義，再看它的數字。
#   「病人內打散」保住了每位病人對兩個部位各貢獻幾顆——這張表置換前後完全一樣。
#   如果某位病人幾乎只出現在一邊，打散就幾乎動不到東西，產生的虛無分布是假的。
cat("\n每位病人對兩個部位的貢獻（置換不會改變這張表）：\n"); print(table(obj$patient, obj$tissue))
if (diff(range(perm.sig)) == 0)
  message("  ⚠ ", n.perm, " 次置換得到同一個數字——代表這個設計裡「病人」與「部位」幾乎重合，\n",
          "    病人內打散動不了什麼。這個置換檢定不能拿來當「訊號是真的」的證據；\n",
          "    它反而是在告訴你：這份資料分不開病人效應與部位效應。")
p <- ggplot(data.frame(n = perm.sig), aes(n)) + geom_histogram(bins = 15, fill = "grey70") +
     geom_vline(xintercept = sum(de$p_val_adj < 0.05, na.rm = TRUE), colour = "#D62728", linewidth = 1) +
     theme_classic() + labs(x = "significant genes under permuted labels", y = "count",
                            title = sprintf("%s: real (red) vs permuted", TYPE))
ggsave("output/figs/06b_permutation.png", p, width = 6, height = 4, dpi = 150, bg = "white")

## >>> 參考答案 ------------------------------------------------------
# 本課這份資料的結果（實跑驗證，TYPE = "OPC"）：
#   樣本結構：BT_S1 2/23、BT_S2 154/22、BT_S4 176/2、BT_S6 21/0（Periphery/Tumor）。
#     兩部位都 ≥ 20 顆的病人只有 BT_S2 一位——配對 pseudobulk 不成立，所以才落到這支腳本。
#   §2 Wilcoxon 83 個顯著、MAST + patient 67 個，前 15 名只重疊 3 個。
#     重疊這麼低本身就是訊息：把病人放進模型之後，排序幾乎重來一次。
#   只在 Wilcoxon 前 30 名的基因裡有 XIST——X 染色體去活化的長鏈非編碼 RNA，
#     它的高低是病人的性別，跟腫瘤核心或邊緣完全無關。這是「病人效應假扮成部位效應」最乾淨的例子。
#   §3 一致性：四位病人只有 BT_S2 可評估（n_evaluable = 1），三十個基因全部「無法判斷」。
#     通過一致性防線的基因是 0 個。表格裡整片 NA 不是壞掉，是這份資料的實情被畫出來了。
#
# ★ §4 的置換結果要特別小心讀。實跑得到「真實 67；置換後中位數 7，範圍 7–7」。
#   二十次置換得到完全相同的數字，這不是巧合，是設計本身的問題：
#   「病人內打散」保住了每位病人對兩個部位各貢獻幾顆，而這份資料裡
#     Tumor 組 47 顆：BT_S1 49%、BT_S2 47%、BT_S4 4%、BT_S6 0%
#     Periphery 組 353 顆：BT_S1 1%、BT_S2 44%、BT_S4 50%、BT_S6 6%
#   Tumor 組有一半是 BT_S1，Periphery 組有一半是 BT_S4——兩組的病人組成天差地遠，
#   而且置換前後一模一樣（打散只在病人內部發生，不會把細胞從一位病人搬到另一位）。任何在 BT_S1/BT_S2 與 BT_S4/BT_S6
#   之間有差異的基因，每一次置換都會被判顯著——所以次次都是 7。
#   結論：**這個置換檢定在這份資料上沒有鑑別力**，不能拿 67 > 7 當作「訊號是真的」。
#   它真正告訴你的是：這個設計裡病人與部位幾乎重合，兩者分不開。腳本會自動偵測並提醒。
#   什麼時候置換才有意義：每位病人在兩個部位都有相當數量的細胞——那才有東西可以打散。
## <<< 參考答案

## ---- 5. deliver-with-caveats ---------------------------------------- Q3 頁 49
# 交付：表上明寫方法與限制。這份表的用途是「挑基因去驗證」，不是結論。
de$method <- "cell-level MAST, patient as latent var; hypothesis-generating (insufficient replicates for pseudobulk)"
write.csv(de, sprintf("output/tables/06b_de_%s_cell_level.csv", gsub("[^A-Za-z0-9]+", "_", TYPE)), row.names = FALSE)
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
#  6b-5 §4 的置換為什麼二十次都得到同一個數字？把 table(obj$patient, obj$tissue) 抄下來，
#       算算看 Tumor 組與 Periphery 組各由哪些病人組成、比例差多少。
#       在這種設計下，「病人內打散」還剩多少自由度？置換檢定還能證明什麼？
#  進階 讀 muscat 的 vignette：它的 mixed model（dream / MM-dream）與本腳本的做法差在哪？
#       什麼情況值得升級到 mixed model？
# =====================================================================
