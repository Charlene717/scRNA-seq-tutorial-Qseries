# =====================================================================
# 06a_pseudobulk_gsea.R — 練習腳本 6a：組成分析、pseudobulk + 配對 DESeq2、每樣本點圖、GSEA
#
# 對應影片：Q3 頁 29–47（§1 組成與 propeller、§2 每型別 pseudobulk + DESeq2、§3 每樣本點圖與火山圖、§4 GSEA / ORA、§5 交付）
# 輸入：output/gbm4_malignant.rds（05_infercnv.R）；若尚未跑 05（例如 JAGS 還沒裝），
#       §0 會退而用 output/gbm4_unintegrated.rds + 作者的 Neoplastic 標籤當替代惡性標籤
# 輸出：output/06_de_<型別>.csv、06_de_summary_by_type.csv、06_gsea_all_types.csv、06_ora_go_all_types.csv、figs/06_*.png
# 時間：約 5–10 分鐘（enrichKEGG 需連網）
# 套件：本版新增 ashr、reshape2、ggrepel（CRAN）與 clusterProfiler、org.Hs.eg.db、enrichplot（Bioc）——
#       請先重跑 00_setup.R（已安裝的會自動略過），或執行下面的檢查提示
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
set.seed(1234)
need <- c("ashr", "reshape2", "ggrepel", "DESeq2", "fgsea", "msigdbr", "speckle", "EnhancedVolcano",
          "clusterProfiler", "org.Hs.eg.db", "enrichplot", "data.table")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("缺少套件：", paste(miss, collapse = ", "),
                       "\n  → 請先重跑 R/00_setup.R（或 install.packages / BiocManager::install 上面這些）", call. = FALSE)
## ---- 0. load ---------------------------------------------------------
# 正式流程用 05 的 inferCNV 結果；05 還沒跑完（如 JAGS 未安裝）時用作者標籤代替，方便先測試 06–08
if (file.exists("output/gbm4_malignant.rds")) {
  gbm4 <- readRDS("output/gbm4_malignant.rds")
} else {
  warning("找不到 output/gbm4_malignant.rds（05 尚未完成）；改用 celltype_author == 'Neoplastic' 當替代惡性標籤。\n",
          "  裝好 JAGS、跑完 05 之後請重跑 06–08，結果會以 inferCNV 版本為準。", call. = FALSE)
  gbm4 <- readRDS("output/gbm4_unintegrated.rds")
  gbm4$malignant <- ifelse(gbm4$celltype_author == "Neoplastic", "malignant", "normal")
}
gbm4$tissue <- factor(gbm4$tissue, levels = c("Periphery", "Tumor"))   # 係數 = Tumor(核心) vs Periphery(邊緣)

## ---- 1. composition ------------------------------------------------ Q3 頁 29–31
gbm4$type <- ifelse(gbm4$malignant == "malignant", "Malignant", gbm4$celltype_author)
prop <- gbm4@meta.data |> dplyr::count(patient, tissue, type) |>
        group_by(patient, tissue) |> mutate(frac = n / sum(n)) |> ungroup()
p <- ggplot(prop, aes(tissue, frac, fill = type)) + geom_col(width = .8) +
     facet_wrap(~ patient, nrow = 1) + theme_classic() + labs(y = "fraction of cells")
ggsave("output/figs/06_composition.png", p, width = 12, height = 4.5, dpi = 150, bg = "white")

# 檢定：比例是組成資料，以「樣本」為單位（8 個），配對設計
library(speckle); library(limma)
sample.id <- paste(gbm4$patient, gbm4$tissue, sep = "_")
## TODO ▶ 組成資料為什麼要轉換？用哪種？（Q3 頁 29–31）
props <- getTransformedProps(clusters = gbm4$type, sample = sample.id, transform = "____")
# 注意：病人 ID 本身含底線（BT_S1），所以「病人」要拿掉最後一段（部位），不能用第一個底線切
sn <- colnames(props$TransformedProps)
design <- model.matrix(~ 0 + tissue + patient,
                       data = data.frame(tissue  = factor(sub(".*_", "", sn), levels = c("Periphery", "Tumor")),
                                         patient = sub("_[^_]*$", "", sn)))
fit <- lmFit(props$TransformedProps, design)
fit <- eBayes(contrasts.fit(fit, makeContrasts(tissueTumor - tissuePeriphery, levels = design)))
topTable(fit, n = Inf)                                # 每種型別：核心 vs 邊緣的比例差異（配對）
# 解讀：Malignant 在核心較高、Oligodendrocyte 在邊緣較高——但記得封閉性：一種變少其他必變多。

## ---- 2. pseudobulk-de（每種細胞型別各自比）-------------------------- Q3 頁 33–36
# 為什麼要「每種型別各自比」：核心 vs 邊緣的差異在惡性細胞、巨噬細胞、寡樹突裡各不相同。
# 把所有細胞混在一起比，得到的是「組成差異」（邊緣的正常腦細胞多），不是任何一種細胞的變化。
library(DESeq2)
## TODO ▶ 每個「病人 × 部位」至少幾顆細胞才算一個 pseudobulk 樣本？（Q3 頁 34）
MIN_CELLS <- ____        # 每個「病人 × 部位」至少要有這麼多顆細胞才算一個可信的 pseudobulk 樣本
MIN_PAIRS <- 2         # 至少要有這麼多位病人兩個部位都有樣本，配對設計才成立

# 把一種細胞型別的 pseudobulk + 配對 DESeq2 包成函式，之後對每種型別呼叫一次
pb_de <- function(obj, type, min.cells = MIN_CELLS, min.pairs = MIN_PAIRS) {
  sub <- subset(obj, cells = colnames(obj)[obj$type == type])   # 用 cells= 避免 subset 的欄名/變數名混淆
  n   <- table(sub$patient, sub$tissue)
  keep <- which(n >= min.cells, arr.ind = TRUE)         # 細胞太少的「病人 × 部位」整格丟掉
  ok.pat <- names(which(rowSums(n >= min.cells) == 2))   # 兩個部位都夠的病人
  if (length(ok.pat) < min.pairs) {
    message(sprintf("  %-16s 跳過：兩部位都 ≥ %d 顆的病人只有 %d 位", type, min.cells, length(ok.pat)))
    return(NULL)
  }
  sub <- subset(sub, cells = colnames(sub)[sub$patient %in% ok.pat])
  pb  <- as.matrix(AggregateExpression(sub, assays = "RNA", group.by = c("patient", "tissue"))$RNA)
  stopifnot("pseudobulk 矩陣不是整數：你加總到了錯的資料層" = sum(pb != round(pb)) == 0)
  coldata <- data.frame(patient = sub("_[^_]*$", "", colnames(pb)),   # 病人 ID 含底線，從最後一段切
                        tissue  = factor(sub(".*_", "", colnames(pb)), levels = c("Periphery", "Tumor")),
                        row.names = colnames(pb))
  ## TODO ▶ 配對設計的公式怎麼寫？（Q3 頁 33–34）
  dds <- DESeqDataSetFromMatrix(pb[rowSums(pb) >= 10, ], coldata, design = ~ ____ + ____)
  dds <- DESeq(dds, quiet = TRUE)
  raw <- results(dds, name = "tissue_Tumor_vs_Periphery")             # 未收縮：stat 給 GSEA 排序、火山圖
  ## TODO ▶ 小 n 的 lfcShrink 用哪種收縮？為什麼不用 apeglm？（Q3 頁 36）
  shr <- lfcShrink(dds, coef = "tissue_Tumor_vs_Periphery", type = "____", quiet = TRUE)
  #   收縮用 ashr 而不是 apeglm：樣本只有 4–8 個時 apeglm 的先驗太強，會把幾乎所有 LFC 壓成同一個值
  #   （火山圖會變成一條直直的柱子）；ashr 對小 n 穩定得多。報告效應量用收縮後的 LFC。
  df <- data.frame(gene = rownames(raw), baseMean = raw$baseMean,
                   log2FC = raw$log2FoldChange, log2FC_shrunk = shr$log2FoldChange,
                   stat = raw$stat, pvalue = raw$pvalue, padj = raw$padj) |> arrange(padj)
  list(type = type, n = n, pb = pb, coldata = coldata, dds = dds, res = df,
       n.sig = sum(df$padj < 0.05, na.rm = TRUE), n.pairs = length(ok.pat))
}

types <- names(which(table(gbm4$type) >= 100))          # 細胞太少的型別連試都不用試
types <- setdiff(types, "Unassigned")
de <- lapply(types, function(ty) pb_de(gbm4, ty)); names(de) <- types
de <- Filter(Negate(is.null), de)
de.summary <- data.frame(type = names(de), n_pairs = sapply(de, `[[`, "n.pairs"),
                         n_sig_padj05 = sapply(de, `[[`, "n.sig"))
print(de.summary)                                         # 每種型別各有幾對病人、幾個顯著基因
for (ty in names(de))
  write.csv(de[[ty]]$res, sprintf("output/06_de_%s.csv", gsub("[^A-Za-z0-9]+", "_", ty)), row.names = FALSE)
stopifnot("Malignant 沒有通過樣本數門檻——請先看 table(gbm4$patient, gbm4$tissue, gbm4$type)" = "Malignant" %in% names(de))
mal.res <- de[["Malignant"]]$res
head(mal.res, 15)

# 對照組（雷本體）：cell-level Wilcoxon，看 p 值膨脹多少（只看惡性細胞）
mal <- subset(gbm4, subset = type == "Malignant"); Idents(mal) <- "tissue"
naive <- FindMarkers(mal, ident.1 = "Tumor", ident.2 = "Periphery", logfc.threshold = 0, min.pct = 0.1)
cat("Malignant  cell-level padj < 0.05：", sum(naive$p_val_adj < 0.05),
    " vs pseudobulk：", de[["Malignant"]]$n.sig, "\n")
# 注意：如果 pseudobulk 反而比 cell-level 多，通常是某些「病人 × 部位」細胞太少、pseudobulk 樣本是噪音
# —— 這就是 MIN_CELLS 存在的理由；也可以把它調高到 30–50 再看一次。

## ---- 3. per-sample-plot + 火山圖 ----------------------------------- Q3 頁 37–38
plot.gene <- function(d, g) {
  cpm <- log2(t(t(d$pb) / colSums(d$pb)) * 1e6 + 1)
  df <- data.frame(expr = cpm[g, ], d$coldata)
  ggplot(df, aes(tissue, expr, group = patient, colour = patient)) +
    geom_line() + geom_point(size = 3) + theme_classic() + ggtitle(g) + labs(y = "log2 CPM (pseudobulk)")
}
top.genes <- head(mal.res$gene[!is.na(mal.res$padj)], 4)
p <- wrap_plots(lapply(top.genes, function(g) plot.gene(de[["Malignant"]], g)), nrow = 1)
ggsave("output/figs/06_per_sample_dots_malignant.png", p, width = 14, height = 4, dpi = 150, bg = "white")
# 看圖：四位病人方向一致嗎？一致的才是可信的差異。

# 火山圖：每種細胞型別一張。x = 未收縮 log2FC（收縮版會把點壓到中間，形狀失真），y = -log10(padj)
library(EnhancedVolcano)
# 配色用標準的紅／藍／灰：紅 = 核心較高且顯著、藍 = 邊緣較高且顯著、灰 = 不顯著
volcano <- function(d, fc = 1, p = 0.05) {
  r  <- d$res
  kv <- ifelse(!is.na(r$padj) & r$padj < p & r$log2FC >  fc, "#D62728",
        ifelse(!is.na(r$padj) & r$padj < p & r$log2FC < -fc, "#1F77B4", "grey70"))
  names(kv) <- ifelse(kv == "#D62728", "Up in core", ifelse(kv == "#1F77B4", "Up in periphery", "NS"))
  ## TODO ▶ 火山圖的 x 用哪種 log2FC？y 用 p 還是 padj？（Q3 頁 38）
  EnhancedVolcano(r, lab = r$gene, x = "____", y = "____",
                  pCutoff = p, FCcutoff = fc, pointSize = 1.2, labSize = 3,
                  colCustom = kv, colAlpha = 0.7,
                  title = paste0(d$type, ": core vs periphery"),
                  subtitle = sprintf("%d pairs; %d genes padj < %.2f", d$n.pairs, d$n.sig, p),
                  legendPosition = "bottom", drawConnectors = TRUE, max.overlaps = 20)
}
vol <- lapply(de, volcano)
ggsave("output/figs/06_volcano_all_types.png", wrap_plots(vol, ncol = 2),
       width = 14, height = 6.5 * ceiling(length(vol) / 2), dpi = 120, bg = "white", limitsize = FALSE)
for (ty in names(de))
  ggsave(sprintf("output/figs/06_volcano_%s.png", gsub("[^A-Za-z0-9]+", "_", ty)), vol[[ty]],
         width = 8, height = 7, dpi = 150, bg = "white")
# 讀火山圖（Q3 頁 38）：右上 = 核心較高且顯著、左上 = 邊緣較高且顯著；中間高高的一根 = 效應小但 p 小，
# 通常是表現量高的基因（baseMean 大），要回頭看每樣本點圖確認四位病人方向是否一致。
# 形狀異常的訊號：所有點擠成一條直柱（LFC 被過度收縮）、或只有正的一側（某個部位樣本幾乎沒有細胞）。

## ---- 4. enrichment：GSEA（Hallmark / GO / KEGG）+ ORA（GO / KEGG）--- Q3 頁 39–45
# 兩種方法回答不同問題：
#   GSEA：「整個排序清單裡，這個基因集是不是偏向一端？」不切閾值，全部基因都參與 —— 小 n 的首選
#   ORA ：「顯著基因裡，這個基因集的比例是不是高於背景？」要先切閾值，背景必須是「被檢定的基因」
library(fgsea); library(msigdbr); library(data.table)
msig <- function(coll, sub = NULL) {                        # 相容新舊版 msigdbr 的參數名
  if ("collection" %in% names(formals(msigdbr))) {
    m <- msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  } else {
    m <- msigdbr(species = "Homo sapiens", category = coll, subcategory = sub)
  }
  split(m$gene_symbol, m$gs_name)
}
gene.sets <- list(Hallmark = msig("H"),
                  GO_BP    = msig("C5", "GO:BP"),
                  KEGG     = tryCatch(msig("C2", "CP:KEGG_LEGACY"),          # 新版 msigdbr 的名稱
                                      error = function(e) msig("C2", "CP:KEGG")))  # 舊版用這個
sapply(gene.sets, length)

rank_stat <- function(res) {                                # 排序統計量：Wald stat（不是 shrunken LFC）
  ## TODO ▶ GSEA 的排序統計量用哪一欄？為什麼不是 shrunken LFC？（Q3 頁 36、41）
  r <- setNames(res$____, res$gene); r <- r[!is.na(r)]
  r <- r + rnorm(length(r), sd = 1e-6)                      # 打散同分
  sort(r, decreasing = TRUE)
}
run_gsea <- function(d, sets, minSize = 15, maxSize = 500) {
  r <- rank_stat(d$res)
  out <- lapply(names(sets), function(nm) {
    ## TODO ▶ fgsea 基因集大小的上下限（Q3 頁 41）
    g <- fgsea(sets[[nm]], r, minSize = ____, maxSize = ____, nproc = 1)   # eps 留預設
    g$collection <- nm; g$type <- d$type; g
  })
  data.table::rbindlist(out)[order(padj)]
}
gsea.all <- data.table::rbindlist(lapply(de, run_gsea, sets = gene.sets))
gsea.all[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]    # list 欄轉字串才能存 CSV
write.csv(gsea.all, "output/06_gsea_all_types.csv", row.names = FALSE)
# 每種型別、每個資料庫各看前 5：
gsea.all[padj < 0.05, .SD[order(-abs(NES))][1:min(5, .N)], by = .(type, collection)][, .(type, collection, pathway, NES, padj)]

# 圖 A：Hallmark NES 熱圖 —— 列 = pathway、欄 = 細胞型別。一眼看出「缺氧只在惡性細胞升」還是「每種細胞都升」
hm.top <- gsea.all[collection == "Hallmark" & padj < 0.05, unique(pathway)]
hm.mat <- reshape2::acast(gsea.all[collection == "Hallmark" & pathway %in% hm.top], pathway ~ type, value.var = "NES")
hm.df  <- reshape2::melt(hm.mat, varnames = c("pathway", "type"), value.name = "NES")
hm.df$sig <- with(hm.df, ifelse(is.na(NES), "", "*"))
p <- ggplot(hm.df, aes(type, gsub("HALLMARK_", "", pathway), fill = NES)) + geom_tile(colour = "white") +
     scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#D62728", na.value = "grey92") +
     theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
     labs(x = NULL, y = NULL, title = "Hallmark GSEA: NES by cell type (core vs periphery)",
          subtitle = "red = higher in core; blue = higher in periphery; grey = not tested / NA")
ggsave("output/figs/06_gsea_hallmark_heatmap.png", p, width = 9, height = 0.28 * length(hm.top) + 2.5, dpi = 150, bg = "white")

# 圖 B：每種型別的 GSEA 條圖（三個資料庫各取 |NES| 最大的前 8）
gsea_bar <- function(ty) {
  g <- gsea.all[type == ty & padj < 0.05][order(-abs(NES))][, .SD[1:min(8, .N)], by = collection]
  if (!nrow(g)) return(NULL)
  g$label <- substr(gsub("^(HALLMARK|GOBP|KEGG)_", "", g$pathway), 1, 45)
  ggplot(g, aes(reorder(label, NES), NES, fill = NES > 0)) + geom_col() + coord_flip() +
    facet_wrap(~ collection, scales = "free_y", ncol = 1) + theme_classic() +
    scale_fill_manual(values = c(`TRUE` = "#D62728", `FALSE` = "#1F77B4"), labels = c("periphery", "core"), name = "higher in") +
    labs(x = NULL, y = "NES", title = paste0(ty, ": GSEA, padj < 0.05"))
}
for (ty in names(de)) { p <- gsea_bar(ty); if (!is.null(p))
  ggsave(sprintf("output/figs/06_gsea_bar_%s.png", gsub("[^A-Za-z0-9]+", "_", ty)), p, width = 8, height = 9, dpi = 150, bg = "white") }

# 圖 C：enrichment plot（running score）—— 一條 pathway 的「證據長什麼樣」
r.mal <- rank_stat(mal.res)
p <- plotEnrichment(gene.sets$Hallmark[["HALLMARK_HYPOXIA"]], r.mal) + labs(title = "Malignant: HALLMARK_HYPOXIA (core vs periphery)")
ggsave("output/figs/06_gsea_hypoxia_malignant.png", p, width = 6, height = 4, dpi = 150, bg = "white")
# 讀法：黑色 tick 是基因集成員在排序中的位置；綠線是 running score；峰值在左 = 富集在「核心較高」那端。
# leading edge = 峰值之前的成員，就是真正在動的基因：
mal.hyp <- gsea.all[type == "Malignant" & pathway == "HALLMARK_HYPOXIA"]
strsplit(mal.hyp$leadingEdge, ";")[[1]][1:15]

# ORA：clusterProfiler 的 GO（離線，org.Hs.eg.db）與 KEGG（要連網）。上調、下調分開做
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
run_ora <- function(d, direction = c("up", "down"), fc = 0.5, p = 0.05) {
  direction <- match.arg(direction)
  sig <- d$res$gene[!is.na(d$res$padj) & d$res$padj < p & (if (direction == "up") d$res$log2FC > fc else d$res$log2FC < -fc)]
  ## TODO ▶ ORA 的背景基因集應該是什麼？（Q3 頁 43、45）
  universe <- ____                                     # ★ 背景 = 被檢定的基因，不是全基因組
  if (length(sig) < 10) { message("  ", d$type, " ", direction, "：顯著基因 < 10，ORA 跳過（改看 GSEA）"); return(NULL) }
  go <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe,
                 pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 15, maxGSSize = 500)
  go <- clusterProfiler::simplify(go, cutoff = 0.7)          # GO 詞條高度重疊，先合併相似的
  #   ↑ 寫全名：igraph（CellChat 會載入）也有 simplify，同一個 R session 跑過 07 會被蓋掉
  ids <- bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  uni <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  kegg <- tryCatch({                                         # enrichKEGG 要連網抓 KEGG 資料庫
    k <- enrichKEGG(ids$ENTREZID, organism = "hsa", universe = uni$ENTREZID, minGSSize = 15)
    if (is.null(k)) NULL else setReadable(k, org.Hs.eg.db, keyType = "ENTREZID")
  }, error = function(e) { message("  KEGG 無結果或無法連網：", conditionMessage(e)); NULL })
  list(type = d$type, direction = direction, n = length(sig), go = go, kegg = kegg)
}
ora <- list()
for (ty in names(de)) for (dr in c("up", "down")) ora[[paste(ty, dr)]] <- run_ora(de[[ty]], dr)
ora <- Filter(Negate(is.null), ora)

ora_dot <- function(o, which = c("go", "kegg"), n = 12) {
  which <- match.arg(which); e <- o[[which]]
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  dotplot(e, showCategory = n, font.size = 9) +
    ggtitle(sprintf("%s, %s in core (%d genes): %s ORA", o$type, o$direction, o$n, toupper(which)))
}
for (nm in names(ora)) {
  f <- gsub("[^A-Za-z0-9]+", "_", nm)
  p <- ora_dot(ora[[nm]], "go");   if (!is.null(p)) ggsave(sprintf("output/figs/06_ora_go_%s.png",   f), p, width = 8, height = 6, dpi = 150, bg = "white")
  p <- ora_dot(ora[[nm]], "kegg"); if (!is.null(p)) ggsave(sprintf("output/figs/06_ora_kegg_%s.png", f), p, width = 8, height = 6, dpi = 150, bg = "white")
}
# 讀 dotplot：x = GeneRatio（顯著基因裡落在此詞條的比例）、點大小 = 命中數、顏色 = padj。
# 只看前幾名、GeneRatio 高又 padj 小的；命中數 3–5 個的詞條再顯著也先別寫進結論。
# 其他常用圖：cnetplot(go, showCategory = 5)（詞條–基因網路）、emapplot(pairwise_termsim(go))（詞條相似網路）
if (!is.null(ora[["Malignant up"]]) && nrow(as.data.frame(ora[["Malignant up"]]$go))) {
  p <- cnetplot(ora[["Malignant up"]]$go, showCategory = 5) +
       ggtitle("Malignant, up in core: GO BP term-gene network")
  ggsave("output/figs/06_ora_cnet_malignant_up.png", p, width = 10, height = 8, dpi = 150, bg = "white")
}
ora.tab <- data.table::rbindlist(lapply(ora, function(o) {
  g <- as.data.frame(o$go); if (!nrow(g)) return(NULL)
  data.frame(type = o$type, direction = o$direction, db = "GO_BP", g[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")])
}), fill = TRUE)
write.csv(ora.tab, "output/06_ora_go_all_types.csv", row.names = FALSE)

# 已知答案驗證：Malignant 的核心端應偏向 HYPOXIA / GLYCOLYSIS（GSEA 正 NES；ORA up 應看到 response to hypoxia）。
# 方向合理才往下讀 leading edge；方向反了，先回去檢查 tissue 的 levels 順序。
sessionInfo()

## ---- 5. deliverables ----------------------------------------------- Q3 頁 46
# 交出去的東西：每種型別一份 DE 表（含 raw 與 shrunken LFC、baseMean、padj）、一張火山圖、
# GSEA 總表（含 leading edge）、ORA 總表；全部由本腳本重生，output/ 裡沒有手工檔。
write.csv(de.summary, "output/06_de_summary_by_type.csv", row.names = FALSE)
ggsave("output/figs/06_volcano_malignant.pdf", vol[["Malignant"]], width = 8, height = 7, bg = "white")
ggsave("output/figs/06_gsea_hypoxia.pdf", plotEnrichment(gene.sets$Hallmark[["HALLMARK_HYPOXIA"]], r.mal) + labs(title = "Malignant: Hypoxia"),
       width = 6, height = 4, bg = "white")
saveRDS(gbm4, "output/gbm4_final.rds")

# =====================================================================
# ▶ 練習 6
#  6-1 把 pb_de 裡的 design 改成 ~ tissue（不配對）重跑 Malignant。padj < 0.05 的基因數變成幾個？
#      為什麼配對設計功效比較高？（提示：看 §3 的每樣本點圖，病人之間的基線差多少）
#  6-2 §2 的 cell-level Wilcoxon 顯著基因裡，有多少在 pseudobulk 也顯著？
#      挑三個「cell-level 極顯著、pseudobulk 不顯著」的基因畫每樣本點圖，它們長什麼樣？
#  6-3 把 MIN_CELLS 從 20 改成 5 再改成 50：de.summary 怎麼變？哪個門檻下 Malignant 的火山圖形狀最「正常」？
#  6-4 比較 Malignant 與 Immune cell 的 Hallmark 熱圖那一欄：哪些 pathway 兩種細胞同方向、哪些只在其中一種？
#      同方向的通常代表什麼（提示：微環境 vs 細胞內在程式）？
#  6-5 ORA 的 universe 改成 NULL（clusterProfiler 預設用全基因組當背景）重跑 Malignant up：
#      顯著詞條多了多少？哪一種背景才對，為什麼？
#  6-6 GSEA 的排序改用 log2FC_shrunk：前五名基因集變了嗎？fgsea 有沒有對 ties 發警告？
#  6-7 用 05 的 cnv.score 當共變量加進 design（~ patient + tissue + cnv）合理嗎？寫下你的判斷。
#  進階 用 CellChat（或 liana）比較核心 vs 邊緣的惡性細胞與免疫細胞之間的配體受體軸；
#      在動手前先寫下：這份資料符合「細胞通訊」那條路的資料前提嗎？（Q3 頁 54）
# =====================================================================
