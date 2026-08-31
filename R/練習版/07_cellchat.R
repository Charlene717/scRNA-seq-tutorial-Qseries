# =====================================================================
# 07_cellchat.R — 練習腳本 7：細胞通訊（CellChat）——每個樣本各跑一次、六種圖、兩條件比較、LIANA 交叉驗證
#
# 對應影片：Q3 頁 54–65（§1 跑一次 CellChat、§2 路徑層級與六種圖、§3 兩條件比較、§4 LIANA）
# 輸入：output/rds/06_gbm4_final.rds（06a_pseudobulk_gsea.R；含 malignant 標籤與 type 欄）
# 輸出：output/rds/07_cellchat/<patient>_<tissue>_min<MIN.CELLS>.rds、output/tables/07_liana_top500.csv、output/figs/07_*.pdf
# 時間：每個樣本約 5–15 分鐘（8 個樣本，建議先跑一位病人）；跑過的樣本會存成 rds，
#       第二次執行由 REUSE.RDS 直接讀回，只有 06 的輸出更新時才重算
# 安裝：devtools::install_github("jinworks/CellChat")；LIANA：remotes::install_github("saezlab/liana")
# =====================================================================
# ---------------------------------------------------------------------
# 【練習版】把 ____ 填上再執行。每個空格上方的「## TODO ▶」寫了要回答的問題與影片頁碼。
# 完整解答在上一層資料夾的同名檔案；建議先自己填，跑不通再對照。
# ---------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(CellChat); library(patchwork)
set.seed(1234)
## TODO ▶ 少於幾顆細胞的群不參與通訊分析？低於門檻的群是整組被移除，不是畫得淡一點（Q3 頁 56）
MIN.CELLS <- ____        # CellChat 建網時的細胞數門檻：低於這個數的群「整組」被移除，不是畫得淡一點
REUSE.RDS <- TRUE      # 已經跑過的樣本直接讀 output/rds/07_cellchat/*.rds（每個樣本 5–15 分鐘，重跑一輪要一小時）
                       # 安全性：只要 06 的輸出比快取新，就自動重跑那個樣本——不會像 inferCNV 那樣默默用舊結果
in.rds <- "output/rds/06_gbm4_final.rds"
gbm4 <- readRDS(in.rds)
gbm4[["RNA"]] <- JoinLayers(gbm4[["RNA"]])                # 04 之後 RNA 是按病人分層的；合併回單一 data 層，
if (!"data" %in% Layers(gbm4[["RNA"]])) gbm4 <- NormalizeData(gbm4)   # 否則 plotGeneExpression/VlnPlot 會拿 counts 畫
dir.create("output/rds/07_cellchat", showWarnings = FALSE, recursive = TRUE)   # 每個樣本一個 CellChat 物件，歸在 rds/ 底下
for (d in c("output/figs", "output/rds", "output/tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 通訊分析用的標籤：CNV 判定後的惡性 + 作者的正常型別；unresolved 排除
gbm4$cc_label <- ifelse(gbm4$malignant == "malignant", "Malignant",
                 ifelse(gbm4$celltype_author == "Immune cell", "Macro/MG",   # 本資料免疫細胞以髓系為主
                        gbm4$celltype_author))
gbm4 <- subset(gbm4, malignant != "unresolved")
table(gbm4$cc_label, paste(gbm4$patient, gbm4$tissue))                      # 每群 ≥ MIN.CELLS 顆才進得了網路

## ---- 1. run-per-sample --------------------------------------------- Q3 頁 55–56
run_cc <- function(obj, label = "cc_label") {
  obj$samples <- factor(paste(obj$patient, obj$tissue, sep = "_"))          # CellChat v2 要求 meta 有 samples 欄
  cc <- createCellChat(object = obj, group.by = label, assay = "RNA")       # 用 data 層（log-normalized）
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")        # 先只看分泌型；熟了可用全部
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  ## TODO ▶ CellChat 用什麼統計量代表一群的表現？（Q3 頁 55–56）
  cc <- computeCommunProb(cc, type = "____", population.size = TRUE)     # 注意大小寫 triMean；群大小校正
  cc <- filterCommunication(cc, min.cells = MIN.CELLS)   # 細胞數不足的群會被整組移除，見下面的存活表
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cc
}
samples <- unique(paste(gbm4$patient, gbm4$tissue, sep = "_"))
cc.all <- list()
for (s in samples) {                                                        # 8 個樣本；先跑一位病人也可以
  f <- paste0("output/rds/07_cellchat/", s, "_min", MIN.CELLS, ".rds")           # 門檻寫進檔名：改了 MIN.CELLS 就是另一份快取，
  if (REUSE.RDS && file.exists(f) && file.mtime(f) > file.mtime(in.rds)) {   # 不會拿舊門檻算出來的網路冒充新的
    cc.all[[s]] <- readRDS(f); cat("  ", s, "讀快取\n"); next
  }
  pt <- sub("_[^_]*$", "", s); ti <- sub(".*_", "", s)    # 病人 ID 含底線，從最後一段切
  obj <- subset(gbm4, patient == pt & tissue == ti)
  # 能不能跑，看的不是總細胞數，是「有幾群過得了 MIN.CELLS」——通訊至少要兩群才成立。
  # 用總數當門檻會誤殺：BT_S6_Tumor 全部只有 157 顆，但 Macro/MG 54、Malignant 90，兩群都夠。
  keep <- names(which(table(obj$cc_label) >= MIN.CELLS))
  if (length(keep) < 2) { cat("  ", s, "只有", length(keep), "群過得了門檻，跳過（通訊至少要兩群）\n"); next }
  cat(sprintf("   %s 重新計算：%d 顆細胞，可用的群 %s\n", s, ncol(obj), paste(keep, collapse = "、")))
  res <- tryCatch(run_cc(obj),                                               # 一個樣本失敗不該讓整支腳本停下來
                  error = function(e) { message("  ", s, " 跑不完：", conditionMessage(e)); NULL })
  if (is.null(res)) next
  cc.all[[s]] <- res; saveRDS(res, f)
}
# 哪幾群撐過 min.cells：後面畫得出什麼圖，看的就是這張表。
# filterCommunication 不是把細胞少的群畫淡一點，是把它整組從網路裡拿掉——
# 所以「某一條互動不見了」有兩種可能：真的沒有訊號，或者那一群根本沒進網路。這兩件事要分得開。
cc_groups <- function(x, min.cells = MIN.CELLS) { n <- table(x@idents); names(n)[n >= min.cells] }
cat("\n每個樣本通過 min.cells =", MIN.CELLS, "的細胞群：\n")
for (s in names(cc.all)) {
  n <- table(cc.all[[s]]@idents); n <- n[n > 0]
  cat(sprintf("  %-18s 保留 %s\n", s,
      paste(sprintf("%s(%d)", names(n)[n >= MIN.CELLS], n[n >= MIN.CELLS]), collapse = "、")))
  if (any(n < MIN.CELLS))
    cat(sprintf("  %-18s 移除 %s\n", "",
        paste(sprintf("%s(%d)", names(n)[n < MIN.CELLS], n[n < MIN.CELLS]), collapse = "、")))
}

if (!length(cc.all)) stop("沒有任何樣本跑完 CellChat：§1 的 ncol(obj) < 200 把八個樣本全擋掉了，把門檻放低再跑一次。")
DEMO <- if ("BT_S2_Tumor" %in% names(cc.all)) "BT_S2_Tumor" else names(cc.all)[1]
cc <- cc.all[[DEMO]]                                                        # 以下六種圖用一個樣本示範
DEMO.p <- sub("_[^_]*$", "", DEMO); DEMO.t <- sub(".*_", "", DEMO)
cat("\n六種圖的示範樣本：", DEMO, "\n")

## ---- 2. six-plots --------------------------------------------------- Q3 頁 57–62
groupSize <- as.numeric(table(cc@idents))
pdf("output/figs/07_1_circle.pdf", width = 10, height = 5); par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cc@net$count,  vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cc@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights")
dev.off()
# 圖 2：熱圖（誰送給誰）
p2 <- netVisual_heatmap(cc, measure = "weight", color.heatmap = "Reds"); pdf("output/figs/07_2_heatmap.pdf", 6, 5); print(p2); dev.off()
# 圖 4：signaling role（先做，找主角）
cc@netP$pathways
p4a <- netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing", height = 8)
p4b <- netAnalysis_signalingRole_heatmap(cc, pattern = "incoming", height = 8)
pdf("output/figs/07_4_roles.pdf", 12, 7); print(p4a + p4b); dev.off()
netAnalysis_signalingRole_scatter(cc); ggsave("output/figs/07_4_scatter.pdf", width = 5, height = 4, bg = "white")
# 圖 3：挑一條路徑拆開（以 SPP1 為例；名稱一定要在 cc@netP$pathways 裡，否則整段會報錯）
pw <- if ("SPP1" %in% cc@netP$pathways) "SPP1" else cc@netP$pathways[1]
if (pw != "SPP1") cat("這個樣本推不出 SPP1，改用", pw, "示範\n")
pdf("output/figs/07_3_pathway.pdf", 10, 5); par(mfrow = c(1, 2))
netVisual_aggregate(cc, signaling = pw, layout = "circle")
netVisual_aggregate(cc, signaling = pw, layout = "chord")
dev.off()
netAnalysis_contribution(cc, signaling = pw); ggsave("output/figs/07_3_contribution.pdf", width = 5, height = 3, bg = "white")
netVisual_heatmap(cc, signaling = pw, color.heatmap = "Reds")
plotGeneExpression(cc, signaling = pw); ggsave("output/figs/07_3_genes.pdf", width = 8, height = 5, bg = "white")   # 驗證表現
netAnalysis_signalingRole_network(cc, signaling = pw, width = 8, height = 2.5)
# 圖 5：bubble（寫進論文的那張）
# 群名一律先跟 cc_groups() 取交集：寫了一個這個樣本裡沒有（或細胞數不足被移除）的群，整行就報錯
g5.src <- intersect(c("Macro/MG", "Malignant"), cc_groups(cc))
g5.tgt <- intersect(c("Malignant", "Macro/MG", "Vascular"), cc_groups(cc))
if (length(g5.src) && length(g5.tgt)) {
  netVisual_bubble(cc, sources.use = g5.src, targets.use = g5.tgt, remove.isolate = TRUE)
  ggsave("output/figs/07_5_bubble.pdf", width = 7, height = 8, bg = "white")
} else cat("這個樣本裡", DEMO, "沒有足夠的 Macro/MG 或 Malignant，跳過圖 5\n")
# 只看幾條路徑：名稱一定要在 cc@netP$pathways 裡（每個樣本推得出的路徑不同，硬寫 MIF 這種名稱會報錯）
sig3 <- head(intersect(c("SPP1", "MIF", "VEGF", "PTN", "TGFb"), cc@netP$pathways), 3)
if (length(sig3) == 0) sig3 <- head(cc@netP$pathways, 3)
netVisual_bubble(cc, signaling = sig3, remove.isolate = TRUE)
ggsave("output/figs/07_5_bubble_paths.pdf", width = 7, height = 6, bg = "white")
# 回 Seurat 驗證具體的一對：CellChat 的機率是推出來的，配體與受體到底表現在誰身上要自己看
feats <- intersect(c("SPP1", "CD44"), rownames(gbm4))
if (length(feats)) VlnPlot(subset(gbm4, patient == DEMO.p & tissue == DEMO.t),
                           features = feats, group.by = "cc_label", pt.size = 0)

## ---- 3. compare-conditions ------------------------------------------ Q3 頁 63–64
PAIR <- "BT_S2"                                                             # 同一位病人的核心 vs 邊緣
need <- paste0(PAIR, c("_Tumor", "_Periphery"))
if (!all(need %in% names(cc.all)))
  stop("§3 需要 ", paste(need, collapse = " 與 "), " 兩個樣本都跑完。\n",
       "  目前 cc.all 裡有：", paste(names(cc.all), collapse = ", "), "\n",
       "  如果你只跑了一位病人，把 §1 的 samples 放寬再跑一次，或把 PAIR 改成你跑過的病人。")
cc.list <- list(Core = cc.all[[need[1]]], Periphery = cc.all[[need[2]]])
cc.m <- mergeCellChat(cc.list, add.names = names(cc.list))
g1 <- compareInteractions(cc.m, show.legend = FALSE, group = c(1, 2))
g2 <- compareInteractions(cc.m, show.legend = FALSE, group = c(1, 2), measure = "weight")
g1 + g2; ggsave("output/figs/07_6_compare.pdf", width = 6, height = 3, bg = "white")
pdf("output/figs/07_6_diff.pdf", 10, 5); par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cc.m, weight.scale = TRUE)                       # 紅：Periphery > Core；藍：反之
netVisual_diffInteraction(cc.m, weight.scale = TRUE, measure = "weight")
dev.off()
rankNet(cc.m, mode = "comparison", stacked = TRUE, do.stat = TRUE); ggsave("output/figs/07_6_rankNet.pdf", width = 5, height = 6, bg = "white")
# 兩條件並排的 bubble 有兩道門檻，順序不能顛倒：
#   ① 來源與目標兩群都要在兩個條件裡「活著」（細胞數 ≥ MIN.CELLS，否則整組被移除）
#   ② 兩邊都要推得出顯著互動（只有一邊有的時候，CellChat 會丟 seq 的 'by' 錯誤，是已知 bug）
# 這份資料卡在第①關：BT_S2 的 Periphery 只有 13 顆 Malignant，低於 20，
# 建網時整群被拿掉（log 會寫「91.4% interactions are removed」）。
# 所以先自己問「兩群都在嗎」，不要等 CellChat 丟錯誤才發現——錯誤訊息說的是「沒有互動」，
# 真正的原因卻是「沒有細胞」，這兩句話在論文裡的意思完全不同。
SRC <- "Macro/MG"; TGT <- "Malignant"
alive <- vapply(cc.list, function(x) all(c(SRC, TGT) %in% cc_groups(x)), logical(1))
cat(sprintf("\n%s → %s：%s\n", SRC, TGT,
    paste(sprintf("%s＝%s", names(alive), ifelse(alive, "兩群都在", "有一群細胞數不足，已被移除")), collapse = "，")))

# 一組「來源 → 目標」的兩條件 bubble：先試合併物件；②那個 bug 一觸發就退回兩張單樣本並排；
# 兩種都畫不出來，才把原因說清楚。共同細胞群那一組用的是同一個函式。
bubble_pair <- function(src, tgt, file, ttl) {
  ok <- vapply(cc.list, function(x) all(c(src, tgt) %in% cc_groups(x)), logical(1))
  one <- function(x, k)
    tryCatch(netVisual_bubble(x, sources.use = src, targets.use = tgt, remove.isolate = TRUE) + ggtitle(k),
             error = function(e) { message("  ", ttl, "／", k, " 畫不出來：", conditionMessage(e)); NULL })
  p <- NULL
  if (all(ok))
    p <- tryCatch(netVisual_bubble(cc.m, sources.use = src, targets.use = tgt, comparison = c(1, 2), angle.x = 45),
                  error = function(e) {
                    message("  ", ttl, "：合併物件的並排失敗（", conditionMessage(e), "）——改畫單樣本對照")
                    NULL })
  if (is.null(p)) {
    ps <- Filter(Negate(is.null), lapply(names(cc.list), function(k) one(cc.list[[k]], k)))
    if (length(ps)) {
      p <- patchwork::wrap_plots(ps, nrow = 1)
      if (length(ps) < length(cc.list))        # 只畫得出一邊時，圖上一定要寫清楚另一邊為什麼是空的——
        p <- p + patchwork::plot_annotation(   # 不然這張圖被單獨貼進投影片，就會被讀成「另一邊沒有通訊」
          subtitle = sprintf("只有 %d/%d 個條件畫得出來；缺的那一邊是細胞數不足（< %d）整群被移除，不是沒有訊號",
                             length(ps), length(cc.list), MIN.CELLS))
    }
  }
  if (is.null(p)) {
    cat(ttl, "：兩個條件都畫不出來。這一格本身就是結論，不是失敗：\n",
        "  細胞數不足的群會被整組移除，圖上的空白代表「沒有東西可以比」，不代表「兩邊沒有差異」。\n",
        "  下一步只有兩條路：降低 MIN.CELLS（代價是機率估計不穩），或換一組兩邊都夠大的來源→目標。\n", sep = "")
    return(invisible(FALSE))
  }
  print(p); ggsave(file, p, width = 9, height = 6, bg = "white"); invisible(TRUE)
}
bubble_pair(SRC, TGT, "output/figs/07_6_bubble_compare.pdf", paste(SRC, "→", TGT))

# 換一組兩個條件都活著的細胞群——這不是補救，是把「這份資料到底能比什麼」講清楚。
# 只剩一群時畫的是自分泌（autocrine）：這份資料真正撐得起核心 vs 邊緣對比的，就只有 Macro/MG 對自己。
both <- Reduce(intersect, lapply(cc.list, cc_groups))
cat("兩個條件都存活的細胞群：", paste(both, collapse = "、"), "\n")
if (length(both) >= 1)
  bubble_pair(both, both, "output/figs/07_6_bubble_shared.pdf", "兩條件共同的細胞群")
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）
# 四位病人一致性（練習 7-4）：先看哪幾位病人兩個部位都跑得出來
pairs.ok <- Filter(function(x) all(paste0(x, c("_Tumor", "_Periphery")) %in% names(cc.all)),
                   unique(gbm4$patient))
cat("兩個部位都有 CellChat 結果的病人：", paste(pairs.ok, collapse = "、"), "\n")
# 對每位病人重複 §3，收集 rankNet 的顯著路徑，取交集
# rank.list <- lapply(patients, function(p) rankNet(mergeCellChat(list(cc.all[[paste0(p,"_Tumor")]], cc.all[[paste0(p,"_Periphery")]]), add.names = c("Core","Periphery")), mode = "comparison", do.stat = TRUE, return.data = TRUE)$signaling.contribution)

## ---- 4. liana-crosscheck --------------------------------------------- Q3 頁 65
if (requireNamespace("liana", quietly = TRUE)) {
  library(liana)
  obj <- subset(gbm4, patient == DEMO.p & tissue == DEMO.t); Idents(obj) <- "cc_label"
  li <- liana_wrap(obj, method = c("natmi", "connectome", "sca", "cellphonedb"), resource = "Consensus")
  li.agg <- liana_aggregate(li)                                              # 共識排名（aggregate_rank 越小越好）
  # dplyr:: 寫全名：org.Hs.eg.db／AnnotationDbi（06a 的 ORA 會載入）也有 select、filter，
  # 同一個 R session 先跑過 06a 再跑這裡會被蓋掉
  print(li.agg |> dplyr::filter(source == "Macro/MG", target == "Malignant") |>
        dplyr::select(source, target, ligand.complex, receptor.complex, aggregate_rank) |> head(10))
  # LIANA 自己也會丟掉 < 5 顆的群（見上面的訊息），所以群名一樣要先取交集再畫
  li.grp <- unique(c(li.agg$source, li.agg$target))
  li.tgt <- intersect(c("Malignant", "Vascular"), li.grp)
  if ("Macro/MG" %in% li.grp && length(li.tgt)) {
    p <- li.agg |> liana_dotplot(source_groups = "Macro/MG", target_groups = li.tgt, ntop = 15)
    print(p); ggsave("output/figs/07_liana_dotplot.pdf", p, width = 10, height = 6, bg = "white")
  } else cat("這個樣本裡 Macro/MG 或目標群被 LIANA 的 5 顆門檻擋掉了，跳過 dotplot\n")
  write.csv(li.agg |> dplyr::select(-starts_with("natmi"), -starts_with("connectome")) |> head(500),
            "output/tables/07_liana_top500.csv", row.names = FALSE)
  # 與 CellChat 的 bubble 對照：SPP1–CD44、MIF–CD74 兩邊都在前段，才寫進結果
## （這裡在解答版有一段參考答案；先自己跑出數字，再回去對照）
}
sessionInfo()

# =====================================================================
# ▶ 練習 7
#  7-1 把 subsetDB 改成全部三類（不加 search），互動數變多少？新增的顯著路徑主要是哪一類？
#  7-2 圖 4：這個樣本裡「既送又收」的樞紐是誰？換另一位病人，樞紐一樣嗎？
#  7-3 SPP1 路徑的 netAnalysis_contribution 前兩對是什麼？用 VlnPlot 確認：配體在來源群的表現比例 > 25% 嗎？
#  7-4 pairs.ok 列出的病人各做一次 rankNet 比較（這份資料是三位，BT_S6 的邊緣湊不出兩群）：
#      哪些路徑三位方向一致？只有一位病人顯著的路徑有幾條？三位一致跟一位顯著，能寫的話一不一樣？
#  7-5 LIANA 的共識前 10 對與 CellChat bubble 的前 10 對重疊幾對？不重疊的原因可能是什麼？
#      再做一件事：把前 10 對的「配體」一個一個查 UniProt 的 subcellular location，
#      有幾個真的是分泌型或單次穿膜的表面蛋白？細胞內的蛋白排進前十，代表什麼？
#  7-6 §1 的存活表：哪些「樣本 × 細胞群」被 MIN.CELLS 擋掉？把 MIN.CELLS 改成 10 重跑 §1–§3，
#      Macro/MG → Malignant 的並排 bubble 畫得出來了嗎？畫得出來的話，那張圖可以寫進論文嗎？
#      （提示：13 顆細胞估出來的 triMean 機率，換一個 seed 或少抽兩顆細胞，還會是同一個數字嗎？）
#      改 MIN.CELLS 會自動存成另一個檔名（..._min10.rds），不會誤用 20 那一輪的快取；
#      代價是那些樣本要重算一次，這是應該付的——門檻變了，網路就是不一樣的網路。
#  進階 用 NicheNet 反推：邊緣惡性細胞相對核心上調的基因（06 的 DE），最能被哪個配體解釋？與 CellChat 的結論一致嗎？
# =====================================================================
