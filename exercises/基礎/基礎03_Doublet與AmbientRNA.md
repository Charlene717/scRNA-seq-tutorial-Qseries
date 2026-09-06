# 基礎 03 · QC 的第二層：Doublet 與 Ambient RNA

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應 Q 系列**：Q2・QC（第二層 QC）、Q2・分群與註解 ｜ **對應練習腳本**：`R/02_cluster.R`、`R/03_annotate.R`、`R/01_qc.R`

## 背景與研究主題

基礎01 的 QC 是第一層：三個閾值切掉爛細胞。但有兩種髒東西是閾值切不掉的——**doublet**（兩顆細胞裝進同一個液滴，變成一個「嵌合轉錄組」，會偽裝成過渡態或新型別）與 **ambient RNA**（破掉細胞漏出的 RNA 漂在液體裡，混進每個液滴，讓 marker「到處都亮一點」）。這一題用 10x 官方的 10k PBMC 練 QC 的第二層。注意一個細節：這份資料**目標是 10k 顆，實際 detected 11,769 顆**——多出來的不是白撿的，上樣越多 doublet 率越高，這個「超收」本身就是本題的討論素材。你的研究主題：

1. 這份資料裡有多少 doublet？它們在 UMAP 上藏在哪、長什麼樣？
2. ambient RNA 污染有多重？哪些 marker 因此「越界發光」？
3. 清完 doublet 與 ambient 之後，分群與 marker 圖變了什麼？哪些結論其實是髒訊號撐起來的？

## 資料集

- **10x 10k PBMCs from a Healthy Donor (v3 chemistry)**，11,769 cells detected（目標 10k）。
- 下載：`https://www.10xgenomics.com/datasets` 搜尋「10k PBMCs from a Healthy Donor (v3 chemistry)」，填 email 後**同時下載 filtered 與 raw 的 h5**（filtered/raw feature-barcode matrix）——raw 這次不是「Q2・QC 之前的世界」，它是估 ambient RNA 的原料。
- 格式注意：兩個 h5 都用 `Read10X_h5()` 讀；raw 含幾十萬個空液滴 barcode，屬正常，別對它做 QC。

## 任務

### 階段 A：基準線（對應 Q2 全流程）

1. 用 filtered 矩陣照基礎01 的骨架跑到註解完成（QC 三件事、標準流程、marker 註解）。這版叫 **baseline**，凍結留檔，之後所有比較都以它為對照。
2. 記下 baseline 的：細胞數、cluster 數、各型別比例、每群 top markers。順手觀察：有沒有哪個群的 marker 組合「不合常理」（例如同時亮 T 與單核球的 marker）？先記下嫌疑名單，不動手。

### 階段 B：doublet 偵測（對應 Q2・QC）

3. 用 `scDblFinder` 跑 doublet 偵測（預設參數起步；它假設的 doublet 率與 10x 的「每千顆約增加 ~0.8%」經驗律有關——查它文件裡 `dbr` 的預設邏輯，寫進筆記）。把 doublet 分類與 score 掛回 Seurat metadata。
4. 交叉檢視，別盲信工具：(a) doublet 在 UMAP 的分布——是散在各群邊緣還是聚成獨立小群？(b) 與 nCount 的關係——doublet 群的 nCount 是否偏高？boxplot 佐證；(c) 與階段 A 嫌疑名單對照——同時亮兩個 lineage marker 的群被抓到了嗎？
5. 回答討論題：目標 10k、實收 11,769，超收對 doublet 率意味著什麼？如果這是你自己的實驗，上樣量你會怎麼決定？（Q2・QC：QC 問題有一半其實是實驗設計問題。）

### 階段 C：ambient RNA 與總清算（對應 Q2・QC、Q2・分群與註解）

6. 用 **SoupX**（需要 raw + filtered + baseline 的 cluster）或 **decontX**（用 raw 當 background）估 ambient 污染分數並產出校正後矩陣。記錄整體污染率估計值。
7. 找證據：挑 2–3 個「該專一卻不專一」的 marker（PBMC 常見嫌疑：LYZ、HBB、PPBP 這類高表達基因的越界瀰漫），畫校正前後的 FeaturePlot/VlnPlot 並排。
8. 總清算：移除 doublet＋用校正後矩陣重跑分群註解，與 baseline 比較——cluster 數、型別比例、你在階段 A 存疑的群還在嗎？每群 top marker 變乾淨了嗎？寫 300–400 字筆記：「QC 第二層」改寫了 baseline 的哪些結論？哪些又幾乎沒變（這也重要——說明訊號穩健）？

## 繳交物

1. 可重跑的 R 專案（baseline 與 cleaned 兩條流程都可重跑）。
2. 圖：doublet score/分類 UMAP、doublet vs nCount boxplot、ambient 校正前後 marker 並排圖、baseline vs cleaned 的 UMAP 與比例比較（英文標籤）。
3. 分析筆記：參數三件事＋doublet 交叉檢視紀錄＋污染率估計。
4. 「QC 第二層」總結筆記（300–400 字）。

## 自我檢核點

- [ ] baseline 有凍結留檔，所有比較都是同一把尺
- [ ] doublet 結果經過至少兩種獨立證據交叉檢視（UMAP 分布、nCount、marker 共表現），不是工具說了算
- [ ] 能解釋「超收 vs doublet 率」的關係，並連回實驗設計
- [ ] ambient 校正有 before/after 的 marker 圖為證，說得出污染率數字與它的意義
- [ ] 總清算筆記同時寫了「被改寫的結論」與「沒變的結論」

## 提示（卡關再看）

<details><summary>提示 1：scDblFinder 上手</summary>
`sce <- scDblFinder(as.SingleCellExperiment(obj))` 後把 `sce$scDblFinder.class` 與 `sce$scDblFinder.score` 塞回 `obj@meta.data`。它靠模擬 doublet 訓練分類器，同型別互撞的 homotypic doublet 本來就難抓——所以交叉檢視才是必修，不是加分題。
</details>

<details><summary>提示 2：SoupX 三件原料</summary>
`SoupChannel(tod = raw, toc = filtered)` → `setClusters()`（用 baseline 的 cluster）→ `autoEstCont()` → `adjustCounts()`。autoEstCont 需要 cluster 資訊才能找「某群該是零卻不是零」的基因來估污染——這就是為什麼要先有 baseline 才能清 ambient。
</details>

<details><summary>提示 3：清完差異不大？</summary>
健康 PBMC 本來就是相對乾淨的樣本，污染率個位數％、doublet 抓到幾百顆是常態。這題教的是流程與判讀，不保證戲劇性翻盤——把「差多少」誠實量化出來就是正確答案。真正戲劇性的場面在腫瘤與核抽取資料裡等你。
</details>

## 進階挑戰

- 用 `DoubletFinder` 再跑一次，與 scDblFinder 的結果做交叉表——兩個工具的分歧細胞是誰？工具間一致性本身就是可靠度的量尺。
- 對 raw 矩陣畫 barcode-rank（knee）圖，理解 Cell Ranger 的 cell calling 在切哪裡，「11,769」這個數字就是這條曲線切出來的。
- 這一層 QC 是所有進階題的隱藏前置：每張進階卡的資料都值得先問「doublet 清了嗎？ambient 估了嗎？」——尤其進階05/14/15 這類 per-sample 腫瘤資料（每個樣本要分開跑），以及任何你打算拿去發表的分析：審稿人現在會問這兩件事。

## 參考文獻

- 10x Genomics Datasets：10k PBMCs from a Healthy Donor (v3 chemistry)（官方資料集頁；報告引用時寫明資料集名稱與下載日期）。
- Germain PL, et al. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* (2021).
- Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience* (2020).
- Yang S, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. *Genome Biology* (2020).
