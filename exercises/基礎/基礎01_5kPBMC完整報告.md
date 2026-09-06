# 基礎 01 · 你的第一份完整報告：5k PBMC

**難度**：★ ｜ **預估時間**：0.5–1 個工作天 ｜ **對應 Q 系列**：Q1、Q2 全流程 ｜ **對應練習腳本**：`R/01_qc.R`、`R/02_cluster.R`、`R/03_annotate.R`

## 背景與研究主題

PBMC（周邊血單核細胞）是單細胞領域的「模式資料」：細胞型別界線清楚、marker 教科書化、任何教學都拿它開場。這一題不是要你發現新生物學，而是把 Q2 全流程 的完整流程**獨立**走一遍，並且把結果寫成一份像論文的報告——這份報告會成為你之後所有題目的基準模板。你的研究主題：

1. 這份 5k PBMC 裡有哪些主要免疫細胞型別？各佔多少比例？
2. 你的每一個參數選擇（QC 閾值、nPC、resolution…），換成別人來問，你答得出「數字、理由、紀錄」三件事嗎？
3. 一份「別人能照著重現」的方法段落，到底要寫到多細？

## 資料集

- **10x 5k Human PBMCs（3' v3.1）**，~5,000 cells，健康捐贈者周邊血。
- 下載：到 `https://www.10xgenomics.com/datasets` 搜尋「5k Peripheral blood mononuclear cells」，填 email（免費）後下載 **Filtered feature-barcode matrix**（HDF5 或 tar.gz 皆可）。
- 注意：抓 *filtered* 不是 *raw*——raw 矩陣含大量空液滴，是 Q2・QC 之前的世界。h5 用 `Read10X_h5()`，tar.gz 解開後用 `Read10X()`。

## 任務

### 階段 A：跑通（對應 Q2 全流程）

1. 建 RStudio Project（renv、`set.seed(1234)`、相對路徑），讀入矩陣，建 Seurat 物件。
2. 按課程順序一路跑到底，先不糾結參數（用課程示範值即可）：QC（Q2・QC）→ NormalizeData（Q2・前處理）→ FindVariableFeatures + ScaleData（Q2・前處理）→ RunPCA（Q2・PCA）→ FindNeighbors/FindClusters + RunUMAP（Q2・分群）。
3. 出一張 UMAP（依 cluster 上色）。此刻的目標只有一個：**流程從頭到尾不報錯**。

### 階段 B：每個參數補三件事（對應 Q2 前半流程）

4. 回頭檢視階段 A 用到的每一個自訂參數，逐一補上「數字、理由、紀錄」：
   - QC 三閾值（nFeature、nCount、percent.mt）：畫 violin/散點圖，說明閾值切在分布的哪裡、為什麼（Q2・QC）。
   - HVG 數量（預設 2000）：換 1000、3000 看 elbow 前段 PC 是否穩定（Q2・前處理）。
   - nPC：ElbowPlot 加上你選擇的理由——「拐點在哪」是理由，「教學用 10 所以我用 10」不是（Q2・PCA）。
   - resolution：至少試 0.4 / 0.8 / 1.2，用 `clustree` 或並排 UMAP 說明你選哪個、為什麼（Q2・分群）。
5. 若任何參數改動後分群明顯變化，記錄「什麼變了、你如何取捨」——這段筆記之後直接進報告。

### 階段 C：註解與完整報告（對應 Q2・註解、Q1）

6. 用經典 marker 註解每個 cluster（T：CD3D/CD3E；CD4 vs CD8：CD4/CD8A；B：MS4A1/CD79A；NK：NKG7/GNLY；單核球：CD14/FCGR3A/LYZ；DC：FCER1A；血小板：PPBP），出 DotPlot + 註解後 UMAP，並做一張型別比例表。
7. 寫成一份完整分析報告（md 或 Rmd/Quarto 輸出皆可），結構模仿論文：
   - **Methods**：從下載到註解，每步的軟體版本、函數、參數與理由——標準是「另一個修完 Q 系列的人能照著重現」。
   - **Results**：圖（英文標籤）＋每張圖 2–3 句描述。
   - **Data availability**：資料來源與下載日期（格式見資料集總覽）。

## 繳交物

1. 可重跑的 R 專案（renv + `set.seed(1234)`，相對路徑）。
2. 圖：QC 分布圖、ElbowPlot、resolution 比較圖、marker DotPlot、註解後 UMAP（英文標籤）。
3. 參數筆記：每個自訂參數的「數字、理由、紀錄」三件事（中文可）。
4. 完整分析報告一份（Methods + Results + Data availability）。

## 自我檢核點

- [ ] 專案能在乾淨環境重跑到底，圖與報告一致
- [ ] 每個 QC 閾值都能指著自己畫的分布圖說出理由，不是抄課程數字
- [ ] resolution 的選擇有比較過至少三個值，且說得出「為什麼不是另外兩個」
- [ ] 每個 cluster 的註解都有 marker 證據（DotPlot 上指得出來），沒有「看位置猜的」
- [ ] Methods 段落拿給同學，對方不用問你任何問題就能重現

## 提示（卡關再看）

<details><summary>提示 1：percent.mt 怎麼訂</summary>
`PercentageFeatureSet(obj, pattern = "^MT-")` 後看 violin 分布。PBMC 常見切在 5–10%，但理由必須來自你這份資料的分布（例如「主峰在 3% 以下，10% 之後是稀疏長尾」），而不是「大家都用 5%」。
</details>

<details><summary>提示 2：cluster 對不上 marker</summary>
先出一張 `FeaturePlot` 把 CD3D、MS4A1、CD14、NKG7 四個大類 marker 攤開看，確定大陸塊的歸屬，再處理小群。若某群什麼 marker 都不亮，回 Q2・QC 想想：會不會是低品質細胞或 doublet 漏網？
</details>

<details><summary>提示 3：Methods 寫多細</summary>
找一篇你領域的論文，讀它的 scRNA-seq Methods 段——注意它寫了版本號、參數值與過濾標準。你的標準只多一條：每個參數多一句理由。
</details>

## 進階挑戰

- 用 `DoubletFinder` 或 `scDblFinder` 加一步 doublet 偵測，比較移除前後的分群，把這步也寫進 Methods。
- 把同一流程對 3k PBMC（10x 教學經典款）重跑一次，比較兩份資料的型別比例——體會「同組織不同批次」的差異有多大。
- 這套單樣本完整流程的下一站是腫瘤檢體：進階01（黑色素瘤）、進階02（頭頸癌）用同一套骨架面對腫瘤，會遇到惡性細胞這個 PBMC 沒有的難題。

## 參考文獻

- 10x Genomics Datasets：5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, 3' v3.1（官方資料集頁，無對應論文；報告引用時寫明資料集名稱與下載日期）。
- Hao Y, et al. Integrated analysis of multimodal single-cell data. *Cell* (2021).（Seurat 流程引用）
