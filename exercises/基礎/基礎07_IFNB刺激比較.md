# 基礎 07 · 刺激 vs 對照：狼瘡 PBMC 的 IFN-β 反應

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應 Q 系列**：Q3・整合與差異表達 ｜ **對應練習腳本**：`R/04_multipatient.R`、`R/06a_pseudobulk_gsea.R`

## 背景與研究主題

Kang 等人（2018）把 8 位狼瘡病人的 PBMC 分成兩半：一半不處理（ctrl），一半用 IFN-β 刺激（stim）。這是「條件比較」設計的教學經典：同一批病人、同一種細胞、兩個條件——問題從「有哪些細胞」變成「**同一種細胞在刺激後改變了什麼**」。這一題也是你第一次正面遭遇 Q3・差異表達 的核心警告：統計單位是病人，不是細胞。你的研究主題：

1. ctrl 與 stim 的細胞需要整合嗎？整合後同型別跨條件對得齊嗎？
2. 各免疫型別對 IFN-β 的反應一樣大嗎？哪個型別反應最強？
3. 同一個 DE 問題，cell-level 與 pseudobulk 的答案差多少？為什麼該信後者？

## 資料集

- **[GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583)**（Kang）；本題直接用 [SeuratData](https://github.com/satijalab/seurat-data) 的 **`ifnb`** 物件：13,999 cells，8 位病人的 ctrl + IFN-β 兩條件（10x 平台，demuxlet 多工分樣）。
- 下載：R 內 `install.packages("SeuratData", repos="https://seurat.nygenome.org")` 後 `InstallData("ifnb")`、`LoadData("ifnb")`。
- 注意：(1) 物件的 `stim` 欄是條件；病人編號在 metadata（不同版本欄名可能是 `donor`/`ind` 之類，先 `colnames(obj[[]])` 確認）——pseudobulk 沒有它就做不成。(2) 物件是未處理的 counts，QC 與標準流程照常要走。

## 任務

### 階段 A：整合 ctrl 與 stim（對應 Q3・整合）

1. QC 後把資料按條件 split，跑 Seurat anchors 或 Harmony 整合（整合方法、維度數等參數：三件事照規矩）。
2. 出兩張 UMAP：按條件上色（檢查混合）、按分群上色。再用 marker 註解主要型別（PBMC 的 marker 你在基礎01 練過了）。
3. 回答 Q3・整合 的必答題：這份資料**為什麼適合整合**？（同一批病人、實驗性刺激、預期同型別跨條件對應——對照腫瘤資料想一想差在哪——那裡的病人間差異可能是生物訊號，不見得該修掉。）

### 階段 B：組成與 DE——統計單位！（對應 Q3・差異表達）

4. 各型別在 ctrl/stim 的細胞比例表與長條圖：刺激有沒有改變組成？（比例的統計檢定注意 Q3・差異表達 講過的組成資料陷阱，描述性呈現即可。）
5. 挑單核球（CD14+ Mono）做 ctrl vs stim 的 DE，**做兩遍**：
   - (a) cell-level：`FindMarkers` 直接以細胞為單位。
   - (b) pseudobulk：以「病人 × 條件」`AggregateExpression` 加總成 16 個樣本，用 DESeq2 或 edgeR 做配對比較。
6. 比較 (a)(b)：各自的顯著基因數、p 值分布、logFC 一致性（散點圖）。回答：cell-level 的 p 值為什麼膨脹？（同一位病人的上千個細胞不是獨立樣本——Q3・差異表達 的核心。）

### 階段 C：富集分析與生物學解讀（對應 Q3・富集）

7. 用 pseudobulk 的顯著基因跑富集（clusterProfiler 的 GO/或 MSigDB hallmark 的 GSEA）。interferon response 相關 pathway 應該明顯跳出——如果沒有，先回頭查 DE 方向與基因 ID 轉換。
8. 對至少兩個型別重複步驟 5(b)，比較 IFN 反應強度（顯著基因數或 hallmark IFN score）。寫 200–400 字小結：哪個型別反應最強？與文獻預期一致嗎？

## 繳交物

1. 可重跑的 R 專案（renv + `set.seed(1234)`）。
2. 圖：整合前後/條件上色 UMAP、組成長條圖、cell-level vs pseudobulk 的 p 值與 logFC 對比圖、富集結果圖（英文標籤）。
3. 分析筆記：參數三件事＋各階段回答。
4. 「統計單位」說明短文（200–400 字）：用你自己的 (a)(b) 結果解釋給沒修過 Q3・差異表達 的同學聽。

## 自我檢核點

- [ ] 能說出這份資料適合整合的理由，而不是「教學都有整合所以整合」
- [ ] 整合後同型別 ctrl/stim 混合良好，且有圖為證
- [ ] pseudobulk 的樣本數是 16（8 病人 × 2 條件），設計矩陣有配對（病人）項
- [ ] 能用自己的圖說明 cell-level p 值膨脹的現象與原因
- [ ] 富集結果的解讀有連回生物學（IFN-β 刺激），不是貼一張 dotplot 了事

## 提示（卡關再看）

<details><summary>提示 1：pseudobulk 實作</summary>
`AggregateExpression(obj, assays="RNA", group.by=c("celltype","donor_id","stim"))` 取 counts 後，對目標型別建 colData（donor、condition），DESeq2 設計式用 `~ donor + condition`——donor 當配對項，這就是「以病人為統計單位」的具體寫法。
</details>

<details><summary>提示 2：IFN 反應太強害你對不齊型別</summary>
stim 組整體轉錄組被 IFN 拉動，未整合時同型別會按條件分開——這正是本題需要整合的原因。註解時以整合後的 cluster 為單位，用不受 IFN 影響的 lineage marker（CD3D、CD14、MS4A1…）判型別。
</details>

<details><summary>提示 3：富集沒跳出 interferon</summary>
三個常見死因：基因 ID 沒轉成 ENTREZ/SYMBOL 一致格式、拿了不顯著基因當前景、logFC 方向弄反（stim vs ctrl 還是 ctrl vs stim）。逐一排查。
</details>

## 進階挑戰

- 對「顯著基因數最多」與「最少」的型別各畫 top 基因的 per-patient pseudobulk heatmap，檢查訊號是所有病人一致還是被一兩位病人拉動——這是 Q3・差異表達「看見病人」的延伸。
- 用 `AddModuleScore` 算 hallmark interferon response score，畫各型別 ctrl/stim 的 score 分布，與 DE 結果互相印證。
- 這套「條件比較 + pseudobulk」的思路，在進階14（黑色素瘤免疫治療 responder vs non-responder）會搬進臨床腫瘤場景。

## 參考文獻

- Kang HM, et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nature Biotechnology* (2018). GEO: GSE96583。
- SeuratData `ifnb`（Seurat 官方整合教學同款資料）。
