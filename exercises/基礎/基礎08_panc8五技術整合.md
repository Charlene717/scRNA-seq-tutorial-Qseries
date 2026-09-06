# 基礎 08 · 整合的標準練習場：八個胰臟資料集、五種技術

**難度**：★★ ｜ **預估時間**：1–1.5 個工作天 ｜ **對應 Q 系列**：Q3・整合 ｜ **前置**：先完成基礎05、基礎06（胰臟型別與 marker 是本題的判準）；建議做過 B 系列 基礎11（跨技術比較的思路） ｜ **對應練習腳本**：`R/04_multipatient.R`

**原題**：B 系列 基礎17（本卡為 Q 系列改寫版，B 系列原檔未更動）

## 背景與研究主題

panc8 是整合方法論文最愛用的展示資料：八個人類胰臟資料集、五種技術（Smart-seq2、Fluidigm C1、CEL-seq、CEL-seq2、inDrops）合在一個 Seurat 物件裡，共 14,892 顆細胞。它好用的原因跟 B 系列 基礎11 的 pbmcsca 一樣——組織相同、技術不同，批次效應巨大而生物學已知；但它多給你一樣東西：**你自己就懂這個組織**。inDrops 與 CEL-seq2 的部分正是你在基礎05、基礎06 親手做過的 Baron 與 Muraro 資料，α/β/δ/PP 的 marker 你可以默寫。這讓你第一次有資格用「生物判準」而不是「目測混合度」來評整合：整合做得好不好，標準是**同型別跨批次合為一群、不同型別仍分得開**。本題把「不整合、CCA anchors、Harmony」三版並排，正面回答 Q3・整合 的問題。你的研究主題：

1. 不整合直接合併，八個資料集在 UMAP 上按什麼聚？同一個 β 細胞族群被切成幾塊？
2. CCA 與 Harmony 各自把型別對齊得如何？計算成本差多少？結果差在哪些型別上？
3. 整合有沒有「校正過頭」——把該保留的生物差異（稀有型別、donor 特徵）也抹掉？

## 資料集

- **SeuratData `panc8`**：14,892 cells，8 個資料集、5 種技術（Smart-seq2、Fluidigm C1、CEL-seq、CEL-seq2、inDrops）。
- 下載：R 內 `install.packages("SeuratData", repos="https://seurat.nygenome.org")` 後 `InstallData("panc8")`、`LoadData("panc8")`。
- 格式注意：技術別、資料集別與作者的型別註解都在 metadata，欄名以 `colnames(obj[[]])` 實際確認為準；作者註解欄留到對答案用。Smart-seq2 與 Fluidigm C1 是無 UMI 的全長技術，counts 尺度與液滴技術不同（B 系列 基礎11 討論過）。

## 任務

### 階段 A：不整合基準版（對應 Q3・整合）

1. metadata 偵查：`table()` 技術欄與資料集欄，確認 8 個資料集怎麼對應到 5 種技術。這裡藏著本題第一個決策——**整合範圍**：你要按「資料集」（8 批）還是按「技術」（5 批）整合？先寫下你的初步選擇與理由，做完階段 B 允許反悔。
2. QC 後直接 merge 跑標準流程到 UMAP（參數三件事照規矩），按技術、按資料集各出一張上色圖。用 GCG/INS/SST/PPY 快查：同一內分泌型別被技術切成幾塊？
3. 把作者註解欄拿出來做 celltype × tech 交叉表，記錄「不整合版」每個型別橫跨幾個 cluster——這張表是階段 C 評分的基線。

### 階段 B：兩條整合路線（對應 Q3・整合）

4. 路線一 CCA anchors：`SelectIntegrationFeatures` → `FindIntegrationAnchors` → `IntegrateData`，用你在第 1 步決定的整合範圍分批。記錄跑了多久、吃了多少記憶體（記憶體吃緊見提示 2）。
5. 路線二 Harmony：merge 物件的 PCA 空間上 `RunHarmony`（分組變數用同一個整合範圍）。同樣記錄成本。
6. 公平比較的紀律：兩版整合後的下游（dims、resolution）盡量用同一組參數，差異才可歸因於整合方法本身。各出 UMAP 兩張（按批次、按型別上色）。

### 階段 C：用胰臟知識當裁判（對應 Q3・整合、Q2・註解）

7. 三版（不整合／CCA／Harmony）並排評分：對每一版做 celltype × cluster 與 celltype × tech 交叉表，量化（a）同型別跨批次是否合為一群（混合）；（b）不同型別是否仍分開（解析度）。α/β/δ/PP、acinar、ductal、stellate 逐型別給結論。
8. 過度校正檢查：挑細胞數最少的 1–2 個型別（用作者註解找），看它們在兩版整合後是自成一群、還是被揉進大群裡；再挑一個你在基礎05 看過 donor 效應的型別，看整合是否連 donor 間可能的真實差異也一併抹平。寫 300–400 字比較短文：這份資料你最終選哪個方法、整合範圍選資料集還是技術、證據是哪幾張表。

## 繳交物

1. 可重跑的 R 專案（三版流程共用一份 QC 與參數設定）。
2. 圖：三版 UMAP（各按批次／型別上色，共 6 張）、關鍵型別的 marker DotPlot（英文標籤）。
3. 分析筆記：參數三件事＋整合範圍決策＋三版交叉表＋計算成本紀錄。
4. 方法比較短文（300–400 字）。

## 自我檢核點

- [ ] 整合範圍（8 批 vs 5 批）有明確決策與理由，且寫在筆記裡
- [ ] 不整合基準版有留檔，三版比較是在同一組下游參數下做的
- [ ] 每個主要型別都有「跨批次合了沒、跟別型別分得開嗎」的逐型別結論，不是整張圖一句「混得不錯」
- [ ] 有做過度校正檢查：至少一個稀有型別在整合後的下落有被追蹤
- [ ] CCA 與 Harmony 的計算成本（時間／記憶體）有實測數字
- [ ] 短文的結論每一句都指得出對應的表或圖

## 提示（卡關再看）

<details><summary>提示 1：metadata 偵查與拆批</summary>
`LoadData("panc8")` 後先 `colnames(obj[[]])`，再 `table()` 候選欄位；物件版本不同欄名可能不同，親眼確認。拆批用 `SplitObject(obj, split.by = "<你的批次欄>")`；Seurat v5 亦可用 `split()` 分 layer 的寫法，選一條路線並記錄版本。
</details>

<details><summary>提示 2：CCA 太慢或記憶體爆</summary>
八批兩兩找 anchor 很貴。省法：改用 `reduction = "rpca"`、或 reference-based（挑細胞數最多的 1–2 批當 reference）。你做的任何簡化都要寫進筆記——「為了跑得動而換方法」本身就是真實世界的參數決策。
</details>

<details><summary>提示 3：怎麼把「混合度」變成表</summary>
不必急著上 LISI/kBET：`table(obj$celltype, obj$seurat_clusters)` 與 `table(obj$celltype, obj$tech)` 兩張交叉表就能回答大半——好的整合裡，一個型別集中在一個 cluster、且該 cluster 內各技術都有貢獻。逐型別數「橫跨幾個 cluster」，三版一比就有排名。
</details>

## 進階挑戰

- 把整合範圍換成另一種（8 批 ↔ 5 批）重跑你選定的方法，看結論會不會翻盤——「範圍」跟「方法」哪個影響大？
- 用 LISI 或 kBET 把混合度做成數字，跟你的交叉表結論互相印證。
- 留一批出來當 query（如 Fluidigm C1），用其餘七批整合後的物件做 label transfer，體驗「整合當 reference」的用法。
- 整合的兩難在腫瘤資料會更兇險（病人間差異可能是生物訊號），見 B 系列 進階02；胰島的疾病版雙資料集互驗見 B 系列 進階54（T2D）；圖譜級整合與 mapping 的終點站是 B 系列 進階30。

## 參考文獻

- Stuart T, Butler A, et al. Comprehensive Integration of Single-Cell Data. *Cell* (2019).（panc8 合集與 CCA anchor 方法）
- Butler A, et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. *Nature Biotechnology* (2018).
- Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* (2019).
