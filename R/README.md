# 練習腳本

Q 系列課程的十二支 R 腳本，從讀檔一路做到反卷積與存活分析。每支都可以獨立閱讀：
開頭註明對應的投影片頁碼、輸入輸出與預估時間，結尾附 4–6 題練習。

## 兩個版本

| 資料夾 | 是什麼 |
|---|---|
| `R/`（這裡） | **完整版**，可以從頭跑到尾，等於解答 |
| [`R/練習版/`](練習版) | 把影片裡講過「為什麼」的關鍵參數挖空成 `____`，空格上方有 `## TODO ▶` 提示與投影片頁碼 |

建議先自己填 `練習版/`，跑不通再對照完整版。

## 執行順序

腳本有先後依賴：後面的讀前面存下的 `.rds`，不能跳著跑。

| 腳本 | 做什麼 | 投影片 | 時間 |
|---|---|---|---|
| [`00_setup.R`](00_setup.R) | 裝套件、下載資料與基因座標檔、建資料夾 | Q2 P7–8 | 10–20 分 |
| [`01_qc.R`](01_qc.R) | SoupX（選做）→ 讀檔 → QC 三指標 → MAD 閾值 → DoubletFinder | Q2 P8–22 | 3–5 分 |
| [`02_cluster.R`](02_cluster.R) | 前處理與週期分數 → nPC → 解析度掃描 → 穩定性檢查 → doublet 群診斷 | Q2 P25–36 | 3–5 分 |
| [`03_annotate.R`](03_annotate.R) | marker 面板 → SingleR → 掛名字 → 免疫亞群 → Neftel 狀態分數 → 交付 | Q2 P38–64 | 10–15 分 |
| [`04_multipatient.R`](04_multipatient.R) | 載入 GSE84465 → 未整合基線 → CCA + Harmony → 正負對照 → LISI | Q3 P8–17 | 5–10 分 |
| [`05_infercnv.R`](05_infercnv.R) | inferCNV → CNV 分數與相關 → 三角驗證 → 與作者標籤對答案 | Q3 P21–28 | 10–30 分 |
| [`06a_pseudobulk_gsea.R`](06a_pseudobulk_gsea.R) | 組成分析 → 每型別 pseudobulk + 配對 DESeq2 → 火山圖 → GSEA / ORA | Q3 P29–51 | 5–10 分 |
| [`06b_cell_level_de.R`](06b_cell_level_de.R) | 不能 pseudobulk 時的備案：MAST + 病人共變量、逐病人一致性、標籤置換 | Q3 P53–54 | 3–5 分 |
| [`07_cellchat.R`](07_cellchat.R) | 每樣本各跑 CellChat → 六種圖 → 兩條件比較 → LIANA 交叉驗證 | Q3 P58–71 | 每樣本 5–15 分 |
| [`08_trajectory.R`](08_trajectory.R) | 單一病人惡性細胞的軌跡：Slingshot → tradeSeq → Monocle 比較 | Q3 P72–75 | 5–10 分 |
| [`09_activity.R`](09_activity.R) | PROGENy 路徑活性 → 以病人為單位配對比較 → SCENIC（選配） | Q3 P76–78 | 2 分 |
| [`10_deconv_survival.R`](10_deconv_survival.R) | MuSiC 反卷積 TCGA-GBM → KM / Cox 存活分析 | Q3 P79–82 | 15 分（含下載） |

用法、作業格式與**常見錯誤對照表**見 [`練習手冊.md`](練習手冊.md)。

## 開始之前

- 新建一個 RStudio Project，把 `R/` 複製進去，從 `00_setup.R` 開始。
- 全程不要用 `setwd()`，所有路徑相對於專案根目錄。
- `05_infercnv.R` 需要**系統層級**的 JAGS 4.x（rjags 只是介面，R 裝不了 JAGS 本體；裝完要重開 R）。
  還沒裝 JAGS 之前，06 之後的腳本可以先用作者標籤當替代惡性標籤測試（見 `06a` §0）。
- 全程 `set.seed(1234)`。同一套環境重跑結果一致。

資料檔與分析輸出不隨倉庫發布（見 `.gitignore`），`00_setup.R` 會自動下載前兩份資料。

課程說明、投影片與影片連結在[倉庫首頁](../README.md)。
