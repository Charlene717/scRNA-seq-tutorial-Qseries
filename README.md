# scRNA-seq Q 系列 · 快速上手篇

> 三集課程、約 2.5 小時、十二支可實跑的 R 腳本，外加 23 題課後實作練習——從「單細胞 RNA 定序在看什麼」到多樣本整合、惡性細胞判定、差異表達與下游分析，完成一份完整的 scRNA-seq 分析。

![code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)
![content: CC BY--NC 4.0](https://img.shields.io/badge/content-CC%20BY--NC%204.0-lightgrey.svg)
![R ≥ 4.3](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3.svg)
![Seurat v5](https://img.shields.io/badge/Seurat-v5-1a6b5a.svg)

這是一門給趕時間的人的 scRNA-seq 快速上手課，內容涵蓋**通用的分析章節**（QC、分群、細胞註釋、多樣本整合、差異表達與富集、細胞通訊、軌跡、反卷積）與**腫瘤特有的章節**（inferCNV 判定惡性、惡性細胞按病人分群的整合兩難、CNV 三角驗證），以膠質母細胞瘤（GBM）的公開資料作為貫穿範例。課程的核心主張是：每個參數都是一個假設，要有數字、有理由、有紀錄；本倉庫的十二支練習腳本全部可以從頭跑到尾（已於 Windows、R 4.4.1、Seurat 5.3.0 實跑驗證），`set.seed(1234)`，結果可重現。

## 課程一覽

| # | 標題 | 投影片 | 練習腳本 |
|---|---|---|---|
| Q1 | 觀念篇：單細胞 RNA 定序在看什麼、怎麼運作、怎麼讀圖 | [中文](slides/Q1_觀念篇_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/Q1_Concepts_EN.pdf) | — |
| Q2 | 實作篇 I：單一樣本標準流程——從原始矩陣到細胞註釋 | [中文](slides/Q2_實作篇I_單一樣本標準流程_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/Q2_Single_Sample_Pipeline_EN.pdf) | [00–03](R/) |
| Q3 | 實作篇 II：多樣本整合與下游分析 | [中文](slides/Q3_實作篇II_多樣本分析_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/Q3_Multi_Sample_Analysis_EN.pdf) | [04–10](R/) |

影片連結：錄製中，上線後補。三集投影片共 228 頁，**中英雙語各一套**（PDF）。

**九句話骨架**——整門課要帶走的東西：

1. 腫瘤是生態系，單細胞把它拆開；分析在 PC 空間，UMAP 只是投影；惡性看基因體，病人是統計單位。（Q1）
2. 每個參數是一個假設：有數字、有理由、有紀錄；註釋先確定譜系，惡性與否留待多重證據共同判定。（Q2）
3. 整合的目的是對齊各樣本的共同細胞型別；惡性細胞帶病人特異的基因體變異，宜逐病人分析；差異表達以病人為統計單位；下游分析各有其適用的資料條件。（Q3）

## 適合誰

| 你想要 | 看 | 做 | 動手時間預估 |
|---|---|---|---|
| 看懂單細胞論文、跟生資合作者對話 | Q1 | — | —（只看影片） |
| 自己分析一份樣本（QC → 分群 → 註釋） | Q1 + Q2 | 腳本 00–03 | 1–2 天 |
| 多樣本比較、判定惡性、下游分析與發表 | Q1–Q3 | 腳本 00–10 | 3–5 天 |
| 把學到的東西練成一份研究 | Q1–Q3 | 腳本 00–10 + [`exercises/`](exercises/) | 挑一題進階練習題做完，2–3 週 |

**時間預估怎麼看**：這裡估的是動手的時間，不含看影片。估計假設你會停下來想每個參數的理由、
把「數字、理由、紀錄」寫下來——只求「跑完不報錯」會快很多，但那不是這門課要教的。
腳本 00–10 光是跑完就要約 2–4 小時機器時間（inferCNV 與 CellChat 佔大半），套件安裝與資料下載的等待另計。
進階練習題一題就是一個完整的再分析專案：讀懂資料、跑完課程教過的全部下游分析，再把訊號推到臨床驗證，抓 7–14 個工作天，第一次接觸單細胞的人會更久一些。十四題是讓你挑一個自己想做的題目，不是要你全部做完。若要做到能投稿，文獻回顧、外部資料集驗證與寫作是另一段路。

前置需求：會開 RStudio、跑過幾行 R。不需要單細胞經驗。

## 各集內容

### Q1 · 觀念篇（47 頁，無腳本）

從一位 58 歲 GBM 病人的切片出發，建立三個心智模型：腫瘤（組織）是什麼、資料從哪裡來、分析流程在做什麼。內容包括：Bulk 定序看到什麼、漏掉什麼；從組織到矩陣（解離、四種平台、GEM、barcode/UMI、稀疏與兩種零、實驗設計三個數字）；九步流程各在做什麼、PC 空間與 UMAP 的關係；腫瘤資料的兩個特有現象（惡性細胞按病人分群、惡性與否需要基因體層級的證據）；怎麼讀圖與讀論文（UMAP 能說與不能說的、dotplot、FeaturePlot、四種常見誤讀、讀論文五問）；三個常見錯誤與三句總結。

### Q2 · 實作篇 I：單一樣本標準流程——從原始矩陣到細胞註釋（80 頁，腳本 00–03）

資料：10x 官方 GBM 5k（一位病人，Chromium 3' v3，5,604 顆細胞）。從讀檔做到「每群有名字、惡性細胞有狀態分數」，每個參數答得出為什麼。內容包括：環境安裝與資料初探（含環境 RNA 與 SoupX）；QC 三指標與 MAD 動態閾值、DoubletFinder 抓 doublet（並與 scDblFinder 對照）；前處理四行與細胞週期分數、LogNormalize vs SCT、nPC 怎麼選；分群解析度掃描與 clustree、群穩定性三檢查；細胞註釋的完整工作流（門牌基因、marker 面板、FindAllMarkers、SingleR、層級式註釋與亞群重跑、免疫亞群面板、Neftel 四狀態三種打分數算法、三層命名規範、常見誤註釋）；方法段模板與交付清單。

### Q3 · 實作篇 II：多樣本整合與下游分析（101 頁，腳本 04–10）

資料：GSE84465（Darmanis et al. 2017；4 位 GBM 病人 × 腫瘤核心/浸潤邊緣，Smart-seq2，3,589 顆，含作者的細胞型別標籤——可以「對答案」）＋ TCGA-GBM Bulk。內容包括：多病人整合的兩難與方法地圖（Seurat CCA 為預設、Harmony 作比較、正負對照與 LISI）；inferCNV 判定惡性（CNV 分數與相關、四象限、三角驗證、與作者標籤對答案）；組成分析（propeller）；差異表達正確做法（pseudoreplication 為什麼錯、每型別各自 pseudobulk + 配對 DESeq2、三種 log2FC 的用途、火山圖讀法、不能 pseudobulk 時的 cell-level 備案與三道防線）；富集分析（ORA vs GSEA、Hallmark/GO/KEGG、NES 熱圖與 dotplot 讀法、五個常見錯）；細胞通訊（CellChat 六種圖的畫法與讀法、兩條件比較、LIANA 交叉驗證）；軌跡分析（Slingshot + tradeSeq、Monocle3/Monocle2 比較、手動選起點）；路徑與 TF 活性（decoupleR/PROGENy、SCENIC）；反卷積與存活分析（MuSiC + TCGA、KM/Cox）。

## 練習腳本（R/）

每支腳本可獨立閱讀：開頭註明對應投影片頁碼、輸入輸出與預估時間，結尾附 4–6 題練習（含進階題）。兩個版本：

- [`R/`](R/) — **完整版（解答）**，從頭跑到尾。
- [`R/練習版/`](R/練習版/) — 投影片裡講過「為什麼」的關鍵參數挖空成 `____`，空格上方有 `## TODO ▶` 提示與投影片頁碼；填完再跑，跑不通對照完整版。

| 腳本 | 做什麼 | 投影片&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 運行時間&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |
|---|---|---|---|
| [`00_setup.R`](R/00_setup.R) | 裝套件（CRAN + Bioconductor + GitHub）、下載兩份資料與基因座標檔、建資料夾 | Q2<br>P7–8 | 10–20&nbsp;分 |
| [`01_qc.R`](R/01_qc.R) | （選做 SoupX）→ Read10X → 初探 → QC 三指標 → MAD 三條線 → DoubletFinder | Q2<br>P8–22 | 3–5&nbsp;分 |
| [`02_cluster.R`](R/02_cluster.R) | 前處理四行 + 週期分數 → PC1 檢查 → nPC → 解析度掃描 + clustree → 穩定性三檢查 | Q2<br>P25–36 | 3–5&nbsp;分 |
| [`03_annotate.R`](R/03_annotate.R) | 門牌基因 → marker 面板 → FindAllMarkers → SingleR → 掛名字 → 免疫亞群重跑 → Neftel 分數 → 品質檢查 → 交付 | Q2<br>P38–64 | 10–15&nbsp;分 |
| [`04_multipatient.R`](R/04_multipatient.R) | 載入 GSE84465 counts + metadata → 未整合基線 → Seurat CCA（預設）+ Harmony（比較）→ 正負對照 → LISI | Q3<br>P8–17 | 5–10&nbsp;分 |
| [`05_infercnv.R`](R/05_infercnv.R) | inferCNV（參考 = 免疫 + 寡樹突）→ CNV 分數與相關 → 三角驗證 → 與作者標籤對答案 | Q3<br>P21–28 | 10–30&nbsp;分 |
| [`06a_pseudobulk_gsea.R`](R/06a_pseudobulk_gsea.R) | propeller 組成 → 每型別 pseudobulk + 配對 DESeq2（ashr 收縮）→ 火山圖 → GSEA（Hallmark/GO/KEGG）+ ORA → NES 熱圖 | Q3<br>P29–51 | 5–10&nbsp;分 |
| [`06b_cell_level_de.R`](R/06b_cell_level_de.R) | 不能 pseudobulk 時的備案：Wilcoxon 基準 → MAST + 病人共變量 → 逐病人一致性 → 標籤置換 → 標註限制交付 | Q3<br>P53–54 | 3–5&nbsp;分 |
| [`07_cellchat.R`](R/07_cellchat.R) | 每樣本各跑 CellChat → 六種圖 → mergeCellChat 兩條件比較 → LIANA 共識交叉驗證 | Q3<br>P58–71 | 每樣本&nbsp;5–15&nbsp;分 |
| [`08_trajectory.R`](R/08_trajectory.R) | 單一病人惡性細胞的軌跡：Slingshot → tradeSeq → Monocle3/Monocle2 比較（含手動選起點）→ 穩健性 | Q3<br>P72–75 | 5–10&nbsp;分 |
| [`09_activity.R`](R/09_activity.R) | decoupleR/PROGENy 路徑活性 → 以病人為單位的配對比較 → SCENIC（選配） | Q3<br>P76–78 | 2&nbsp;分 |
| [`10_deconv_survival.R`](R/10_deconv_survival.R) | MuSiC 以單細胞為參考反卷積 TCGA-GBM → KM / Cox 存活分析 | Q3<br>P79–82 | 15&nbsp;分（含下載） |

使用說明、作業格式與**常見錯誤對照表**見 [`R/練習手冊.md`](R/練習手冊.md)。

## 資料

三份公開資料，前兩份由 `00_setup.R` 自動下載（失敗時腳本內附手動下載位置）：

| 資料 | 內容 | 用在&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 大小&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |
|---|---|---|---|
| 10x GBM 5k<br>（[10x Datasets](https://www.10xgenomics.com/datasets)） | 一位 GBM 病人，10x 3' v3，5,604 顆 | Q2：01–03 | ≈&nbsp;30&nbsp;MB |
| [GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)<br>（Darmanis et al. 2017） | 4 病人 × 核心/邊緣，Smart-seq2，3,589 顆，含作者標籤 | Q3：04–09 | ≈&nbsp;20&nbsp;MB |
| [TCGA-GBM](https://portal.gdc.cancer.gov/projects/TCGA-GBM)<br>（GDC Data Portal） | Bulk RNA-seq + 臨床存活，`10_deconv_survival.R` 經 TCGAbiolinks 下載 | Q3：10 | ≈&nbsp;1&nbsp;GB |

本倉庫不含任何資料檔與分析輸出（見 `.gitignore`）。

## 環境需求

- R ≥ 4.3（實測 4.4.1）+ RStudio，Seurat v5（實測 5.3.0）。其餘套件由 [`R/00_setup.R`](R/00_setup.R) 安裝。
- `05_infercnv.R` 需要**系統層級**的 [JAGS 4.x](https://sourceforge.net/projects/mcmc-jags/)（rjags 只是介面，R 裝不了 JAGS 本體；裝完要重開 R）。還沒裝 JAGS 前，06 之後的腳本可用作者標籤當替代惡性標籤先行測試（見 `06a` §0）。
- 選配：`08` §2–3 需要 monocle3 + SeuratWrappers（與 Bioconductor 的 monocle）；`09` §2 需要 Python 的 pySCENIC。安裝指令都在 `00_setup.R` 的選配區塊。
- Windows 使用者：`10` 的 TCGA 下載已內建處理 260 字元路徑上限（下載到 `C:/GDCdata`）；`08` §3 已內建 Monocle2 對 igraph ≥ 2.1 的相容補丁。

建議做法：新建一個 RStudio Project，把 `R/` 複製進去，從 `00_setup.R` 開始。全程不要用 `setwd()`，所有路徑相對於專案根目錄。

## 自測題庫

33 題單選（Q1 7 題、Q2 12 題、Q3 14 題），中英對照，含影片內自測題與常見錯誤、腳本練習的延伸題。互動題庫在 [`quiz/`](quiz/)：把資料夾下載到電腦後，點擊 `index.html` 即可開啟作答——即點即答、附解析、可切換中英文。

## 實作練習題（exercises/）

課程之外的動手關卡，共 23 題，全部放在 [`exercises/`](exercises/)：

- **基礎 9 題**——把課程裡的某一段流程換一份資料獨立走完，半天到一天一題。
- **進階 14 題**——多病人、多樣本的真實公開資料，從讀檔到臨床驗證走完整條線，一題約 2–3 週。**挑一題你想做的做完就好**；
  每一題都要走完「階段 D」：pseudobulk 差異表達 → cell-level DE 對照 → GSEA → CellChat →
  軌跡 → PROGENy/SCENIC 活性 → 反卷積與存活驗證，把課程教過的下游分析全部用上。

**練習題用的 23 個公開資料集，沒有一個是課程示範用的那兩份**——重跑一次示範學不到東西。
每張題卡都標明對應課程的哪一段與哪一支腳本；資料集的 accession、規模、格式與下載連結整理在
[`exercises/資料集總覽.md`](exercises/資料集總覽.md)。

## 延伸課程

本課程是快速上手篇；想繼續深入，同系列還有三個延伸課程（連結整理中，將陸續補上）：

- **A 系列 · 入門篇**——觀念與判讀：實驗設計、圖表解讀、從結果到研究敘事。全系列不需撰寫程式，適合建立整體框架與跨領域溝通。
- **B 系列 · 實作篇**——以 PBMC 公開資料完整實作：從 FASTQ 與 Cell Ranger 到 Seurat 標準流程與各項下游分析，隨集附練習腳本。本課程投影片各頁的「深入：Bx」即為對應索引。
- **C 系列 · 數學篇**——方法背後的數學與統計：PCA 與降維、變異數穩定化、統計推論（pseudobulk 的理論基礎）、整合演算法的原理。

## 目錄結構

```
scRNA-seq-tutorial-Qseries/
├── README.md
├── LICENSE                  # 程式碼：MIT
├── slides/                  # 投影片 PDF，中英各三份
├── R/                       # 十二支練習腳本（完整版）+ 練習手冊
│   └── 練習版/              # 關鍵參數挖空版（## TODO ▶ 提示）
├── exercises/               # 課後實作練習題 23 題（基礎 9、進階 14）
│   ├── 基礎/                # 換一份資料把課程流程走完
│   └── 進階/                # 發表導向的完整再分析專案
└── quiz/                    # 互動自測題庫（index.html + 題目資料）
```

## 授權

- **程式碼**（`R/`、`quiz/`）：[MIT License](LICENSE)——可自由使用、修改、再散布。
- **教材**（`slides/` 的投影片與 `exercises/` 的題卡，及其中的圖表文字）：[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/deed.zh-hant)——註明出處、非商業使用；商業授權請聯絡作者。

## 資料與工具引用

使用本課程材料發表時，請引用對應的原始資料與工具：

- **資料**：Darmanis et al. (2017) *Cell Reports*（GSE84465）；10x Genomics 公開資料集；TCGA Research Network。
- **主要工具**：Seurat v5（Hao et al. 2024）、Harmony（Korsunsky et al. 2019）、inferCNV（Broad Institute）、DoubletFinder（McGinnis et al. 2019）、scDblFinder、SoupX、SingleR、DESeq2（Love et al. 2014）、ashr、edgeR/limma 生態的 propeller/speckle、fgsea、msigdbr、clusterProfiler（Wu et al. 2021）、MAST、CellChat v2（Jin et al. 2021）、LIANA、Slingshot（Street et al. 2018）、tradeSeq、Monocle 2/3（Trapnell/Qiu/Cao et al.）、decoupleR、PROGENy、SCENIC、MuSiC（Wang et al. 2019）、TCGAbiolinks。
- **判定與狀態框架**：Tirosh et al. (2016)（CNV 相關法）、Neftel et al. (2019)（GBM 四狀態）。

---

問題與勘誤歡迎開 [Issue](../../issues)。祝分析順利——記住：每個參數是一個假設，有數字、有理由、有紀錄。
