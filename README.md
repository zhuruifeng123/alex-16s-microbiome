# 菌群观察室

一个面向学生的中文 16S rRNA 肠道菌群探索工具，使用 Python 和 Streamlit 开发。打开应用即可查看模拟数据，也可以上传自己的 OTU/ASV 丰度表及分组表，比较菌群组成和样本内的 α 多样性，并下载用于汇报的图表和数据。

Alex 的示例包含 18 个 ASV、12 个模拟小鼠样本，对照组、藜麦组、燕麦组各 4 个。组名来自分组表，支持换成其他实验组。示例用于演练操作，其数值不能作为饮食干预的实验结论。

## 本地启动

先解压完整项目，进入包含 `app.py` 的目录。推荐 Python 3.12 或 3.13，本项目已在 Python 3.13 上验证。首次安装依赖需要联网。

Windows PowerShell（已安装 Python 3.13）：

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

macOS 或 Linux（已安装 Python 3.12）：

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run app.py
```

在浏览器打开终端显示的 Local URL，默认是 `http://localhost:8501`。以后启动时只需运行最后一行。在终端按 Ctrl+C 可以停止服务。Linux 上若图中文字显示为方框，请安装 Noto CJK 字体；Debian/Ubuntu 对应的系统包为 `fonts-noto-cjk`。

## 上传自己的数据

两张表都以英文逗号分隔，推荐保存为 CSV UTF-8。在侧栏可下载两张已填示例的模板。保留列名，用自己的样本名、组名和丰度替换数据。两张表的样本名区分大小写，分组表行序可以不同。

**丰度表**每行一个 OTU 或 ASV。前三个字段是 `feature_id`、`Phylum` 和 `Genus`，其余各列都是样本。每个特征 ID 唯一。未知门或属可以留空，其丰度归为“未分类”。下面的短表仅用于说明格式。

```csv
feature_id,Phylum,Genus,C01,Q01
ASV001,Bacteroidota,Bacteroides,18,27
ASV002,Bacillota,,2,1
```

**分组表**每个样本只出现一次，并填写组名。

```csv
sample_id,group
Q01,藜麦组
C01,对照组
```

在侧栏选择“上传 CSV”，放入两张表并声明输入类型。原始计数须为非负整数。相对丰度可使用数字比例（如 0.25）或数字百分数（如 25），同一样本内使用相同单位。空丰度不会自动填零；确认是未检出时再填 0。应用检查重复 ID、缺表头、非数值、空值、负值、总丰度为零的样本和样本名错配，并给出更正位置。任一文件更新后会重新校验。

## 如何阅读和导出结果

| 视图 | 计算口径 | 可以从图中观察什么 |
| --- | --- | --- |
| 门或属组成图 | 在各样本中先换算相对丰度，再按当前分类层级汇总。已注释类群在所有展示样本中的平均占比严格低于阈值时合并为 Others。“未分类”始终保留。 | 比较同一类群在各样本中的占比。每根柱子合计 100%。阈值设为 0 可关闭合并。相对丰度本身不能确定绝对数量的变化。 |
| Shannon | H = −Σpᵢ ln pᵢ，零丰度项贡献为 0。 | 同时观察样本内的丰富度和均匀度。 |
| Simpson | D = 1 − Σpᵢ²。 | 衡量随机独立抽取两次得到不同特征的概率。该值的大小方向与 Σpᵢ² 版本相反。 |
| Chao1 | 偏差修正估计 Sobs + F₁(F₁−1) / [2(F₂+1)]。Sobs 为有观测的特征数，F₁、F₂ 为恰好出现一次、两次的特征数。 | 从稀有特征估计样本丰富度。F₂ 为 0 时分母为 2；F₁ 为 0 或 1 时结果等于 Sobs。只有原始计数能提供这些观测频数。 |

α 指数在原始 OTU/ASV 层级逐样本计算，切换门、属或调整合并阈值不会改变指数。箱体为第 25–75 百分位，中线为中位数，散点代表所有样本。输入相对丰度时，页面解释 Chao1 的限制，并保留 Shannon、Simpson 和组成图。

在“α 多样性”中打开“标注探索性组间比较”，选定指标、两组和双侧检验。Welch t 检验比较均值，允许方差不同。Mann–Whitney U 检验比较秩分布；解释为中位数差异还需要分布形状等条件。图注给出组名、样本量、方法和原始 p 值，并标明未做多重比较校正。本工具要求每组至少 2 个有效样本；样本不足、两组方差均为零而无法计算 Welch 检验等情形，会显示原因。重复测量、配对或同笼关联需要单独的统计设计。

每张图下方都有 **PNG（300 dpi）**、**SVG（矢量）** 和 **绘图数据 CSV** 下载按钮。组成图导出当前层级和阈值下的百分数长表，字段为 `sample_id,group,rank,taxon,relative_abundance_pct`。多样性图导出每个样本的组名和当前指标，用于逐点复核。统计图注随 PNG 和 SVG 保存。CSV 使用带 BOM 的 UTF-8 编码，便于在 Excel 中查看中文。

## 三个首版验收操作

1. **载入示例并导出。** 点击“载入示例”，确认显示 12 个样本、3 个分组和 18 个 OTU/ASV。进入“菌群组成”，切换至属水平，阈值保留为 1%。图例应有 Others 和“未分类”。下载 PNG、SVG、CSV；在 CSV 中按 `sample_id` 汇总百分数，每个样本应约为 100%。再到“α 多样性”选择 Shannon，开启 Welch 检验，选“对照组”和“藜麦组”。两组 n 均为 4，示例原始 p 值约为 0.159。
2. **故意改错一个样本名。** 复制 `data/example_groups.csv`，将其中的 `C01` 改为 `C01_typo`。选择“上传 CSV”，上传原丰度表与改动的分组表。页面应指出仅存在于丰度表的 `C01` 和仅存在于分组表的 `C01_typo`。改回正确名称并替换文件，图表应恢复。
3. **上传真正的相对丰度。** 在项目目录内，用上面安装了依赖的 Python 解释器运行下面代码，将生成 `relative_example.csv`。在“上传 CSV”中上传此表和原分组表，选择“相对丰度”。组成图仍可查看；“α 多样性”仅提供 Shannon 和 Simpson，页面会解释 Chao1 需要原始计数。按相同样本核对，前两个指标应与原计数版本一致（浮点舍入误差除外）。

```python
import pandas as pd

table = pd.read_csv("data/example_abundance.csv")
samples = table.drop(columns=["feature_id", "Phylum", "Genus"])
relative = samples.div(samples.sum(axis=0), axis="columns")
pd.concat([table[["feature_id", "Phylum", "Genus"]], relative], axis=1).to_csv(
    "relative_example.csv", index=False, encoding="utf-8-sig"
)
```

可重复运行核心计算回归测试。在 Windows 项目目录内执行：

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

## GitHub 仓库组织

把以下内容放在同一个仓库根目录，包括以点开头的 `.streamlit` 目录。将来上传自己的研究数据时可通过页面上传；仓库内的 `data/` 用于存放公开演示样例。

```text
Alex16S/
├── app.py                     # Streamlit 入口、交互和内置中文说明
├── analysis_core.py           # 文件校验、指数、组间检验、组成汇总
├── charts.py                  # 绘图、中文字体和 PNG/SVG 导出
├── requirements.txt           # Python 依赖
├── packages.txt               # 云端 Linux 中文字体
├── .streamlit/config.toml      # 主题
├── .gitignore
├── data/
│   ├── example_abundance.csv
│   └── example_groups.csv
├── tests/test_analysis.py
└── README.md
```

## 部署到 Streamlit Community Cloud

1. 新建 GitHub 仓库，将上述目录中的文件提交到默认分支，确认在 GitHub 打开仓库时能直接看到 `app.py`、`requirements.txt` 和 `packages.txt`。
2. 登录 Streamlit Community Cloud，连接可访问该仓库的 GitHub 账号，选择创建应用。填写仓库、分支和入口 `app.py`。如果你把整个文件夹放在已有仓库的子目录，请相应填写入口路径，并按 Cloud 文档放置依赖文件。[Community Cloud 部署说明](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app)
3. 在高级设置中选择 Python 3.12 或 3.13，开始部署。平台会安装 Python 依赖和 `packages.txt` 中的中文字体。应用使用模拟数据即可启动，无需配置密钥。
4. 部署成功后打开平台分配的应用网址，用上面的三个操作复查线上版本。需要对外展示时，在应用分享设置中选择合适的访问权限。

## 1–2 分钟演示提纲

| 时间 | 屏幕操作与口述重点 |
| --- | --- |
| 0:00–0:15 | 点击“载入示例”。介绍研究问题和两张输入表，指出屏幕上是三组共 12 个模拟样本。 |
| 0:15–0:40 | 打开属水平组成图，指出一个优势属、“未分类”和 Others；调整阈值，解释为什么柱高仍为 100%。 |
| 0:40–1:10 | 切到 Shannon 箱线图，指着一个点说明它代表一只小鼠。比较对照组与藜麦组，读出 n 和原始 p，说明是探索性比较。 |
| 1:10–1:35 | 下载图片和当前绘图数据。展示本地准备好的相对丰度表，指出切换输入类型后 Chao1 的限制说明。 |
| 1:35–1:50 | 打开“使用说明”，说明完整研究仍需原始数据质控、合适的统计设计和后续验证。 |

## 分析边界及后续扩展

本工具适合已整理数据的学生探索与汇报。测序深度影响特征检出和多样性估计，首版直接使用上传值，不进行稀释抽平、原始序列质控或差异丰度建模。正式研究仍需结合实验设计，选用 QIIME 2、MicrobiomeAnalyst 等工具及恰当的统计方法。

可在后续版本增加基于 Bray–Curtis 距离的 PCoA、按组别筛选样本、饮食或抗生素干预实验模板，以及一键生成包含文字和图表的 Markdown/PDF 简报。
