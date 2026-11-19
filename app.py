"""Streamlit entry point for the student 16S microbiome explorer."""

from hashlib import sha256
from pathlib import Path
from threading import RLock

import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt

from analysis_core import alpha_table, compare_groups, composition_table, read_inputs
from charts import PALETTES, alpha_figure, composition_figure, figure_bytes


DATA = Path(__file__).resolve().parent / "data"
VIEWS = ("菌群组成", "α 多样性", "使用说明")
KINDS = {"原始计数": "counts", "相对丰度": "relative"}
METHODS = {"Welch t 检验": "welch", "Mann–Whitney U 秩和检验": "mannwhitney"}
METRIC_HELP = {
    "Shannon": "同时反映丰富度和均匀度。样本内各类群的丰度越均匀，其他条件相同时，指数越高。",
    "Simpson": "这里采用 1 − Σpᵢ²，表示从同一样本随机独立抽取两次，得到不同 OTU/ASV 的概率。",
    "Chao1": "利用只观测到一次和两次的 OTU/ASV 数量，估计样本内的总丰富度。",
}
RELATIVE_NOTE = "当前输入为相对丰度。Chao1 需要原始计数来识别恰好出现 1 次和 2 次的 OTU/ASV。此处仅提供 Shannon、Simpson 和菌群组成。"


def csv_bytes(frame):
    return frame.to_csv(index=False, lineterminator="\n").encode("utf-8-sig")


def load_demo():
    st.session_state["source"] = "示例数据"


def input_panel():
    demo_abundance = (DATA / "example_abundance.csv").read_bytes()
    demo_groups = (DATA / "example_groups.csv").read_bytes()
    with st.sidebar:
        st.header("数据输入")
        source = st.radio("数据来源", ("示例数据", "上传 CSV"), key="source", horizontal=True)
        if source == "示例数据":
            st.caption("模拟小鼠数据 · 对照组、藜麦组、燕麦组")
            inputs = (demo_abundance, demo_groups, "counts")
        else:
            abundance = st.file_uploader("丰度表 CSV", type="csv", key="abundance_file")
            groups = st.file_uploader("分组表 CSV", type="csv", key="groups_file")
            declared = st.radio("丰度表中的数值", tuple(KINDS), key="input_kind",
                                help="请选择测序原始计数，或已经换算的比例或百分数。")
            inputs = (abundance.getvalue() if abundance is not None else None,
                      groups.getvalue() if groups is not None else None, KINDS[declared])
        st.button("载入示例", on_click=load_demo, width="stretch")
        st.divider()
        st.subheader("CSV 模板")
        st.caption("两张模板使用相同的样本名，可替换为自己的数据。")
        st.download_button("丰度表模板", demo_abundance, "abundance_template.csv", "text/csv",
                           width="stretch", on_click="ignore")
        st.download_button("分组表模板", demo_groups, "groups_template.csv", "text/csv",
                           width="stretch", on_click="ignore")
    return source, inputs


@st.cache_resource
def plot_lock():
    # Matplotlib's renderer and font manager are shared across Streamlit sessions.
    return RLock()


def render_chart(signature, build, plotted_data, file_stem, caption):
    # Only retain the most recent image pair, in the current student's session.
    saved = st.session_state.get("_last_chart")
    if saved is None or saved[0] != signature:
        with plot_lock():
            fig = build()
            try:
                png, svg = figure_bytes(fig, "png"), figure_bytes(fig, "svg")
            finally:
                plt.close(fig)
        st.session_state["_last_chart"] = (signature, png, svg)
    else:
        _, png, svg = saved
    st.image(png, width="stretch")
    st.caption(caption)
    buttons = st.columns(3)
    for column, (label, content, ext, mime) in zip(buttons, [
        ("下载 PNG · 300 dpi", png, "png", "image/png"),
        ("下载 SVG · 矢量图", svg, "svg", "image/svg+xml"),
        ("下载绘图数据 · CSV", csv_bytes(plotted_data), "csv", "text/csv"),
    ]):
        column.download_button(label, content, f"{file_stem}.{ext}", mime, width="stretch",
                               on_click="ignore", key=f"{file_stem}_{ext}")
    with st.expander("查看当前图的数值"):
        st.dataframe(plotted_data, hide_index=True, width="stretch")


def show_composition(dataset, token):
    st.subheader("菌群组成")
    st.write("每根柱子代表一个样本。沿同一种颜色横向比较，可以看到该类群在不同样本中所占的比例。")
    rank_col, threshold_col, palette_col = st.columns([1, 1.4, 1])
    rank_label = rank_col.selectbox("分类层级", ["门（Phylum）", "属（Genus）"])
    rank = "Phylum" if rank_label.startswith("门") else "Genus"
    threshold = threshold_col.slider("低丰度合并阈值（%）", 0.0, 100.0, 1.0, 0.1,
                                     help="所有样本的平均相对丰度低于此值时合并。设为 0 可关闭合并。")
    palette = palette_col.selectbox("图表配色", list(PALETTES), key="composition_palette")
    wide = composition_table(dataset, rank, threshold)
    # Construct the long export from the same matrix that draws the bars.
    long = wide.rename_axis(index="sample_id", columns="taxon").stack().rename(
        "relative_abundance_pct").reset_index()
    long.insert(1, "group", long["sample_id"].map(dataset.metadata["group"]))
    long.insert(2, "rank", rank)
    rule = (f"每个样本归一到 100%。已注释类群按所有展示样本的平均相对丰度判断，"
            f"低于 {threshold:g}% 的合入 Others。缺失注释单列为“未分类”。")
    render_chart((token, rank, threshold, palette),
                 lambda: composition_figure(wide, dataset.metadata, rank, threshold, palette),
                 long, f"composition_{rank}", rule)


def show_alpha(dataset, kind, token):
    st.subheader("α 多样性")
    st.write("每个点代表一个样本。先看同一组内各点的分散程度，再比较不同组的中位数。")
    if kind == "relative":
        st.info(RELATIVE_NOTE)
    table = alpha_table(dataset, kind)
    available = list(METRIC_HELP) if kind == "counts" else ["Shannon", "Simpson"]
    left, right = st.columns(2)
    metric = left.selectbox("多样性指标", available)
    palette = right.selectbox("图表配色", list(PALETTES), key="alpha_palette")
    st.caption(METRIC_HELP[metric])
    show_test = st.toggle("标注探索性组间比较", value=False,
                          help="选择两个独立组进行一次双侧检验。箱线图仍展示全部组。")
    caption = ""
    if show_test:
        groups = table["group"].drop_duplicates().tolist()
        if len(groups) < 2:
            caption = "仅有一个分组。至少需要两个独立组才能比较。"
        else:
            a, b, m = st.columns([1, 1, 1.5])
            group_a = a.selectbox("第一组", groups)
            group_b = b.selectbox("第二组", [group for group in groups if group != group_a])
            method = m.selectbox("双侧检验方法", list(METHODS))
            _, caption = compare_groups(table, metric, group_a, group_b, METHODS[method])
    description = caption or "箱体表示第 25–75 百分位数，中线为中位数。散点显示全部样本。"
    render_chart((token, metric, palette, caption),
                 lambda: alpha_figure(table, metric, caption, palette),
                 table[["sample_id", "group", metric]], f"alpha_{metric}", description)


def show_guide():
    st.subheader("如何准备两张 CSV")
    st.markdown("""
1. **丰度表**第一行写 `feature_id,Phylum,Genus`，随后每列是一个样本。每行一个 OTU/ASV，样本丰度填非负数。未检出填 0，未知分类可以留空。
2. **分组表**包含 `sample_id,group`。样本名与丰度表逐一对应，顺序可以不同。两个文件均以逗号分隔，推荐从 Excel 另存为 CSV UTF-8。
3. 在侧栏选择“上传 CSV”，放入两张表并声明数值类型。相对丰度可填比例或百分数，使用数值即可，例如 0.125 或 12.5。同一样本使用同一单位。
""")
    st.subheader("读懂三个 α 多样性指标")
    st.write("计算单位是单个样本中的原始 OTU/ASV。先用各特征丰度除以样本总丰度得到 pᵢ。门、属切换和低丰度合并只影响组成图。")
    st.markdown("**Shannon**　H = −Σ pᵢ ln pᵢ，只对 pᵢ > 0 的项求和。它同时反映丰富度和均匀度。")
    st.markdown("**Simpson**　D = 1 − Σ pᵢ²。单一类群时为 0；优势类群集中时通常较低。")
    st.markdown("**Chao1**　采用偏差修正公式。Sobs 为有计数的特征数，F₁ 和 F₂ 分别为恰好出现 1 次和 2 次的特征数。")
    st.latex(r"\widehat S_{\rm Chao1}=S_{\rm obs}+\frac{F_1(F_1-1)}{2(F_2+1)}")
    st.write("当 F₂ = 0 时，分母为 2；当 F₁ 为 0 或 1 时，估计值等于 Sobs。总计数为零的样本需在输入阶段处理。相对丰度无法还原单次或双次观测，因此无法估计 Chao1。")
    st.subheader("组成图和组间检验各回答什么")
    st.markdown("""
**组成图**展示每个样本在门或属水平的丰度比例，不能直接说明细菌绝对数量的变化。Others 汇总低于当前阈值的已注释类群；“未分类”保留缺失注释的丰度。

**Welch t 检验**比较两个独立组的均值，允许方差不同。**Mann–Whitney U 检验**比较秩分布；将结果解释为位置或中位数差异还需要分布形状等条件。两者均为双侧探索性检验。检验之前需确认样本独立；重复测量、同笼效应等需要与研究设计相应的处理。这里每组至少需有 2 个有效样本，样本少时请结合每个散点解读。

**原始 p 值**是在检验假设下衡量观测差异的一个指标。反复尝试指标或组别会增加偶然发现差异的机会。页面展示原始 p 值，未做多重比较校正。
""")
    st.subheader("适用范围")
    st.write("本工具从已整理的丰度表开始，适合学生探索及汇报。示例全部为模拟数据。测序深度会影响特征检出和丰富度；这里没有做稀释抽平。正式研究还需根据设计完成质控、混杂因素处理及多重检验等，可配合 QIIME 2、MicrobiomeAnalyst 等工具。")


def main():
    st.set_page_config(page_title="菌群观察室 · 16S rRNA", page_icon="🧫", layout="wide",
                       initial_sidebar_state="auto")
    st.markdown("""<style>
    .stApp {font-family: "Inter","Microsoft YaHei",sans-serif;}
    [data-testid="stMainBlockContainer"] {max-width: 1140px; padding-top: 2.4rem;}
    [data-testid="stSidebar"] {border-right: 1px solid #dce6e1;}
    [data-testid="stSidebarContent"] h2 {font-size: 1.42rem;}
    .st-key-page_nav {border-bottom: 1px solid #d8e3de; padding-bottom: .35rem;}
    .st-key-page_nav [role="radiogroup"] {gap: .4rem;}
    .st-key-page_nav [role="radiogroup"] label {padding: .28rem .6rem; border-radius: .45rem;}
    .st-key-page_nav label:has(input:checked) {background: #e4f1eb;}
    @media (max-width: 700px) {[data-testid="stMainBlockContainer"] {padding-top: 1.2rem;}}
    </style>""", unsafe_allow_html=True)
    source, (abundance, metadata, kind) = input_panel()
    st.title("菌群观察室")
    st.write("从 16S rRNA 丰度表查看样本的菌群组成与多样性。")
    with st.container(key="page_nav"):
        view = st.radio("分析视图", VIEWS, horizontal=True, label_visibility="collapsed")
    if view == "使用说明":
        show_guide()
        return
    if abundance is None or metadata is None:
        st.info("请在数据输入区上传丰度表和分组表，或点击“载入示例”查看模拟数据。")
        return
    try:
        dataset = read_inputs(abundance, metadata, kind)
    except ValueError as error:
        st.error(str(error))
        return
    sample_count = dataset.counts.shape[1]
    group_count = dataset.metadata["group"].nunique()
    label = "模拟示例" if source == "示例数据" else "上传数据"
    kind_label = "原始计数" if kind == "counts" else "相对丰度"
    st.caption(f"{label}　·　{sample_count} 个样本　·　{group_count} 个分组　·　"
               f"{len(dataset.counts)} 个 OTU/ASV　·　{kind_label}")
    if dataset.warnings:
        with st.expander(f"数据检查提示（{len(dataset.warnings)}）"):
            for warning in dataset.warnings:
                st.write(warning)
    if kind == "relative" and view != "α 多样性":
        st.info(RELATIVE_NOTE)
    token = (sha256(abundance).digest(), sha256(metadata).digest(), kind)
    if view == "菌群组成":
        show_composition(dataset, token)
    else:
        show_alpha(dataset, kind, token)


if __name__ == "__main__":
    main()
