"""Matplotlib charts shared by the Streamlit view and downloadable figures."""

from __future__ import annotations

from colorsys import hls_to_rgb, rgb_to_hls
from functools import lru_cache
from io import BytesIO
from math import ceil
from pathlib import Path
import re
from threading import RLock
from unicodedata import east_asian_width

import matplotlib as mpl
import numpy as np
from matplotlib import font_manager
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_hex, to_rgb
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator, PercentFormatter


PALETTES = {
    "湖蓝与陶土": ("#237E82", "#CD8052", "#6688AE", "#AC6F83", "#81A26B", "#AD9650", "#826D9B", "#63A5AE"),
    "色觉友好": ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#8B7600", "#4B535A"),
    "森林与浆果": ("#33765E", "#A85C7B", "#85944B", "#75678D", "#BF775A", "#477C99", "#A38D69", "#A675A5"),
}
_INK = "#213B40"
_MUTED = "#51686B"
_GRID = "#DFE7E5"
_FIGURE_LOCK = RLock()


@lru_cache(maxsize=1)
def _font_file():
    names = ("Microsoft YaHei", "PingFang SC", "Noto Sans CJK SC",
             "Noto Sans SC", "Source Han Sans CN", "Noto Sans CJK JP",
             "SimHei", "WenQuanYi Zen Hei")
    for name in names:
        try:
            return font_manager.findfont(FontProperties(family=name), fallback_to_default=False)
        except ValueError:
            continue
    for path in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
                 "C:/Windows/Fonts/msyh.ttc", "/System/Library/Fonts/PingFang.ttc"):
        if Path(path).is_file():
            return path
    return None


def _font(size, weight="normal"):
    return FontProperties(fname=_font_file(), size=size, weight=weight)


def _width(text):
    return sum(2 if east_asian_width(c) in ("W", "F") else 1 for c in str(text))


def _wrap(text, limit):
    lines = []
    for paragraph in str(text).split("\n"):
        line, used = [], 0
        # Keep ordinary English words and data IDs intact across Chinese prose.
        tokens = re.findall(r"[A-Za-z0-9_.+-]+|.", paragraph)
        for token in tokens:
            parts = list(token) if _width(token) > limit else [token]
            for piece in parts:
                size = _width(piece)
                if line and used + size > limit and piece not in "。，、；：！？）】”’.,)!%":
                    lines.append("".join(line).rstrip())
                    line, used = [], 0
                if line or not piece.isspace():
                    line.append(piece)
                    used += size
        lines.append("".join(line).rstrip())
    return "\n".join(lines)


def _colors(names, palette, is_taxonomy=False):
    base = PALETTES[palette]
    regular = 0
    colors = []
    for label in names:
        if is_taxonomy and label in ("Others", "未分类"):
            colors.append("#AEBAB9" if label == "Others" else "#D8D1BC")
            continue
        cycle, slot = divmod(regular, len(base))
        color = base[slot]
        if cycle:
            h, lightness, saturation = rgb_to_hls(*to_rgb(color))
            lightness = .33 + ((lightness - .33 + .17 * cycle) % .40)
            color = to_hex(hls_to_rgb((h + .035 * cycle) % 1, lightness, saturation))
        colors.append(color)
        regular += 1
    return colors


def _frame(width, plot_height, tick_room, legend_room, title, footer, y_label):
    left, right, top = .88, .34, .70
    wrapped = _wrap(footer, max(35, int((width - left - right) / .078)))
    footer_height = .23 + .225 * len(wrapped.splitlines())
    bottom = footer_height + legend_room + tick_room
    height = top + plot_height + bottom
    fig = Figure(figsize=(width, height), dpi=110, facecolor="white")
    FigureCanvasAgg(fig)
    ax = fig.add_axes([left / width, bottom / height, (width - left - right) / width, plot_height / height])
    fig.text(left / width, 1 - .23 / height, title, va="top", color=_INK,
             fontproperties=_font(16, "bold"))
    fig.text(left / width, (footer_height - .04) / height, wrapped, va="top",
             color=_MUTED, fontproperties=_font(9), linespacing=1.35)
    ax.set_ylabel(y_label, labelpad=8, color=_MUTED, fontproperties=_font(10))
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=_GRID, linewidth=.7)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_GRID)
    ax.tick_params(axis="both", length=0, pad=6, colors=_MUTED)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
    return fig, ax, footer_height


def _tick_style(ax):
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(_font(9))
        label.set_color(_MUTED)
    ax.yaxis.get_offset_text().set_fontproperties(_font(9))


def alpha_figure(alpha, metric, caption, palette) -> Figure:
    groups = alpha["group"].drop_duplicates().tolist()
    width = min(16, max(8.8, 4 + .74 * len(groups)))
    title = {"Shannon": "Shannon 多样性", "Simpson": "Simpson 多样性",
             "Chao1": "Chao1 丰富度估计"}[metric]
    cols = max(10, min(30, int((width - 1.3) * 12 / max(1, len(groups)))))
    data = [alpha.loc[alpha["group"] == g, metric].dropna().to_numpy(dtype=float) for g in groups]
    labels = [f"{_wrap(g, cols)}\n(n={len(values)})" for g, values in zip(groups, data)]
    # Sample-rich projects may have many groups. Keep all boxes and thin only labels.
    stride = max(1, ceil(len(groups) / 20))
    visible = list(range(0, len(groups), stride))
    tick_room = .24 + .19 * max((len(labels[i].splitlines()) for i in visible), default=1)
    footer = "每个点代表一个样本。箱体为第 25–75 百分位，中线为中位数，须线取 1.5 倍四分位距内的最远观测值。"
    if caption:
        footer += "\n" + str(caption)
    if stride > 1:
        footer += f"\n组名每隔 {stride} 组标注一次。"
    fig, ax, _ = _frame(width, 2.8, tick_room, 0, title, footer, metric)
    rng = np.random.default_rng(629)
    for i, (values, color) in enumerate(zip(data, _colors(groups, palette)), start=1):
        if not len(values):
            continue
        ax.boxplot([values], positions=[i], widths=.53, patch_artist=True, showfliers=False,
                   boxprops={"facecolor": mpl.colors.to_rgba(color, .22), "edgecolor": color, "linewidth": 1.45},
                   medianprops={"color": _INK, "linewidth": 1.65},
                   whiskerprops={"color": color, "linewidth": 1.25},
                   capprops={"color": color, "linewidth": 1.25})
        offsets = rng.permutation(np.linspace(-.15, .15, len(values))) if len(values) > 1 else [0]
        ax.scatter(i + np.asarray(offsets), values, s=26 if len(values) < 80 else 12,
                   c=color, edgecolors="white", linewidths=.55, alpha=.87, zorder=4, clip_on=False)
    ax.set_xlim(.5, max(1, len(groups)) + .5)
    ax.set_xticks([i + 1 for i in visible])
    ax.set_xticklabels([labels[i] for i in visible])
    lower, upper = ax.get_ylim()
    ax.set_ylim(max(0., lower), upper)
    _tick_style(ax)
    return fig


def _legend_layout(names, width):
    # Keep exports usable for a zero threshold on a very diverse genus table.
    names = names[:72]
    max_width = max((_width(n) for n in names), default=1)
    cell_width = max(1.30, min(3.25, .41 + .068 * max_width))
    cols = max(1, int((width - 1.22) / cell_width))
    cell_width = (width - 1.22) / cols
    max_chars = max(10, int((cell_width - .31) / .075))
    rows, height = [], .12
    for start in range(0, len(names), cols):
        row = [_wrap(name, max_chars) for name in names[start:start + cols]]
        row_height = .08 + .19 * max(len(label.splitlines()) for label in row)
        rows.append((row, row_height))
        height += row_height
    return rows, height, cell_width


def composition_figure(wide, metadata, rank, threshold_pct, palette) -> Figure:
    samples = list(wide.index)
    names = list(wide.columns)
    colors = _colors(names, palette, is_taxonomy=True)
    width = min(16., max(8.8, 4 + .38 * len(samples)))
    rows, legend_height, cell_width = _legend_layout(names, width)
    stride = max(1, ceil(len(samples) / 40))
    tick_indices = list(range(0, len(samples), stride))
    show_group = len(samples) <= 12
    labels = []
    for i in tick_indices:
        # Wrap unusually long sample IDs so they still fit inside the exported image.
        sample_label = _wrap(samples[i], 13 if show_group else 25)
        if show_group:
            sample_label += "\n" + _wrap(metadata.at[samples[i], "group"], 13)
        labels.append(sample_label)
    rotate = not show_group
    if rotate:
        tick_room = .28 + .071 * max((_width(line) for s in labels for line in s.splitlines()), default=1)
    else:
        tick_room = .24 + .19 * max((len(s.splitlines()) for s in labels), default=1)
    rank_zh = {"Phylum": "门", "Genus": "属"}[rank]
    footer = (f"每根柱代表一个样本，合计 100%。按所有展示样本的平均相对丰度，将低于 {threshold_pct:g}% 的已注释类群合并为 Others。"
              "缺少该层级注释的丰度保留在“未分类”中。")
    if stride > 1:
        footer += f" 横轴每隔 {stride} 个样本标注一次。"
    if len(names) > 72:
        footer += f" 共 {len(names)} 个类别，图例列出前 72 项，全部名称及数值见配套 CSV。"
    fig, ax, footer_height = _frame(width, 2.95, tick_room, legend_height,
                                    f"{rank_zh}水平菌群组成", footer, "相对丰度")
    baseline = np.zeros(len(samples))
    x = np.arange(len(samples), dtype=float)
    for name, color in zip(names, colors):
        values = wide[name].to_numpy(dtype=float)
        top = baseline + values
        mask = values > 0
        left, right = x[mask] - .39, x[mask] + .39
        bottom, upper = baseline[mask], top[mask]
        vertices = np.stack((np.column_stack((left, bottom)), np.column_stack((right, bottom)),
                             np.column_stack((right, upper)), np.column_stack((left, upper))), axis=1)
        bars = PolyCollection(vertices, facecolors=[color], edgecolors="none", zorder=2)
        ax.add_collection(bars)
        baseline = top
    ax.set_xlim(-.54, max(.54, len(samples) - .46))
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_xticks(tick_indices)
    ax.set_xticklabels(labels, rotation=90 if rotate else 0, ha="center")
    _tick_style(ax)
    height = fig.get_figheight()
    legend_top = footer_height + legend_height - .03
    color_index = 0
    for row, row_height in rows:
        for col, label in enumerate(row):
            xpos = .88 + cell_width * col
            fig.add_artist(Rectangle((xpos / width, (legend_top - .16) / height), .13 / width, .13 / height,
                                     facecolor=colors[color_index], edgecolor="none", transform=fig.transFigure))
            fig.text((xpos + .19) / width, legend_top / height, label, va="top", color=_MUTED,
                     fontproperties=_font(9), linespacing=1.2)
            color_index += 1
        legend_top -= row_height
    return fig


def figure_bytes(fig, fmt) -> bytes:
    """Export actual rendered data and captions, with portable SVG glyph outlines."""
    fmt = fmt.lower()
    if fmt not in {"png", "svg"}:
        raise ValueError("图表导出格式须为 PNG 或 SVG。")
    output = BytesIO()
    with _FIGURE_LOCK, mpl.rc_context({"svg.fonttype": "path"}):
        fig.savefig(output, format=fmt, dpi=300, facecolor="white")
    return output.getvalue()
