"""读取 16S 丰度表，计算样本多样性、分类组成和探索性组间比较。"""

from __future__ import annotations

import csv
import io
import warnings
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class Dataset:
    counts: pd.DataFrame  # 原始输入丰度，行是 feature_id，列是样本。
    taxonomy: pd.DataFrame
    metadata: pd.DataFrame
    warnings: tuple[str, ...]
    input_kind: str


def _names(items) -> str:
    items = list(items)
    preview = "、".join(f"“{value}”" for value in items[:6])
    return preview + (f" 等 {len(items)} 项" if len(items) > 6 else "")


def _read_csv(data: bytes, label: str, required: tuple[str, ...]) -> pd.DataFrame:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        try:
            text = data.decode("gb18030")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{label}编码无法识别，请另存为 CSV UTF-8 后上传。") from exc

    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    try:
        header = [name.strip() for name in next(reader, [])]
        if not header:
            raise ValueError(f"{label}没有表头，请下载并填写模板。")
        if "" in header:
            raise ValueError(f"{label}第 {header.index('') + 1} 列缺少表头，请填入列名或移除空列。")
        duplicates = [key for key, count in Counter(header).items() if count > 1]
        if duplicates:
            raise ValueError(f"{label}表头重复：{_names(duplicates)}。请为每列填写唯一名称。")
        missing = [name for name in required if name not in header]
        if missing:
            raise ValueError(f"{label}缺少必需列 {_names(missing)}。请按模板拼写列名，并用逗号分隔。")

        rows = []
        for row in reader:
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"{label}第 {reader.line_num} 行有 {len(row)} 列，表头有 {len(header)} 列。"
                    "请补齐字段，含逗号的名称应放在英文双引号内。"
                )
            rows.append([value.strip() for value in row])
    except csv.Error as exc:
        raise ValueError(f"{label}第 {reader.line_num} 行 CSV 格式错误，请检查引号或按模板重新导出。") from exc
    if not rows:
        raise ValueError(f"{label}只有表头，请至少填写一条记录。")
    return pd.DataFrame(rows, columns=header)


def _index_by_id(frame: pd.DataFrame, key: str, label: str) -> pd.DataFrame:
    ids = frame[key]
    if ids.eq("").any():
        position = int(np.flatnonzero(ids.eq("").to_numpy())[0]) + 1
        raise ValueError(f"{label}第 {position} 条数据的 {key} 为空，请补充唯一 ID。")
    repeated = ids[ids.duplicated()].unique()
    if len(repeated):
        raise ValueError(f"{label}的 {key} 重复：{_names(repeated)}。请核对记录，让每个 ID 仅出现一次。")
    return frame.set_index(key)


def _numeric_sample(values: pd.Series, sample: str, input_kind: str) -> pd.Series:
    def reject(mask, message: str) -> None:
        if np.any(mask):
            feature = values.index[int(np.flatnonzero(mask)[0])]
            raise ValueError(f"丰度表样本“{sample}”、特征“{feature}”{message}")

    reject(values.eq("").to_numpy(), "的丰度为空。请补充数值，确认未检出时可填写 0。")
    if input_kind == "relative":
        try:
            # 保留以长小数字符串表示的罕见类群，避免短精度十进制解析提前舍入至零。
            numeric = values.astype(float)
        except (ValueError, OverflowError):
            numeric = pd.to_numeric(values, errors="coerce")
    else:
        numeric = pd.to_numeric(values, errors="coerce")
    array = numeric.to_numpy(dtype=float)
    reject(~np.isfinite(array), "的丰度不是有限数值。请填写数字，去掉百分号、文字或无穷大。")
    reject(array < 0, "的丰度为负值。请核对原表，填写非负丰度。")
    if input_kind == "counts":
        reject(array != np.floor(array), "含小数。原始计数须为非负整数；比例或百分数请选“相对丰度”。")
        # Decimal 只用于带小数点/指数的计数，避免浮点舍入把接近整数的比例判作整数。
        unusual = values[values.str.contains(r"[.eE]", regex=True)]
        for feature, raw in unusual.items():
            exact = Decimal(raw)
            if exact != exact.to_integral_value() or exact < 0:
                raise ValueError(
                    f"丰度表样本“{sample}”、特征“{feature}”含非整数计数。"
                    "请核对原表；比例或百分数请选“相对丰度”。"
                )
    if not np.any(array > 0):
        raise ValueError(f"样本“{sample}”总丰度为零。请核对计数，或在两张表中同时移除这个空样本。")
    return numeric


def read_inputs(abundance: bytes, metadata: bytes, input_kind: str) -> Dataset:
    """校验丰度与分组并对齐样本。相对丰度支持比例或百分数。"""
    if input_kind not in ("counts", "relative"):
        raise ValueError("请选择输入类型：原始计数或相对丰度。")
    abundance_keys = ("feature_id", "Phylum", "Genus")
    raw = _index_by_id(_read_csv(abundance, "丰度表", abundance_keys), "feature_id", "丰度表")
    meta = _index_by_id(_read_csv(metadata, "分组表", ("sample_id", "group")), "sample_id", "分组表")
    samples = [name for name in raw.columns if name not in ("Phylum", "Genus")]
    if not samples:
        raise ValueError("丰度表缺少样本列。请在分类列后添加样本 ID，并在分组表填写相同的 sample_id。")
    missing_groups = [name for name in samples if name not in meta.index]
    sample_set = set(samples)
    missing_counts = [name for name in meta.index if name not in sample_set]
    if missing_groups or missing_counts:
        details = []
        if missing_groups:
            details.append(f"仅在丰度表中出现 {_names(missing_groups)}")
        if missing_counts:
            details.append(f"仅在分组表中出现 {_names(missing_counts)}")
        raise ValueError("样本名不匹配：" + "；".join(details) + "。请统一两个文件中的样本 ID（含大小写），再重新上传。")
    meta = meta.loc[samples, ["group"]]
    empty_groups = meta.index[meta["group"].eq("")]
    if len(empty_groups):
        raise ValueError(f"分组表样本 {_names(empty_groups)} 的 group 为空，请补充组名。")
    counts = pd.DataFrame({sample: _numeric_sample(raw[sample], sample, input_kind) for sample in samples})

    notes = []
    taxonomy = raw.loc[:, ["Phylum", "Genus"]].copy()
    missing_labels = {"", "na", "n/a", "nan", "null", "none", "unknown", "unassigned", "unclassified", "未分类"}
    for rank in taxonomy:
        missing = taxonomy[rank].str.lower().isin(missing_labels)
        if missing.any():
            taxonomy.loc[missing, rank] = "未分类"
            notes.append(f"{rank} 有 {int(missing.sum())} 个特征缺少分类注释，其丰度保留在“未分类”中。")
    if input_kind == "relative":
        with np.errstate(over="ignore"):
            totals = counts.to_numpy(dtype=float).sum(axis=0)
        unusual = ~(np.isclose(totals, 1, rtol=0.02) | np.isclose(totals, 100, rtol=0.02))
        if np.any(unusual):
            names = [sample for sample, flag in zip(samples, unusual) if flag]
            notes.append(f"样本 {_names(names)} 的丰度合计偏离 1 或 100 超过 2%。已按各样本自己的合计归一化，请核对输入单位与表格完整性。")
    return Dataset(counts, taxonomy, meta, tuple(notes), input_kind)


def _proportions(dataset: Dataset) -> pd.DataFrame:
    # 先除以各列最大值，防止直接求和出现浮点溢出。
    scaled = dataset.counts.astype(float).div(dataset.counts.max(axis=0), axis="columns")
    return scaled.div(scaled.sum(axis=0), axis="columns")


def alpha_table(dataset: Dataset, input_kind: str) -> pd.DataFrame:
    """计算原始特征层级的 Shannon、Gini–Simpson 和偏差修正 Chao1。"""
    if input_kind != dataset.input_kind:
        raise ValueError("输入类型已改变，请按当前类型重新校验两张 CSV。")
    props = _proportions(dataset).to_numpy()
    logs = np.zeros_like(props)
    np.log(props, out=logs, where=props > 0)
    table = dataset.metadata.copy()
    table["Shannon"] = np.maximum(-(props * logs).sum(axis=0), 0)
    table["Simpson"] = np.clip(1 - (props * props).sum(axis=0), 0, 1)
    if input_kind == "counts":
        counts = dataset.counts
        observed = counts.gt(0).sum(axis=0)
        singletons = counts.eq(1).sum(axis=0).astype(float)
        doubletons = counts.eq(2).sum(axis=0)
        table["Chao1"] = observed + singletons * (singletons - 1) / (2 * (doubletons + 1))
    else:
        table["Chao1"] = np.nan
    return table.reset_index()


def composition_table(dataset: Dataset, rank: str, threshold_pct: float) -> pd.DataFrame:
    """返回样本 × 分类单元的百分数表，列名与图例一致。"""
    if rank not in ("Phylum", "Genus"):
        raise ValueError("分类层级应选 Phylum（门）或 Genus（属）。")
    if not 0 <= threshold_pct <= 100:
        raise ValueError("低丰度阈值须在 0 至 100% 之间。")
    taxa = dataset.taxonomy[rank].copy()
    if taxa.eq("Others").any():
        # 真实注释中也可能有名为 Others 的标签，避免覆盖自动合并的类群。
        original = "Others（原始注释）"
        existing = set(taxa)
        while original in existing:
            original += "·"
        taxa = taxa.replace("Others", original)
    grouped = _proportions(dataset).groupby(taxa, sort=True).sum().T * 100
    means = grouped.mean(axis=0).sort_values(ascending=False, kind="stable")
    # 为恰好落在阈值上的分类保留极小的浮点误差余量。
    small = (means < threshold_pct) & ~np.isclose(means, threshold_pct, rtol=1e-12, atol=0)
    low = [name for name in means.index[small] if name != "未分类"]
    low_names = set(low)
    kept = [name for name in means.index if name not in low_names and name != "未分类"]
    if "未分类" in grouped.columns:
        kept.append("未分类")
    result = grouped[kept].copy()
    if low:
        result["Others"] = grouped[low].sum(axis=1)
    result.index.name = "sample_id"
    result.columns.name = None
    return result


def compare_groups(
    alpha: pd.DataFrame, metric: str, group_a: str, group_b: str, method: str
) -> tuple[float | None, str]:
    """双侧独立样本比较；无有效结果时返回 None 和可直接用于图注的原因。"""
    methods = {"welch": "Welch t 检验", "mannwhitney": "Mann–Whitney U 秩和检验"}
    if method not in methods:
        raise ValueError("请选择 Welch t 检验或 Mann–Whitney U 秩和检验。")
    if metric not in ("Shannon", "Simpson", "Chao1"):
        raise ValueError("请选择 Shannon、Simpson 或 Chao1 指标。")
    a = alpha.loc[alpha["group"].eq(group_a), metric].dropna().to_numpy(dtype=float)
    b = alpha.loc[alpha["group"].eq(group_b), metric].dropna().to_numpy(dtype=float)
    label = f"{metric} | {group_a}（n={len(a)}）与{group_b}（n={len(b)}） | 双侧 {methods[method]}"

    def unavailable(reason: str) -> tuple[None, str]:
        return None, f"{label} | {reason}，无可报告的 p 值。未做多重比较校正。"

    if group_a == group_b:
        return unavailable("请选择两个不同的独立组")
    if min(len(a), len(b)) < 2:
        return unavailable("样本量不足，每组至少需要 2 个有效样本")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return unavailable("样本指标存在非有限数值")
    if method == "welch" and np.ptp(a) == 0 and np.ptp(b) == 0:
        return unavailable("两组组内方差均为零，无法估计标准误")
    if method == "mannwhitney" and np.all(a == a[0]) and np.all(b == a[0]):
        return unavailable("所有观测值均相同，秩次无法区分")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            if method == "welch":
                p = float(stats.ttest_ind(a, b, equal_var=False, alternative="two-sided").pvalue)
            else:
                p = float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue)
    except (ArithmeticError, ValueError) as exc:
        return unavailable(f"检验无法计算（{exc}）")
    if not np.isfinite(p):
        return unavailable("统计量无有效数值，请检查样本的变异程度")
    display = f"{p:.5g}" if p > 0 else "0（低于数值精度）"
    return p, f"{label} | 原始 p = {display}，未做多重比较校正。"
