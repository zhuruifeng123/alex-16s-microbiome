"""在项目根目录运行 python -m unittest discover -s tests。"""

import io
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analysis_core import alpha_table, compare_groups, composition_table, read_inputs


DATA = Path(__file__).resolve().parents[1] / "data"
SMALL_COUNTS = (
    "feature_id,Phylum,Genus,S1,S2\n"
    "f1,门甲,属甲,2,0\nf2,门甲,属甲,1,4\nf3,门乙,,1,0\n"
).encode("utf-8")
SMALL_GROUPS = "sample_id,group\nS2,乙\nS1,甲\n".encode("utf-8")


class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.abund = (DATA / "example_abundance.csv").read_bytes()
        cls.groups = (DATA / "example_groups.csv").read_bytes()

    def test_demo_joins_by_identifier_not_position(self):
        dataset = read_inputs(self.abund, self.groups, "counts")
        alpha = alpha_table(dataset, "counts").set_index("sample_id")
        self.assertEqual(len(alpha), 12)
        self.assertEqual(alpha.at["C01", "group"], "对照组")
        self.assertEqual(alpha.at["Q02", "group"], "藜麦组")
        self.assertEqual(alpha.at["O04", "group"], "燕麦组")
        self.assertTrue(np.isfinite(alpha[["Shannon", "Simpson", "Chao1"]]).all().all())
        for rank in ("Phylum", "Genus"):
            with self.subTest(rank=rank):
                wide = composition_table(dataset, rank, 1.0)
                self.assertEqual(set(wide.index), set(alpha.index))
                np.testing.assert_allclose(wide.sum(axis=1), 100)
        genera = composition_table(dataset, "Genus", 1.0)
        self.assertIn("Others", genera.columns)
        self.assertIn("未分类", genera.columns)

    def test_wrong_sample_id_identifies_mismatch(self):
        wrong = self.groups.replace(b"Q02,", b"Q02_typo,", 1)
        with self.assertRaises(ValueError) as caught:
            read_inputs(self.abund, wrong, "counts")
        self.assertIn("Q02_typo", str(caught.exception))

    def test_correct_alpha_uses_features_before_merging_genera(self):
        dataset = read_inputs(SMALL_COUNTS, SMALL_GROUPS, "counts")
        alpha = alpha_table(dataset, "counts").set_index("sample_id")
        self.assertAlmostEqual(alpha.at["S1", "Shannon"], 1.5 * math.log(2))
        self.assertAlmostEqual(alpha.at["S1", "Simpson"], 0.625)
        self.assertAlmostEqual(alpha.at["S1", "Chao1"], 3.5)
        self.assertEqual(alpha.at["S2", "Chao1"], 1)
        wide = composition_table(dataset, "Genus", 90)
        self.assertAlmostEqual(wide.at["S1", "未分类"], 25)
        np.testing.assert_allclose(wide.sum(axis=1), 100)

    def test_relative_keeps_diversity_and_composition_and_omits_chao1(self):
        table = pd.read_csv(io.BytesIO(self.abund))
        samples = table.columns[3:]
        relative = pd.concat([table.iloc[:, :3], table[samples].div(table[samples].sum())], axis=1)
        dataset = read_inputs(relative.to_csv(index=False).encode("utf-8"), self.groups, "relative")
        alpha = alpha_table(dataset, "relative").set_index("sample_id")
        baseline = read_inputs(self.abund, self.groups, "counts")
        raw = alpha_table(baseline, "counts").set_index("sample_id").loc[alpha.index]
        self.assertTrue(alpha["Chao1"].isna().all())
        np.testing.assert_allclose(alpha[["Shannon", "Simpson"]], raw[["Shannon", "Simpson"]])
        np.testing.assert_allclose(
            composition_table(dataset, "Genus", 1.0),
            composition_table(baseline, "Genus", 1.0), atol=1e-10,
        )

    def test_fractional_count_negative_empty_and_text_identify_column(self):
        for invalid in (b"1.5", b"-1", b"", b"typo"):
            with self.subTest(value=invalid):
                table = SMALL_COUNTS.replace(b",2,0\n", b"," + invalid + b",0\n", 1)
                with self.assertRaisesRegex(ValueError, "S1"):
                    read_inputs(table, SMALL_GROUPS, "counts")

    def test_zero_sum_duplicate_ids_and_missing_headers(self):
        cases = [
            (SMALL_COUNTS.replace(b",1,4\n", b",1,0\n"), SMALL_GROUPS, "S2"),
            (SMALL_COUNTS.replace(b"f2,", b"f1,"), SMALL_GROUPS, "f1"),
            (SMALL_COUNTS.replace(b",S1,S2", b",S1,S1"), SMALL_GROUPS, "S1"),
            (SMALL_COUNTS, SMALL_GROUPS.replace(b"S2,", b"S1,"), "S1"),
            (SMALL_COUNTS.replace(b"Phylum", b""), SMALL_GROUPS, ""),
        ]
        for table, groups, expected in cases:
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(ValueError, expected):
                    read_inputs(table, groups, "counts")

    def test_pairwise_tests_report_raw_p_and_skip_tiny_groups(self):
        alpha = alpha_table(read_inputs(self.abund, self.groups, "counts"), "counts")
        left = alpha.loc[alpha["group"] == "对照组", "Shannon"]
        right = alpha.loc[alpha["group"] == "藜麦组", "Shannon"]
        expected = {
            "welch": stats.ttest_ind(left, right, equal_var=False).pvalue,
            "mannwhitney": stats.mannwhitneyu(left, right, alternative="two-sided").pvalue,
        }
        for method, target in expected.items():
            with self.subTest(method=method):
                p, caption = compare_groups(alpha, "Shannon", "对照组", "藜麦组", method)
                self.assertAlmostEqual(p, target)
                self.assertIn("未做多重比较校正", caption)
                self.assertIn("对照组", caption)
                self.assertIn("藜麦组", caption)
                tiny = alpha_table(read_inputs(SMALL_COUNTS, SMALL_GROUPS, "counts"), "counts")
                skipped, reason = compare_groups(tiny, "Shannon", "甲", "乙", method)
                self.assertIsNone(skipped)
                self.assertTrue(reason)


if __name__ == "__main__":
    unittest.main()
