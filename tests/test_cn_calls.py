import math
import unittest

from facetsAPI.cn_calls import CN_CALL_TABLE, get_cn_call


class CnCallTableTests(unittest.TestCase):
    def test_lcn_keyed_basics(self):
        self.assertEqual(get_cn_call(False, 2, 1), "DIPLOID")
        self.assertEqual(get_cn_call(False, 2, 0), "CNLOH")
        self.assertEqual(get_cn_call(False, 1, 0), "HETLOSS")
        self.assertEqual(get_cn_call(False, 0, 0), "HOMDEL")
        self.assertEqual(get_cn_call(False, 3, 1), "GAIN")
        self.assertEqual(get_cn_call(False, 4, 2), "TETRAPLOID")
        self.assertEqual(get_cn_call(False, 6, 3), "AMP (BALANCED)")
        self.assertEqual(get_cn_call(True, 4, 2), "TETRAPLOID")
        self.assertEqual(get_cn_call(True, 2, 1), "DOUBLE LOSS AFTER")

    def test_impossible_lcn_absent(self):
        # lcn cannot exceed tcn/2; the old table had these as (tcn, major) keys
        self.assertIsNone(get_cn_call(False, 1, 1))
        self.assertIsNone(get_cn_call(False, 3, 2))
        self.assertNotIn((False, 1, 1), CN_CALL_TABLE)

    def test_na_handling(self):
        self.assertEqual(get_cn_call(False, 1, None), "HETLOSS")
        self.assertEqual(get_cn_call(False, 1, float("nan")), "HETLOSS")
        self.assertEqual(get_cn_call(False, 1, -1), "HETLOSS")   # FacetsSegment NA sentinel
        self.assertIsNone(get_cn_call(False, 2, None))           # ambiguous by design
        self.assertIsNone(get_cn_call(False, float("nan"), 1))
        self.assertIsNone(get_cn_call(False, -1, -1))

    def test_float_inputs(self):
        self.assertEqual(get_cn_call(0, 2.0, 1.0), "DIPLOID")
        self.assertEqual(get_cn_call("", 4.0, 2.0), "TETRAPLOID")

    def test_every_key_is_lcn_consistent(self):
        for (wgd, tcn, lcn) in CN_CALL_TABLE:
            self.assertIsInstance(wgd, bool)
            if lcn is not None:
                self.assertLessEqual(lcn * 2, tcn, msg=str((wgd, tcn, lcn)))


if __name__ == "__main__":
    unittest.main()
