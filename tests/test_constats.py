import unittest

from mafm.constants import CHROM_LIST, ColName


class TestConstants(unittest.TestCase):

    def test_colname_mandatory_columns(self):
        self.assertEqual(ColName.CHR, "CHR")
        self.assertEqual(ColName.BP, "BP")
        self.assertEqual(ColName.RSID, "rsID")
        self.assertEqual(ColName.EA, "EA")
        self.assertEqual(ColName.NEA, "NEA")
        self.assertEqual(ColName.P, "P")
        self.assertEqual(ColName.BETA, "BETA")
        self.assertEqual(ColName.SE, "SE")
        self.assertEqual(ColName.EAF, "EAF")

    def test_colname_optional_columns(self):
        self.assertEqual(ColName.MAF, "MAF")
        self.assertEqual(ColName.N, "N")
        self.assertEqual(ColName.Z, "Z")
        self.assertEqual(ColName.INFO, "INFO")

    def test_colname_loci_columns(self):
        self.assertEqual(ColName.START, "START")
        self.assertEqual(ColName.END, "END")
        self.assertEqual(ColName.LEAD_SNP, "LEAD_SNP")
        self.assertEqual(ColName.LEAD_SNP_P, "LEAD_SNP_P")
        self.assertEqual(ColName.LEAD_SNP_BP, "LEAD_SNP_BP")

    def test_colname_snpid_column(self):
        self.assertEqual(ColName.SNPID, "SNPID")

    def test_colname_cojo_columns(self):
        self.assertEqual(ColName.COJO_P, "COJO_P")
        self.assertEqual(ColName.COJO_BETA, "COJO_BETA")
        self.assertEqual(ColName.COJO_SE, "COJO_SE")

    def test_colname_posterior_probability_columns(self):
        self.assertEqual(ColName.PP_FINEMAP, "PP_FINEMAP")
        self.assertEqual(ColName.PP_ABF, "PP_ABF")
        self.assertEqual(ColName.PP_PAINTOR, "PP_PAINTOR")
        self.assertEqual(ColName.PP_CAVIARBF, "PP_CAVIARBF")
        self.assertEqual(ColName.PP_SUSIE, "PP_SUSIE")
        self.assertEqual(ColName.PP_POLYFUN_FINEMAP, "PP_POLYFUN_FINEMAP")
        self.assertEqual(ColName.PP_POLYFUN_SUSIE, "PP_POLYFUN_SUSIE")

    def test_colname_ordered_columns(self):
        self.assertEqual(
            ColName.mandatory_cols, ["CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P"]
        )
        self.assertEqual(
            ColName.sumstat_cols,
            ["CHR", "BP", "rsID", "EA", "NEA", "P", "BETA", "SE", "EAF", "MAF"],
        )
        self.assertEqual(
            ColName.loci_cols, ["CHR", "START", "END", "LEAD_SNP", "LEAD_SNP_P", "LEAD_SNP_BP"]
        )

    def test_chrom_list(self):
        self.assertEqual(CHROM_LIST, list(range(1, 24)))


if __name__ == "__main__":
    unittest.main()
