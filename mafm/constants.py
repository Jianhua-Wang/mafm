"""Define constants used in the package."""


class ColName:
    """Define column names."""

    # mandatory columns
    CHR = "CHR"
    BP = "BP"
    RSID = "rsID"
    EA = "EA"
    NEA = "NEA"
    P = "P"
    BETA = "BETA"
    SE = "SE"
    EAF = "EAF"

    # optional columns
    MAF = "MAF"
    N = "N"
    Z = "Z"
    INFO = "INFO"

    # columns for loci
    START = "START"
    END = "END"
    LEAD_SNP = "LEAD_SNP"
    LEAD_SNP_P = "LEAD_SNP_P"
    LEAD_SNP_BP = "LEAD_SNP_BP"

    # unique snpid, chr-bp-sorted(EA,NEA)
    SNPID = "SNPID"

    # COJO results
    COJO_P = "COJO_P"
    COJO_BETA = "COJO_BETA"
    COJO_SE = "COJO_SE"

    # posterior probability
    PP_FINEMAP = "PP_FINEMAP"
    PP_ABF = "PP_ABF"
    PP_PAINTOR = "PP_PAINTOR"
    PP_CAVIARBF = "PP_CAVIARBF"
    PP_SUSIE = "PP_SUSIE"
    PP_POLYFUN_FINEMAP = "PP_POLYFUN_FINEMAP"
    PP_POLYFUN_SUSIE = "PP_POLYFUN_SUSIE"

    # ordered columns
    mandatory_cols = [CHR, BP, EA, NEA, EAF, BETA, SE, P]
    sumstat_cols = [CHR, BP, RSID, EA, NEA, P, BETA, SE, EAF, MAF]
    loci_cols = [CHR, START, END, LEAD_SNP, LEAD_SNP_P, LEAD_SNP_BP]


# only support autosomes
CHROM_LIST = [i for i in range(1, 24)]
