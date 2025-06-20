from bio_data_explorer import vcf_parser as vcf_p


def test_is_low_confidence_variant_if_filter_not_none():
    assert (
        vcf_p._is_low_confidence_variant("LowQual;DepthFilter", 80.0, 40, "identifier")
        == True
    )


def test_is_low_confidence_variant_if_qual_below_threshold():
    assert vcf_p._is_low_confidence_variant("None", 10.0, 40, "identifier") == True


def test_is_low_confidence_variant_if_info_dp_below_threshold():
    assert vcf_p._is_low_confidence_variant("None", 50.0, 10, "identifier") == True


def test_is_low_confidence_variant_returns_false_for_passing_variants():
    assert vcf_p._is_low_confidence_variant("None", 50.0, 50, "identifier") == False
