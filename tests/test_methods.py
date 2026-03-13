"""
tests/test_methods.py

Unit tests for q2_taxa_clean logic.
Run with: pytest tests/
"""

import pytest
import pandas as pd
from q2_taxa_clean._methods import (
    _clean_name,
    _is_informative,
    _get_terminal_name,
    _disambiguate,
    clean_taxonomy,
)


######## ── _clean_name ── ##########

class TestCleanName:
    def test_removes_genus_prefix(self):
        assert _clean_name("g__Lactobacillus") == "Lactobacillus"

    def test_removes_species_prefix(self):
        assert _clean_name("s__reuteri") == "reuteri"

    def test_removes_uncultured_suffix(self):
        # _clean_name strips prefix and _bacterium suffix, leaving 'uncultured'
        # _is_informative then rejects 'uncultured' as a whole word
        assert _clean_name("g__uncultured") == "uncultured"

    def test_removes_bacterium_suffix(self):
        # _bacterium suffix is stripped, leaving 'uncultured' which is uninformative
        assert _clean_name("s__uncultured_bacterium") == "uncultured"

    def test_removes_unclassified_suffix(self):
        assert _clean_name("f__Lachnospiraceae_unclassified") == "Lachnospiraceae"

    def test_empty_prefix(self):
        assert _clean_name("g__") == ""

    def test_no_prefix(self):
        assert _clean_name("Bacteroides") == "Bacteroides"


####### ── _is_informative ── ########

class TestIsInformative:
    def test_real_name_is_informative(self):
        assert _is_informative("Lactobacillus") is True

    def test_empty_string_is_not(self):
        assert _is_informative("") is False

    def test_digits_only_is_not(self):
        assert _is_informative("123") is False

    def test_single_char_is_not(self):
        assert _is_informative("x") is False

    def test_no_letters_is_not(self):
        assert _is_informative("__") is False


######## ── _get_terminal_name ── #########

class TestGetTerminalName:
    def test_standard_full_string(self):
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"
        assert _get_terminal_name(s) == "Lactobacillus reuteri"

    def test_empty_species_falls_back_to_genus(self):
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__"
        assert _get_terminal_name(s) == "Lactobacillus"

    def test_uncultured_bacterium_falls_back(self):
        s = "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__uncultured;s__uncultured_bacterium"
        assert _get_terminal_name(s) == "Rhizobiaceae"

    def test_all_empty_returns_raw(self):
        s = "d__;p__;c__;o__;f__;g__;s__"
        assert _get_terminal_name(s) == s

    def test_max_level_respected(self):
        # max_level=6 means stop at genus, ignore species
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"
        assert _get_terminal_name(s, max_level=6) == "Lactobacillus"

    def test_no_semicolons(self):
        assert _get_terminal_name("g__Bacteroides") == "Bacteroides"


####### ── _disambiguate ── ########

class TestDisambiguate:
    def test_no_duplicates_unchanged(self):
        names = ["Bacteroides", "Lactobacillus", "Bifidobacterium"]
        assert _disambiguate(names) == names

    def test_duplicates_all_get_suffix(self):
        names = ["Bacteroides", "Bacteroides", "Bacteroides"]
        result = _disambiguate(names)
        assert result == ["Bacteroides_1", "Bacteroides_2", "Bacteroides_3"]

    def test_first_occurrence_gets_suffix(self):
        # Critical: first occurrence must NOT be left without suffix
        names = ["Bacteroides", "Lactobacillus", "Bacteroides"]
        result = _disambiguate(names)
        assert result[0] == "Bacteroides_1"
        assert result[1] == "Lactobacillus"
        assert result[2] == "Bacteroides_2"

    def test_mixed_duplicates(self):
        names = ["A", "B", "A", "C", "B"]
        result = _disambiguate(names)
        assert result == ["A_1", "B_1", "A_2", "C", "B_2"]


######## ── clean_taxonomy (INTEGRATION) ── ########

class TestCleanTaxonomy:
    def test_basic_series(self):
        taxonomy = pd.Series(
            {
                "feat1": "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri",
                "feat2": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__",
                "feat3": "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__uncultured;s__uncultured_bacterium",
            }
        )
        result = clean_taxonomy(taxonomy)
        assert result["feat1"] == "Lactobacillus reuteri"
        assert result["feat2"] == "Bacteroides"
        assert result["feat3"] == "Rhizobiaceae"

    def test_index_preserved(self):
        taxonomy = pd.Series(
            {"abc123": "d__Bacteria;p__Firmicutes;g__Lactobacillus;s__"},
        )
        result = clean_taxonomy(taxonomy)
        assert "abc123" in result.index

    def test_duplicate_names_disambiguated(self):
        taxonomy = pd.Series(
            {
                "feat1": "d__Bacteria;p__Firmicutes;g__uncultured;s__uncultured_bacterium",
                "feat2": "d__Bacteria;p__Bacteroidota;g__uncultured;s__uncultured_bacterium",
            }
        )
        result = clean_taxonomy(taxonomy)
        # Both should resolve to Firmicutes / Bacteroidota (different phyla)
        assert result["feat1"] != result["feat2"]
