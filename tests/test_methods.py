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
    _get_truncated_taxonomy,   
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

    def test_stacked_suffixes(self):
        assert _clean_name("f__Lachnospiraceae_unclassified_group") == "Lachnospiraceae"


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

    def test_uncultured_prefix_is_not_informative(self):
        assert _is_informative("uncultured_Ruminococcus") is False

    def test_unclassified_prefix_is_not_informative(self):
        assert _is_informative("unclassified_Lachnospiraceae") is False

    def test_candidatus_prefix_is_not_informative(self):
        assert _is_informative("candidatus_Saccharimonas") is False

    def test_bracketed_name_is_informative(self):
        # brackets are SILVA reclassification flags, not uninformative prefixes
        assert _is_informative("[Clostridium]_innocuum") is True


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

    def test_uncultured_species_uncertain_genus_walks_to_family(self):
        # uncultured_ prefix on both genus and species — should walk up to family
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__uncultured_Ruminococcus;s__uncultured_sp"
        assert _get_terminal_name(s) == "Ruminococcaceae"

    def test_valid_species_uncertain_genus_walks_to_family(self):
        # species epithet is real but genus is uncertain — bare epithet must never be returned
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__uncultured_Ruminococcus;s__reuteri"
        assert _get_terminal_name(s) == "Ruminococcaceae"

    def test_genus_equals_species_returns_genus(self):
        # Ruminococcus in species slot, identity check fires continue, walker resolves to genus
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Ruminococcus;s__Ruminococcus"
        assert _get_terminal_name(s) == "Ruminococcus"


####### ── _get_truncated_taxonomy ── ########

class TestGetTruncatedTaxonomy:
    def test_standard_full_string_preserves_path(self):
        # full informative string — should come back intact to species level
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"
        assert _get_truncated_taxonomy(s) == "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__reuteri"

    def test_empty_species_truncates_at_genus(self):
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__"
        assert _get_truncated_taxonomy(s) == "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus"

    def test_uncultured_genus_truncates_at_family(self):
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__uncultured;s__uncultured_bacterium"
        assert _get_truncated_taxonomy(s) == "d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae"

    def test_valid_species_good_genus_returns_full_path_to_species(self):
        # both genus and species informative — full path preserved to s__ level
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"
        result = _get_truncated_taxonomy(s)
        assert result.endswith("s__reuteri")

    def test_uncertain_genus_skips_species_truncates_at_family(self):
        # species present but genus is uncertain — bare epithet must never be returned
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__uncultured_Ruminococcus;s__reuteri"
        assert _get_truncated_taxonomy(s) == "d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae"

    def test_genus_equals_species_truncates_at_genus(self):
        # Ruminococcus in species slot — identity check fires continue, truncates at g__ level
        s = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Ruminococcus;s__Ruminococcus"
        assert _get_truncated_taxonomy(s) == "d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Ruminococcus"

    def test_max_level_truncates_at_genus(self):
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"
        result = _get_truncated_taxonomy(s, max_level=6)
        assert result == "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus"

    def test_all_empty_returns_raw(self):
        s = "d__;p__;c__;o__;f__;g__;s__"
        assert _get_truncated_taxonomy(s) == s

    def test_no_semicolons(self):
        # single-level string, no semicolons — branch handled separately in function
        assert _get_truncated_taxonomy("g__Bacteroides") == "g__Bacteroides"


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

    def test_default_returns_truncated_strings(self):
        taxonomy = pd.Series(
            {
                "feat1": "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri",
                "feat2": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__",
                "feat3": "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__uncultured;s__uncultured_bacterium",
            }
        )
        result = clean_taxonomy(taxonomy)
        assert result["feat1"] == "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__reuteri"
        assert result["feat2"] == "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides"
        assert result["feat3"] == "d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae"

  
    def test_flat_labels_returns_readable_names(self):
        taxonomy = pd.Series(
            {
                "feat1": "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri",
                "feat2": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__",
                "feat3": "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__uncultured;s__uncultured_bacterium",
            }
        )
        result = clean_taxonomy(taxonomy, flat_labels=True)
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
                "feat2": "d__Bacteria;p__Firmicutes;g__uncultured;s__uncultured_bacterium",
            }
        )
        result = clean_taxonomy(taxonomy, flat_labels=True)  # FIXED: flat_labels=True required for disambiguation
        assert result["feat1"] == "Firmicutes_1"
        assert result["feat2"] == "Firmicutes_2"

    
    def test_max_level_passed_through(self):
        taxonomy = pd.Series(
            {"feat1": "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__reuteri"}
        )
        result = clean_taxonomy(taxonomy, max_level=6)
        assert result["feat1"] == "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus"
