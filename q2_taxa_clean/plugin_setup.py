"""
q2_taxa_clean/plugin_setup.py
Register q2-taxa-clean as a QIIME 2 plugin.

Plugin exposes one Method:

    qiime taxa-clean clean-taxonomy
        --i-taxonomy  FeatureData[Taxonomy]
        --p-max-level INT (1-7, default 7)
        --p-flat-labels BOOL (default False)
        --o-cleaned-taxonomy FeatureData[Taxonomy]

"""
import importlib
from qiime2.plugin import (
    Plugin,
    Int,
    Bool,          
    Range,
    Citations,
)
from q2_types.feature_data import FeatureData, Taxonomy
import q2_taxa_clean
from q2_taxa_clean._methods import clean_taxonomy

######### ── Plugin object ── ############

plugin = Plugin(
    name="taxa-clean",
    version=q2_taxa_clean.__version__,
    website="https://github.com/ssaundy/q2-taxa-clean",
    package="q2_taxa_clean",
    description=(
        "A QIIME 2 plugin for cleaning and simplifying taxonomy strings. "
        "Automatically resolves uninformative terminal levels (i.e. g__, "
        "s__uncultured_bacterium) by walking back up the hierarchy to the "
        "most specific named level."
    ),
    short_description="Clean and simplify QIIME 2 taxonomy strings.",
)

######## ── Register Method ─── ###########

plugin.methods.register_function(
    function=clean_taxonomy,
    inputs={
        "taxonomy": FeatureData[Taxonomy],
    },
    parameters={
        "max_level": Int % Range(1, 7, inclusive_end=True),
        # flat_labels registered as a Bool parameter.
        # Default False means the output (truncated taxonomy
        # strings) is what users get unless they explicitly opt in to flat labels.
        "flat_labels": Bool,
    },
    outputs=[
        ("cleaned_taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={
        "taxonomy": (
            "Taxonomy strings to clean, as produced by qiime feature classifier "
            "classify sklearn or similar. Each string should be semicolon delimited "
            "with standard QIIME 2 prefixes (d__, p__, c__, o__, f__, g__, s__)."
        ),
    },
    parameter_descriptions={
        "max_level": (
            "Maximum taxonomy depth to consider when searching for the most "
            "specific informative name. 1=domain, 2=phylum, 3=class, 4=order, "
            "5=family, 6=genus, 7=species. Default is 7 (species level)."
        ),
      
        "flat_labels": (
            "If False (default), return truncated taxonomy strings with the full "
            "hierarchical path preserved up to the most specific informative level "
            "output is valid FeatureData[Taxonomy] and compatable with downstream QIIME 2 "
            "tools such as taxa barplot and taxa collapse. "
            "If True, return single readable taxon names only "
            "(e.g. 'Ruminococcaceae', 'Lactobacillus reuteri'). Useful for figure "
            "legends and axis labels. Duplicate names are disambiguated with _1, _2 "
            "suffixes. Output is not valid FeatureData[Taxonomy] and should not be "
            "passed to downstream QIIME 2 tools."
        ),
    },
    output_descriptions={

        "cleaned_taxonomy": (
            "Cleaned taxonomy strings where each feature is resolved to its most "
            "specific informative taxonomic level. By default, the full hierarchical "
            "path is preserved up to that level (e.g. 'd__Bacteria; p__Firmicutes; "
            "f__Ruminococcaceae'), keeping output valid as FeatureData[Taxonomy]. "
            "When --p-flat-labels is set to True, returns single readable names "
            "with duplicate names disambiguated using _1, _2 suffixes."
        ),
    },
    name="Clean taxonomy strings",
    description=(
        "Automatically clean QIIME 2 taxonomy strings by resolving each feature "
        "to its most specific informative taxonomic level. Unlike qiime taxa collapse "
        "(which merges all features at a fixed level) or rescript edit-taxonomy "
        "(which requires manual find-and-replace strings), this method walks the "
        "taxonomy hierarchy from the terminal level upward to find the most specific "
        "named level for each feature individually, preserving the full path. "
        "Handles uninformative labels such as g__, s__uncultured_bacterium, "
        "uncultured prefixes (e.g. g__uncultured_Ruminococcus), and "
        "s__uncultured_organism. Use --p-flat-labels for single readable names "
        "suitable for figure labels and visualisation."
    ),
)
