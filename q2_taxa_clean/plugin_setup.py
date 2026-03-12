"""
q2_taxa_clean/plugin_setup.py

Register q2-taxa-clean as a QIIME 2 plugin.

Plugin exposes one Method:
    qiime taxa-clean clean-taxonomy
        --i-taxonomy  FeatureData[Taxonomy]
        --p-max-level INT (1-7, default 7)
        --o-cleaned-taxonomy FeatureData[Taxonomy]
"""

import importlib

from qiime2.plugin import (
    Plugin,
    Int,
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
        "Automatically resolves uninformative terminal levels (e.g. g__, "
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
    },
    outputs=[
        ("cleaned_taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={
        "taxonomy": (
            "Taxonomy strings to clean, as produced by qiime feature-classifier "
            "classify-sklearn or similar. Each string should be semicolon-delimited "
            "with standard QIIME 2 prefixes (d__, p__, c__, o__, f__, g__, s__)."
        ),
    },
    parameter_descriptions={
        "max_level": (
            "Maximum taxonomy depth to consider when searching for the most "
            "specific informative name. 1=domain, 2=phylum, 3=class, 4=order, "
            "5=family, 6=genus, 7=species. Default is 7 (species level)."
        ),
    },
    output_descriptions={
        "cleaned_taxonomy": (
            "Cleaned taxonomy strings where each feature is labelled with its "
            "most specific informative taxonomic name. Uninformative levels "
            "(empty, uncultured, unclassified) are resolved by walking back up "
            "the hierarchy. Duplicate names are disambiguated with _1, _2 suffixes."
        ),
    },
    name="Clean taxonomy strings",
    description=(
        "Automatically clean QIIME 2 taxonomy strings by resolving each feature "
        "to its most specific informative taxonomic name. Unlike qiime taxa collapse "
        "(which merges features at a fixed level) or rescript edit-taxonomy "
        "(which requires manual find-and-replace strings), this method automatically "
        "walks the taxonomy hierarchy from the terminal level upward to find the "
        "most specific named level for each feature individually. "
        "Handles uninformative labels such as g__, s__uncultured_bacterium, "
        "and s__uncultured_organism."
    ),
)
