"""
q2_taxa_clean/_methods.py

Core logic for cleaning QIIME2 taxonomy strings.
These functions are registered as QIIME2 Methods in plugin_setup.py.
"""

import re
import pandas as pd


SUFFIXES_TO_REMOVE = [
    "_unclassified", "_uncultured", "_bacterium", "_sp.",
    "_strain", "_group", "_cluster", "_clade", "_lineage",
    " unclassified", " uncultured", " bacterium", " sp.",
]

# complete word names that are uninformative after prefix stripping
UNINFORMATIVE_WORDS = {
    "uncultured", "unclassified", "bacterium", "organism",
    "metagenome", "environmental", "clone", "sp", "spp",
    "unknown", "unidentified",
}


def _clean_name(raw: str) -> str:
    """Strip taxonomic prefix (g__, s__ etc.) and uninformative suffixes."""
    name = re.sub(r"^[a-z]__", "", raw).strip()
    for suffix in SUFFIXES_TO_REMOVE:
        name = re.sub(re.escape(suffix), "", name, flags=re.IGNORECASE)
    return name.strip("_").strip()


def _is_informative(name: str) -> bool:
    """Return True if a cleaned name carries useful biological meaning."""
    if not name or name.isdigit() or len(name) <= 1:
        return False
    if not re.search(r"[a-zA-Z]", name):
        return False
    if name.lower() in UNINFORMATIVE_WORDS:
        return False
    return True


def _get_terminal_name(taxonomy_string: str, max_level: int = 7) -> str:
    """
    Return the most specific informative name from a QIIME2 taxonomy string

    Walks from the terminal level back toward root until an informative
    name is found. Falls back to the raw string if nothing useful exists
    """
    parts = [p.strip() for p in taxonomy_string.split(";")] \
        if ";" in taxonomy_string else [taxonomy_string.strip()]

    for part in reversed(parts[:max_level]):
        name = _clean_name(part)
        if _is_informative(name):
            return name

    return taxonomy_string  # fallback


def _disambiguate(index: list) -> list:
    """
    Append _1, _2, _3 to ALL occurrences of duplicate names
    First occurrence also gets a suffix for consistency
    """
    from collections import Counter
    duplicates = {n for n, c in Counter(index).items() if c > 1}
    occurrence = {}
    result = []
    for name in index:
        if name in duplicates:
            occurrence[name] = occurrence.get(name, 0) + 1
            result.append(f"{name}_{occurrence[name]}")
        else:
            result.append(name)
    return result


def clean_taxonomy(
    taxonomy: pd.Series,
    max_level: int = 7,
) -> pd.Series:
    """
    QIIME2 Method: clean taxonomy strings to their most specific named level

    Takes a pd.Series of taxonomy strings (as returned from
    FeatureData[Taxonomy]) and returns a cleaned Series

    Parameters
    ----------
    taxonomy : pd.Series
        Taxonomy strings indexed by feature ID.
        e.g. "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;
               f__Lactobacillaceae;g__Lactobacillus;s__"
    max_level : int
        Maximum taxonomy depth to consider (1=domain, 7=species).

    Returns
    -------
    pd.Series
        Cleaned taxonomy strings indexed by feature ID.
    """
    cleaned = taxonomy.apply(
        lambda x: _get_terminal_name(str(x), max_level=max_level)
    )
    cleaned_list = _disambiguate(cleaned.tolist())
    return pd.Series(cleaned_list, index=taxonomy.index, name=taxonomy.name)
