"""
q2_taxa_clean/_methods.py

Operations logic for cleaning QIIME2 taxonomy strings.
These functions are registered as QIIME2 Methods in plugin_setup.py.
"""

import re
import pandas as pd


SUFFIXES_TO_REMOVE = [
    "_unclassified", "_uncultured", "_bacterium", "_sp.",
    "_strain", "_group", "_cluster", "_clade", "_lineage",
    " unclassified", " uncultured", " bacterium", " sp.",
    " uncultured_organism", "_uncultured_organism",
    " uncultured organism",
    " gut_metagenome", "_gut_metagenome",
    " human_gut", "_human_gut",
]

# Complete word names that are uninformative after prefix stripping
UNINFORMATIVE_WORDS = {
    "uncultured", "unclassified", "bacterium", "organism",
    "metagenome", "environmental", "clone", "sp", "spp",
    "unknown", "unidentified", "rumen", "gut", "human_gut",
    "uncultured_organism", "gut_metagenome",
}

# Names starting with these prefixes are uninformative at that level.
# Rather than stripping the prefix and returning the remainder, the whole
# level is rejected and the walker moves up to the next rank. This preserves
# the integrity of classifier decisions and allows for clean informative labels for vis, 
# if SILVA wrote uncultured_Ruminococcus at species level, the resolution is genus-level Ruminococcus from g__

UNINFORMATIVE_PREFIXES = (
    "uncultured_",
    "uncultured ",
    "unclassified_",
    "unclassified ",
    "candidatus_",
    "candidatus ",
)


def _clean_name(raw: str) -> str:
    """Strip taxonomic prefix (g__, s__ etc.) and uninformative suffixes."""
    name = re.sub(r"^[a-z]__", "", raw).strip()
    for suffix in SUFFIXES_TO_REMOVE:
        name = re.sub(re.escape(suffix), "", name, flags=re.IGNORECASE)
    return name.strip("_").strip()


def _is_informative(name: str) -> bool:
    """Return True if a cleaned name carries taxonomic clarity.

    Rejects names that are empty, purely numeric, single characters,
    bare uninformative words, or prefixed with qualifiers that indicate
    the classifier could not resolve this level (i.e. uncultured_Ruminococcus).
    In those prefix cases the whole level is rejected so the walker moves up
    to the next rank rather than reconstructing a name from a failed assignment.
    """
    if not name or name.isdigit() or len(name) <= 1:
        return False
    if not re.search(r"[a-zA-Z]", name):
        return False
    if name.lower() in UNINFORMATIVE_WORDS:
        return False
    if any(name.lower().startswith(p) for p in UNINFORMATIVE_PREFIXES):
        return False
    return True


def _get_terminal_name(taxonomy_string: str, max_level: int = 7) -> str:
    """
    Return the most specific informative name from a QIIME2 taxonomy string.

    Walks from the terminal level back toward root until an informative
    name is found. Falls back to the raw string if nothing useful exists.

    If the most specific informative level is species (s__),
    appends the genus name to give proper genus species nomenclature e.g.
    "Lactobacillus reuteri" rather than the epithet alone. If no informative
    genus exists, the species level is skipped and the walker continues up
    to the nearest clean higher rank. A bare species epithet is never
    returned. Name assembly is skipped when genus and species names are identical (i.e.
    "Ruminococcus Ruminococcus" which occurs when a genus name appears in the species slot).
    """
    parts = [p.strip() for p in taxonomy_string.split(";")] \
        if ";" in taxonomy_string else [taxonomy_string.strip()]

    parts = parts[:max_level]

    for i, part in enumerate(reversed(parts)):
        name = _clean_name(part)
        if _is_informative(name):
            original_idx = len(parts) - 1 - i
            is_species = parts[original_idx].strip().startswith("s__")

            if is_species and max_level >= 6 and original_idx >= 1:
                genus_raw = parts[original_idx - 1].strip()
                genus = _clean_name(genus_raw)
                if _is_informative(genus) and genus.lower() != name.lower():
                    return f"{genus} {name}"
                else:
                    continue  # genus uncertain, walk up

            return name  # only reached for non-species levels

    return taxonomy_string  # fallback


def _disambiguate(index: list) -> list:
    """
    Append _1, _2, _3 to ALL occurrences of duplicate names.
    First occurrence also gets a suffix for consistency.
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
    QIIME2 Method: clean taxonomy strings to their most specific named level.

    Takes a pd.Series of taxonomy strings (as returned from
    FeatureData[Taxonomy]) and returns a cleaned Series.

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
