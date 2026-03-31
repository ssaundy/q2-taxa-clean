"""
q2_taxa_clean/_methods.py

Logic for cleaning QIIME2 taxonomy strings.
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

# Complete words that are uninformative after prefix stripping
UNINFORMATIVE_WORDS = {
    "uncultured", "unclassified", "bacterium", "organism",
    "metagenome", "environmental", "clone", "sp", "spp",
    "unknown", "unidentified", "rumen", "gut", "human_gut",
    "uncultured_organism", "gut_metagenome",
}

# Names starting with these prefixes are uninformative at that level.
# Rather than stripping the prefix and returning the remainder, the whole
# level is rejected and the walker moves up to the next rank. This preserves
# the integrity of classifier decisions, if SILVA wrote g__uncultured_Ruminococcus,
# genus level is rejected and the resolution is f__Ruminococcaceae

UNINFORMATIVE_PREFIXES = (
    "uncultured_",
    "uncultured ",
    "unclassified_",
    "unclassified ",
    "candidatus_",
    "candidatus ",
)


def _clean_name(raw: str) -> str:
    """Strip taxonomic prefix (g__, s__ etc.) and uninformative suffixes.

    Suffixes are anchored to end of the string ($) so they only fire at the
    terminal position, preventing mid name substring matches e.g. _group must not
    fire inside a name that happens to contain it elsewhere.

    To avoid the fragility stripping in a single pass we iterate, removing
    terminal suffixes until the name stabilises. Handles stacked
    suffixes such eg _unclassified_group regardless of the
    order entries appear in SUFFIXES_TO_REMOVE.
    """
    name = re.sub(r"^[a-z]__", "", raw).strip()
    changed = True
    while changed:
        prev = name
        for suffix in SUFFIXES_TO_REMOVE:
            name = re.sub(re.escape(suffix) + r"$", "", name, flags=re.IGNORECASE)
        changed = name != prev
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
    Returns the most specific informative name from a QIIME2 taxonomy string.

    Walks from the terminal level back towards the root until an informative
    name is found. Falls back to the raw string if nothing useful exists.

    If the most specific informative level is species (s__),
    appends the genus name to give proper genus species nomenclature e.g.
    "Lactobacillus reuteri" rather than the epithet alone. If no informative
    genus exists, the species level is skipped and the walker continues up
    to the nearest clean higher rank. A bare species epithet is never
    returned. Name assembly is skipped when genus and species names are identical (i.e.
    "Ruminococcus Ruminococcus" which occurs when a genus name appears in the species slot).

    ***called as a post-processing step only when flat_labels=True
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
                    continue  # genus uncertain or matches species name, return genus on next iteration

            return name  # only reached for non-species levels

    return taxonomy_string  # fallback



def _get_truncated_taxonomy(taxonomy_string: str, max_level: int = 7) -> str:
    """
    Return the taxonomy string truncated at the most specific informative level.

    Walks from the terminal level back toward root using the same rules as
    _get_terminal_name, but returns the full ancestry path cut to the
    informative level rather than a single extracted name.

    The returned string is valid FeatureData[Taxonomy] format and interoperable
    with core downstream QIIME 2 tools.

    Falls back to the raw string if nothing useful exists.

    Params
    ---------
    taxonomy_string : str
        A single QIIME2 taxonomy string, semicolon-delimited.
    max_level : int
        Maximum taxonomy depth to consider (1=domain, 7=species).


    Returns
    -----------
    str
        Taxonomy string truncated at the deepest informative level,
        with levels separated by "; ".
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
                    # species and genus both informative — include both (full string to s__ level)
                    return "; ".join(parts[:original_idx + 1])
                else:
                    continue  # genus uncertain, walk up — do not return bare epithet

            # non-species level, or species was skipped — return up to and including this level
            return "; ".join(parts[:original_idx + 1])

    return taxonomy_string  # fallback: return raw string unchanged



def _disambiguate(index: list) -> list:
    """
    Append _1, _2, _3 to ALL occurrences of duplicate names.
    First occurrence also gets a suffix for consistency.
    Only called when flat_labels=True.
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
    flat_labels: bool = False,
) -> pd.Series:
    """
    QIIME2 Method: clean taxonomy strings to their most specific named level.

    Takes a pd.Series of taxonomy strings (as returned from
    FeatureData[Taxonomy]) and returns a cleaned Series.

    By default returns truncated taxonomy strings with full hierarchical paths
    preserved, valid FeatureData[Taxonomy] format, full compatability with downstream
    QIIME 2 tools.

    When flat_labels=True, returns simple readable names without the full path, suitable for quick figure
    labels and visualisation, but invalid for tools expecting FeatureData[Taxonomy] format specifically.

    Parameters
    ----------
    taxonomy : pd.Series
        Taxonomy strings indexed by feature ID.
        e.g. "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;
               f__Lactobacillaceae;g__Lactobacillus;s__"
    max_level : int
        Maximum taxonomy depth to consider (1=domain, 7=species).
    flat_labels : bool
        If False (default), return truncated taxonomy strings.
        If True, return correct nomenclature taxon names (e.g. "Lactobacillus reuteri").

    Returns
    -------
    pd.Series
        Cleaned taxonomy strings indexed by feature ID.
    """
    if flat_labels:
        # Derive flat labels from the truncated string, not a separate walking pass.
        # This guarantees flat labels are consistent with truncated output since
        # _get_terminal_name applied to a pre truncated string will always
        # return the deepest level name.
        truncated = taxonomy.apply(
            lambda x: _get_truncated_taxonomy(str(x), max_level=max_level)
        )
        cleaned = truncated.apply(
            lambda x: _get_terminal_name(str(x), max_level=max_level)
        )
        # Disambiguate duplicate flat labels since two different features can resolve
        # to the same readable name (i.e. both become "Ruminococcaceae").
        ## Not applied to truncated strings where identical output = identical taxon.
        disambiguated = _disambiguate(cleaned.tolist())
        return pd.Series(disambiguated, index=taxonomy.index, name=taxonomy.name)
    else:
        cleaned = taxonomy.apply(
            lambda x: _get_truncated_taxonomy(str(x), max_level=max_level)
        )
        return pd.Series(cleaned.tolist(), index=taxonomy.index, name=taxonomy.name)
