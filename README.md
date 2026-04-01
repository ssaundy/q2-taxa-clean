# q2-taxa-clean

A QIIME 2 plugin for automatically cleaning taxonomy strings to their most specific informative level.

---

## The problem

After running taxonomic classification in QIIME 2, feature tables commonly contain useless terminal labels:

```
d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__uncultured; s__uncultured_bacterium
```

The classifier could not resolve genus or species... but existing tools give you two bad options:

- **`qiime taxa collapse`**  flattens all features to a single fixed level, merging taxa in the process
- **`rescript edit-taxonomy`**  find and replace, manual operation, string by string

---

## What this plugin does

For each feature, the plugin walks back up the taxonomy hierarchy from the terminal level, applying consistent rules to identify the deepest level that carries taxonomic information. It returns the string truncated at that level, with the full path intact.

The example above becomes:

```
d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae
```

The output is valid `FeatureData[Taxonomy]` which is fully compatable with core Q2 tools like`qiime taxa barplot`, `qiime taxa collapse`, differential abundance tools, and any other downstream step that expects Q2 taxonomy strings.

### What counts as uninformative

The walker rejects a level if the name (after stripping the rank prefix) is:

1. Empty (`g__`, `s__`)
2. A bare uninformative word i.e. `uncultured`, `unclassified`, `bacterium`, `metagenome`, `organism`, `gut`, `rumen`, and others
3. Prefixed with a qualifier indicating the classifier failed at that level: `uncultured_`, `unclassified_`, `candidatus_` then the whole level is rejected. If SILVA wrote `g__uncultured_Ruminococcus`, L6 is rejected and the resolution is family (L5) `Ruminococcaceae` from `f__`, rather than a reconstructed name from a failed assignment.

### Species handling

If the deepest informative level is species (`s__`), the plugin checks whether the genus level is also informative. If it is, both are included in the truncated string, giving a complete binomial. If the genus is uncertain, the species level is skipped entirely and a bare species epithet is never returned.

---

## Usage

### Default: truncated taxonomy strings (recommended)

```bash
qiime taxa-clean clean-taxonomy \
  --i-taxonomy taxonomy.qza \
  --p-max-level 7 \
  --o-cleaned-taxonomy cleaned-taxonomy.qza
```

Output is valid `FeatureData[Taxonomy]`. Pass directly to `taxa barplot`, `taxa collapse`, or any downstream tool.

### Flat labels for figures and exported tables

```bash
qiime taxa-clean clean-taxonomy \
  --i-taxonomy taxonomy.qza \
  --p-max-level 7 \
  --p-flat-labels True \
  --o-cleaned-taxonomy flat-labels-taxonomy.qza
```

Returns single readable names (`Ruminococcaceae`, `Lactobacillus reuteri`). Useful for axis labels and legend entries in R/Python figures. Duplicate names are disambiguated with `_1`, `_2` suffixes. **Not valid `FeatureData[Taxonomy]`** therefore do not pass to downstream QIIME 2 tools.

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--p-max-level` | Int (1–7) | 7 | Maximum taxonomy depth to consider. 1=domain, 7=species. |
| `--p-flat-labels` | Bool | False | Return single readable names instead of truncated strings. |

---

## Installation

Activate your QIIME 2 conda environment first, then:

```bash
pip install git+https://github.com/ssaundy/q2-taxa-clean.git
```

For development:

```bash
git clone https://github.com/ssaundy/q2-taxa-clean.git
cd q2-taxa-clean
pip install -e . --no-deps
```

---

## How it differs from existing tools

| Tool | Behaviour |
|---|---|
| `qiime taxa collapse` | Collapses all features to one fixed rank, merging distinct taxa |
| `rescript edit-taxonomy` | Manual find-and-replace, one pattern at a time |
| **q2-taxa-clean** | Walks each feature individually to its deepest informative rank, preserving the full path |
