q2-taxa-clean
A QIIME 2 plugin for automatically cleaning and simplifying taxonomy strings.

After running taxonomy classification in QIIME 2, feature tables commonly contain uninformative labels like:
d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__uncultured;s__uncultured_bacterium

q2-taxa-clean automatically resolves each feature to its most specific informative name by walking back up the taxonomy hierarchy. 

The above example becomes: Rhizobiaceae

How it works
For each feature, the plugin:

1. Splits the taxonomy string by ;
2. Starts at the terminal level (e.g. species) and walks back toward the root
3. Strips taxonomic prefixes (g__, s__ etc.) and uninformative labels (uncultured, unclassified, bacterium etc.)
4. Returns the first informative name it finds
5. Disambiguates any duplicate names with _1, _2, _3 suffixes

Installation
Activate your QIIME 2 conda environment first, then:
pip install git+https://github.com/ssaundy/q2-taxa-clean.git

Or for development:
git clone https://github.com/ssaundy/q2-taxa-clean.git
cd q2-taxa-clean
pip install -e . --no-deps
