# Species annotation

The supplied expression matrix and the analysis annotation are separate. PCA, differential analysis, trajectories and clinical models use the feature IDs retained after the selected cleaning steps. Changing the annotation settings requires rebuilding the study. Original uploaded matrix values and metadata remain downloadable.

Feature-by-feature concordance and the Integration effect-comparison module require matching native feature IDs between layers. DIABLO instead requires matched sample IDs and can use different feature sets in each layer. Annotation supports gene-set analyses and cross-study feature comparison; it does not rename expression rows. Cross-layer pathway scores can use each layer's separately declared annotation.

## Built-in and automatic mappings

**Use built-in mapping** is the default in each layer's annotation settings. The resource hint updates with the source organism and identifier type. Built-in gene-set collections do not require a GMT upload.

The kidney atlas carries upstream human-symbol-compatible annotation. The historical `ortholog_map.csv` remains the source of strict human-pig atlas comparisons; it is restricted to that atlas, not a general ID converter.

New pig, mouse and rat uploads use `species_orthologs.rds`, derived from the complete three-species records in **babelgene 22.9 / HGNC HCOP**. Symbol, Ensembl gene and Entrez identifiers are supported. Human identifier lookups use the installed **org.Hs.eg.db** version. Ensembl version suffixes are removed for lookup only. UniProt automatic lookup is human-only.

A non-human pair is eligible only if the complete species snapshot has one distinct Entrez-gene partner in each direction and the record has at least three supporting resources. Relationship counts are calculated before filtering to the uploaded genes or applying the support threshold. This is uniqueness within the HCOP snapshot, not an Ensembl-defined homology class. No highest-scoring candidate is chosen among ambiguous pairs. Multiple measured rows mapping to the same human symbol are excluded from the annotated view instead of being counted as separate genes. Original rows remain available for feature-level analysis.

The build recipe is `tools/build_orthologs.R`; babelgene is only needed to regenerate the bundled snapshot, not to run the app. Source/version, rule and SHA-256 are recorded in `species_orthologs_provenance.json`; attribution is in `species_orthologs_LICENSE.txt`. This versioned resource is not a live or complete annotation for every current gene identifier.

## Custom mapping to human

Select the source organism (Other requires its scientific name), choose the input ID type and **Upload a custom mapping**, then supply a CSV or TSV with:

```csv
source_id,target_id
my_feature_1,TP53
my_feature_2,EGFR
```

`source_id` must exactly match the retained expression feature ID. `target_id` must be a human gene symbol compatible with the selected collection. Supply the mapping source and version. If a `relationship` column is included, only `one_to_one` rows are eligible. The app also checks bidirectional uniqueness across the complete supplied mapping, including rows absent from the expression matrix. It does not independently establish biological orthology for a user-supplied table. Multiple-target mappings must be resolved upstream from appropriate biological evidence, not from the observed significance of features.

## Native-ID GMT

For a species or assay without a suitable human mapping, choose **Native IDs only**. Open **Advanced: custom gene sets**, check **Use a custom GMT** and supply a GMT. Each tab-separated line contains a unique set name, a description (or URL), then its member feature IDs:

```text
pathway_A<TAB>source/version<TAB>feature_1<TAB>feature_2<TAB>feature_3
```

Use actual tabs, not the literal `<TAB>` text. Members refer to original feature IDs after cleaning, even when the same layer also has a human mapping. Select **Custom GMT** in enrichment, pathway activity or related modules. Cross-study custom-GMT comparison requires the same declared organism and identical set definitions. Human collections remain separate from the custom collection; one upload never changes another layer's gene sets.

Switching to built-in mapping ignores any previously selected custom mapping file on the next build. Unchecking **Use a custom GMT** likewise excludes that GMT. Selecting a custom option without its required file reports an error and leaves the existing study unchanged.

## Statistical universe and exports

Human collections are projected onto the eligible original IDs. GSEA ranks only eligible mapped features; ORA uses eligible measured features as its background, including those not assigned to a pathway. Native GMT uses all retained measured features as its background. Set-size filters still apply; sparse mapping does not relax them. Native-species gene sets and ortholog-projected human sets answer different annotation questions.

Downloads contain `mapping.csv` with source/query IDs, candidate target IDs, status, relationship and source/version. Status distinguishes mapped, unmapped, ambiguous, duplicate target and insufficient support. The supplied custom table is retained separately. `annotation.json` records settings; `custom.gmt`, input matrices and metadata are included when present. Matrices are exported as CSV values, not as byte-identical copies of the original formatting.

References: [babelgene documentation](https://igordot.github.io/babelgene/articles/babelgene-intro.html), [MSigDB guidance on species and identifier conversion](https://docs.gsea-msigdb.org/MSigDB/MSigDB_FAQ/).
