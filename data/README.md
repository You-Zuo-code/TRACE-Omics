# Data sources

## Kidney perfusion reference

`atlas.rds` contains processed transcriptome and proteome matrices, sample metadata and precomputed annotations for 30 pigs (60 paired kidneys) and 24 human donor kidneys. The porcine series has six visits including an in vivo baseline; the human series has five perfusion visits and 33 clinical variables.

These data accompany Lantinga et al., published online in Kidney International on 22 September 2026, [doi:10.1016/j.kint.2026.08.017](https://www.kidney-international.org/article/S0085-2538(26)00792-1/fulltext). The animal study was approved by the Danish Animal Experiments Inspectorate, reference 2020-15-0201-00624. The human material consists of deceased-donor kidneys discarded for transplantation and made available for research in the United States under the Uniform Anatomical Gift Act.

RNA values are TPM, analysed as log2(TPM+1). Protein values are on the supplied log2 intensity scale. The software examples retain the matrices provided with the application; they are not a rerun of the source paper's complete processing and exclusion pipeline.

The software's MIT licence does not apply to the biological data. Cite the companion study when using them.

## Upload examples

The upload page offers human mRNA alone or matched human mRNA and protein, both copied from the curated reference. There are no external GEO or T-cell examples loaded by these buttons. External validation inputs accompany the separate manuscript-reproduction materials.

## Auxiliary tables

The bundled pathway collections match MSigDB 2025.1.Hs. See `genesets_provenance.json` for exact membership fingerprints and `genesets_LICENSE.txt` for collection-specific terms. These resources are not covered by the software's MIT licence.

`ortholog_map.csv` and `ortholog_map_provenance.json` record the supplied human-pig map from babelgene 22.9 (HCOP-derived). `runtime_table.csv` contains the measurements shown under **Typical module runtimes** on the application's **Reproducibility** page; they are not performance guarantees for other computers.

`species_orthologs.rds` provides the separate, offline pig/mouse/rat mapping for new uploads. Its version, conservative eligibility rule, rebuild script and custom mapping/GMT formats are described in [annotation.md](annotation.md).
