# rsomics-unifrac

Pairwise **UniFrac** phylogenetic beta-diversity between microbiome samples —
unweighted, weighted, and weighted-normalized — from a feature count table and a
rooted Newick tree. Emits a square scikit-bio-compatible `DistanceMatrix` TSV.

```
rsomics-unifrac table.tsv --tree tree.nwk                         # unweighted
rsomics-unifrac table.tsv --tree tree.nwk --weighted              # weighted
rsomics-unifrac table.tsv --tree tree.nwk --weighted --normalized -o dm.tsv
```

The count table is feature-by-sample: first column is the OTU/taxon ID (matching
a tip name in the tree), the header row names the samples, each cell is a count.
Reads stdin when the table is `-` or omitted; pass `--csv` for comma-separated.

## Method

Counts are pushed to the tips and summed up the tree (postorder reduction), so
every branch carries the total count of the leaves below it. For a pair of
samples `u`, `v` with tip totals `A_T`, `B_T` and per-branch descendant counts
`A_i`, `B_i`:

- **unweighted** = (Σ branch lengths spanned by exactly one sample) / (Σ branch
  lengths spanned by either) — branch `i` counts when `A_i > 0` differs from
  `B_i > 0`.
- **weighted** = Σ_i `bl_i · |A_i/A_T − B_i/B_T|`.
- **weighted-normalized** divides the weighted sum by
  Σ_tips `d_tip · (A_tip/A_T + B_tip/B_T)`, where `d_tip` is the tip-to-root
  distance, putting the distance in `[0, 1]`.

The unordered upper triangle is evaluated in parallel over sample pairs.

## Origin

This crate is an independent Rust implementation based on:
- Lozupone & Knight 2005, unweighted UniFrac (DOI 10.1128/AEM.71.12.8228-8235.2005)
- Lozupone, Hamady, Kelley & Knight 2007, weighted UniFrac (DOI 10.1128/AEM.01996-06)
- Hamady, Lozupone & Knight 2010, the array-based formulation (DOI 10.1038/ismej.2009.97)

It reproduces `skbio.diversity.beta_diversity("unweighted_unifrac" |
"weighted_unifrac" | "weighted_normalized_unifrac", …)`. scikit-bio is
BSD-3-Clause; its `skbio.diversity.beta._unifrac` and the `_phylogenetic`
postorder/tip-distance routines were read to match exact semantics (the NaN root
length is treated as 0; an all-zero pair returns 0). Compat tests diff against
scikit-bio value-exact (~1e-9).

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio https://scikit-bio.org (BSD-3-Clause).
