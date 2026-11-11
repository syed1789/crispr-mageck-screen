# CRISPR Screen Hit-Calling with MAGeCK

**Author:** Syed Nurul Hasan, Ph.D.
**Contact:** nayanchoton1789@gmail.com | [ORCID: 0000-0002-4564-872X](https://orcid.org/0000-0002-4564-872X)

---

## Why this project

Several computational-biology roles I've applied to (functional genomics, oncology target discovery) call for hands-on CRISPR screen analysis with MAGeCK — sgRNA/gene-level hit-calling, QC, and functional enrichment on top hits. That work is real (see Professional Experience in my CV), but it lived entirely inside institutional projects with no public artifact behind it. This repo closes that gap with a small, fully reproducible, end-to-end MAGeCK run on a real published screen, runnable from a clean checkout in minutes.

**Scope note, stated plainly:** this is pooled-screen, sequencing-count-based hit-calling — not optical/imaging-based CRISPR screen analysis, which is separate territory I haven't worked in.

---

## Pipeline

```
sgRNA count table (real published screen, HL60 + KBM7 leukemia lines)
        │
        ▼
1. MAGeCK RRA test (Docker: biocontainers/mageck)
   ├── Compare T-final vs. T-initial per cell line
   ├── sgRNA-level and gene-level robust rank aggregation (RRA)
   └── Output: gene_summary.txt, sgrna_summary.txt
        │
        ▼
2. Select top negatively-selected ("essential-gene candidate") hits
   by neg|fdr
        │
        ▼
3. Functional enrichment (Docker: clusterProfiler + org.Hs.eg.db)
   ├── GO Biological Process (enrichGO)
   └── KEGG pathways (enrichKEGG)
        │
        ▼
   go_bp_enrichment.tsv · kegg_enrichment.tsv · go_bp_dotplot.png
```

## Dataset

`data/sgrna_counts_hl60_kbm7.txt` — sgRNA read-count table from the genome-scale
CRISPR-Cas9 knockout screen in HL60 and KBM7 human leukemia cell lines
(Wang et al., *Science* 2014, "Genetic Screens in Human Cells Using the
CRISPR-Cas9 System"), fetched via `scripts/download_data.sh` from
[davidliwei/mageck2-demo](https://github.com/davidliwei/mageck2-demo) — the
canonical MAGeCK RRA tutorial dataset, distributed under BSD-3-Clause.
Four columns: `HL60.initial`, `KBM7.initial` (T0) vs. `HL60.final`,
`KBM7.final` (screen endpoint), 999 sgRNAs across 100 genes.

## Results (this run)

Top negatively-selected genes (lowest `neg|fdr` — sgRNAs depleted between T0
and endpoint in both cell lines, i.e. essential-gene candidates in this
screen) include **ACLY** (ATP citrate lyase — core lipogenesis enzyme) and
**AHCY** (S-adenosylhomocysteine hydrolase — essential one-carbon
metabolism), both independently known essential genes in DepMap
pan-cancer screens — a reasonable sanity check that the RRA hit-calling on
this small demo table is behaving as expected. Full ranked tables:
[`results/hl60_kbm7.gene_summary.txt`](results/hl60_kbm7.gene_summary.txt),
[`results/hl60_kbm7.sgrna_summary.txt`](results/hl60_kbm7.sgrna_summary.txt).

GO:BP and KEGG enrichment on the top 30 negatively-selected genes independently
recovers TCA-cycle/citrate-cycle metabolism (`hsa00020`) and glycerolipid/
glycerophospholipid metabolism (`hsa00561`, `hsa00564`) — directly consistent
with ACLY and AHCY surfacing as top hits, since both sit in those pathways.
Full tables: [`enrichment/go_bp_enrichment.tsv`](enrichment/go_bp_enrichment.tsv),
[`enrichment/kegg_enrichment.tsv`](enrichment/kegg_enrichment.tsv),
[`enrichment/go_bp_dotplot.png`](enrichment/go_bp_dotplot.png).

## Running it

```bash
# 1. Fetch the count table
bash scripts/download_data.sh

# 2. MAGeCK RRA hit-calling (Docker)
docker run --rm -v "$(pwd)":/work -w /work/results \
  quay.io/biocontainers/mageck:0.5.9.5--py312hf731ba3_8 \
  mageck test -k /work/data/sgrna_counts_hl60_kbm7.txt \
    -t HL60.final,KBM7.final -c HL60.initial,KBM7.initial -n hl60_kbm7

# 3. Functional enrichment on top hits (Docker)
docker build -t mageck-enrichment:latest enrichment/
docker run --rm -v "$(pwd)":/work mageck-enrichment:latest \
  Rscript /work/enrichment/run_enrichment.R
```

## Key Features

- Real published screen data, not synthetic counts
- MAGeCK RRA `test` mode — sgRNA- and gene-level hit-calling, matching the
  standard two-condition pooled-screen workflow
- GO:BP + KEGG functional enrichment (clusterProfiler) on top hits
- Fully containerized — no local MAGeCK, R, or Bioconductor install required

## References

1. Wang T, Wei JJ, Sabatini DM, Lander ES. Genetic Screens in Human Cells
   Using the CRISPR-Cas9 System. *Science*. 2014;343(6166):80-84.
2. Li W, Xu H, Xiao T, et al. MAGeCK enables robust identification of
   essential genes from genome-scale CRISPR/Cas9 knockout screens.
   *Genome Biology*. 2014;15:554.
3. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing
   biological themes among gene clusters. *OMICS*. 2012;16(5):284-287.

## License

MIT — see [LICENSE](LICENSE).

## Running on HPC (SLURM) and AWS

Includes `run_analysis.sh` (Docker or Apptainer), a SLURM job script, and a one-shot EC2 launcher. See [HPC_AWS.md](HPC_AWS.md) for steps and what was verified.
