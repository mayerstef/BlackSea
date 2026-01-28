```markdown
# Black Sea Fish Biodiversity — range mapping + diversity indices (R)

Reproducible R workflow to (1) build species distribution products for Black Sea fishes from multi-source data, (2) aggregate outputs to a common grid, and (3) compute and analyze taxonomic, functional, and phylogenetic diversity patterns by bioregion.

Developed as part of an IMBRSea M.Sc. thesis project (2023) in the Modelling for Aquatic Systems (MAST) group (University of Liège).

---

## What this repository does

- Builds depth-constrained species range polygons from occurrence records with bias-aware preprocessing and bathymetry constraints.
- Integrates complementary distribution sources when available (e.g., IUCN/FAO layers; fisheries-derived proxies).
- Aggregates final distributions to a **0.1° × 0.1°** grid to generate a presence/absence community matrix.
- Computes taxonomic, functional, and phylogenetic diversity indices and produces maps/figures summarizing spatial patterns across Black Sea bioregions.

---

## Repository structure

- `Scripts/` — main workflow scripts (Part 1–3) and supporting functions
- `Rasters/` — intermediate raster layers / grids
- `MAYER_MasterThesis.Rproj` — RStudio project

---

## How to run

1. Open `MAYER_MasterThesis.Rproj` in RStudio.
2. Run the workflow scripts in `Scripts/` in order:
   - `Master_part1` — distribution/range products (multi-source mapping workflow)
   - `Master_part2` — functional and phylogenetic indices
   - `Master_part3` — statistical analyses and figure generation

---

## Data & access

This workflow draws on open biodiversity and fisheries sources (e.g., GBIF/OBIS occurrence data, IUCN/FAO distribution layers when available, and fisheries datasets), plus environmental layers for spatial constraints/analysis.

API keys/credentials are not stored in this repository. If an API-backed source is used (e.g., IUCN Red List), set your key via environment variables (e.g., `~/.Renviron`).

---

## Citation

Mayer, S. (2023). *Linking patterns in phylogeny, traits, and space for Black Sea fish* (M.Sc. thesis, IMBRSea / University of Liège).

## Contact

Stefanie Mayer (https://stefanieseas.wordpress.com/)
```

