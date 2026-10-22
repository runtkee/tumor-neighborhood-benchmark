# Third-party data

The MIT license covers the code. The data below keeps its own license.

| What | Where in this repo | Source | License |
|---|---|---|---|
| Per-patient summary statistics derived from IMMUcan IMC panel 1 (179 patients) | `src/bio/kit/reference/IMMUcan_IMC_panel1.csv`, embedded in `docs/index.html` | IMMUcan consortium, Zenodo [10.5281/zenodo.12912567](https://zenodo.org/records/12912567) | CC BY 4.0 |
| Per-patient summary statistics derived from CODEX colorectal cancer (35 patients); 3-patient test subset of the cell table | `src/bio/kit/reference/CODEX_CRC_Schurch2020.csv`, embedded in `docs/index.html`; `tests/fixtures/codex_subset_3pt.csv` | Schürch et al., *Cell* 2020; Mendeley Data [10.17632/mpjzbtfgfr.1](https://data.mendeley.com/datasets/mpjzbtfgfr/1) | CC BY 4.0 |
| Per-patient summary statistics derived from the SPATCH platform benchmark (3 patients, 4 runs) | `src/bio/kit/reference/SPATCH_*.csv` | SPATCH, *Nat Commun* 2025, doi:10.1038/s41467-025-64292-3; spatch.pku-genomics.org | Data license not stated; only derived summary statistics are included |
| Per-patient summary statistics derived from the colorectal cancer IMC table of Marteau et al. (37 patient IDs, intratumoral and margin regions) | `src/bio/kit/reference/MARTEAU_CRC_IMC_*.csv`, embedded in `docs/index.html` | Marteau et al., *Cancer Cell* 2025 (doi:10.1016/j.ccell.2025.12.003); BioImage Archive [S-BIAD2208](https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD2208) | CC BY 4.0 |
| Per-patient summary statistics derived from the colorectal cancer Xenium sections of Marteau et al. (15 patients; tumor core, invasive margin and adjacent normal sections; cell types called by this kit from the 10x clusters) | `src/bio/kit/reference/MARTEAU_XENIUM_*.csv`, embedded in `docs/index.html` | Marteau et al., *Cancer Cell* 2025; BioImage Archive [S-BIAD2208](https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD2208) | CC BY 4.0 |
| Per-patient summary statistics derived from Xenium 5K sections of the six-cancer platform comparison (colorectal, lung, breast, ovarian; one section each) | `src/bio/kit/reference/CERVILLA_*.csv`, embedded in `docs/index.html` | Cervilla et al. 2025; Zenodo [10.5281/zenodo.18000256](https://zenodo.org/records/18000256) | CC BY 4.0 |
| Tiny HALO and inForm example exports | `tests/fixtures/tiny_halo.csv`, `tests/fixtures/tiny_inform.txt` | SPIAT R package, [TrigosTeam/SPIAT](https://github.com/TrigosTeam/SPIAT) `inst/extdata` | Artistic-2.0 |

Changes: all public datasets were re-processed with this project's locked recipe (v3p.1). Cell-type labels were mapped to roles as listed in `kit_manifests/`. CODEX subsets keep only the columns the recipe uses.
