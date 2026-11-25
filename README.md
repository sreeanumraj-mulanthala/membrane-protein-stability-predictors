# Membrane protein stability predictors: data and code

Data and code for the preprint "Benchmarking &Delta;&Delta;G Predictors on Membrane Proteins: Data Provenance, Held-Out Evaluation and Lipid Composition" (bioRxiv, DOI to be added).

## Contents

| File | Content |
|---|---|
| `data/annotation_table.csv` | Information for each of the 379 records contained in the MutTMP Stab-pred files; it includes identifier, selection step, cohort membership (LEGACY, PRIMARY or PRIMARY_SV), the source of the value and the eligibility class. No measured values are included. |
| `data/predictions_corpus.csv` | Predictions from DDGun, ESM-2, RaSP and ThermoMPNN for the 271 records in the LEGACY cohort. |
| `data/predictions_bama.csv` | Predictions from DDGun, ESM-2, RaSP, ThermoMPNN, mCSM-membrane and MutTMB Stab-pred for the 88 variants of BamA. |
| `structures/AF-P0A940-F1-model_v6_21-810_C690S_C700S.pdb` | Model of the BamA protein used in structure-based predictions. |
| `analysis.py` | Script to compute various statistics mentioned in the preprint (Spearman rho, Pearson r, mean absolute error, bootstrap intervals, leave-one-protein-out, predictor and fold contrasts, BamA correlations) using Python 3.9 or higher.|
| `heldout_ladder.py` | Refits the MutTMB Stab-pred model holding out one unit at a time, which relates to Figure 2.  Requires numpy, pandas and scikit-learn. |

All values follow the MPTherm convention: negative ddG indicates destabilization. Scores obtained from RaSP and ThermoMPNN predictions are multiplied by -1. ESM-2 scores are unitless log-likelihood ratios (log P(mutant) - log P(wild type)) and are used in rank statistics (Spearman rho) only.

## Not included, and where to get it

- **Measured ddG values of the 379 records.** Download `alpha.csv` and `beta.csv` from the MutTMP Stab-pred website and then place them in one folder. `analysis.py` will join them on `Mutation_ID`.
- **Measured BamA values.** Tables S8 (PC:PE) and S9 (PC:PG) of George et al., Chem. Sci. 2026, 17, 3103-3121 (DOI 10.1039/d5sc07027a). For the BamA part using `analysis.py`, transcribe these values into a CSV with columns `variant`, `condition` (`PC:PE` or `PC:PG`) and `ddG_F`, as printed.
- **The MutTMB Stab-pred reproducibility package**, for `heldout_ladder.py`, from the MutTMP Stab-pred website.
- **Corpus structures.** The PDB entries and AlphaFold models are listed in `annotation_table.csv` (`Structure_ID`, `PDB_Code`) and were used as downloaded.

## Usage

```
python analysis.py --muttmp-dir PATH_TO_FOLDER --out results
python analysis.py --muttmp-dir PATH_TO_FOLDER --bama-measured bama.csv --out results
python heldout_ladder.py --pkg PATH_TO/MutTMB_Stab_pred --out ladder.csv
```

`analysis.py` takes a few minutes (10,000 bootstrap replicates, fixed seed). The values in the paper for `heldout_ladder.py` were obtained with Python 3.13.5, scikit-learn 1.8.0, numpy 2.3.5 and pandas 2.2.3.

## BamA structure

AlphaFold DB model AF-P0A940-F1, version 6 (Jumper et al., Nature 2021, Varadi et al., Nucleic Acids Res. 2022), licensed CC BY 4.0. The first 20 residues, which represent the signal peptide, were removed. Serine was used to replace each of the cysteines at positions 690 and 700, and no other modifications were made to the protein.

## Licence

- Code (`analysis.py`, `heldout_ladder.py`): MIT, see `LICENSE`.
- Data files (`data/`, `structures/`): CC BY-NC 4.0, see `LICENSE-DATA`.

## Citation

DOI to be added.
