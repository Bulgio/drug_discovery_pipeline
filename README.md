# drug_discovery_pipeline

Automated in-silico screening pipeline for protein targets. Given a target and known active compounds, it retrieves structurally similar molecules from public databases, predicts complex structure and binding affinity with Boltz-2, ranks candidates, and re-evaluates the top hits with GROMACS molecular dynamics and MM-PBSA (gmx_MMPBSA).

Developed for the MSc thesis "Targeting BACE1 Inhibitors: In Silico Drug Discovery Approach" (Politecnico di Torino), supervised by Prof. Jacek A. Tuszynski. Applied to BACE1, a target in Alzheimer's disease.

## Pipeline overview

1. **Retrieval** (`pipeline/Screening/`): seed compounds from DrugBank (10 compounds; requires a DrugBank account), then similar compounds from ChEMBL and PubChem (PubChem 2D similarity search, threshold [90]%).
2. **Prediction** (`pipeline/Boltz_sim/`): Boltz-2 structure and affinity prediction for 3,819 compounds, in batches of up to 500.
3. **Ranking** (`pipeline/Boltz_sim/post-sim/`): top 27 selected by Boltz-2 confidence and affinity rank, aggregated into a combined score, with a structural-diversity criterion. 

Combined_score = Confidence_placing*0.2+Affinity_score*0.8 (Arbitrary coefficients decided to balance confidence and affinity)

4. **Re-evaluation** (`pipeline/GROMACS/`): 10 ns MD per candidate (single replica), then MM-PBSA binding energy.

Flowcharts: `flowchart/`.

## Results on BACE1

- 3,819 compounds screened with Boltz-2; 27 selected and simulated.
- Full table: [`results/top27.csv`]. Figures: `results/images/`.
- MM-PBSA ΔG: mean -22.7 kJ/mol across the 27 (range -33.2 to -14.5).
- Agreement between Boltz-2 ranking metrics and MM-PBSA ΔG (n = 27):

| Boltz-2 metric  | Pearson r (p) | Spearman ρ (p) |
| Combined score  | -0.37 (0.055) | -0.36 (0.068)  |
| Affinity rank   | 0.39 (0.046)  | 0.31 (0.114)   |
| Confidence rank | 0.30 (0.133)  | 0.24 (0.235)   |

Exploratory only: none of these survive correction for testing three metrics. Sign convention: higher combined score, or lower rank number, means predicted stronger binding; more negative ΔG means stronger binding.

## Repository structure

```
drug_discovery_pipeline/
├── pipeline/
│   ├── Screening/      database retrieval and similarity search
│   ├── Boltz_sim/      Boltz-2 input generation, batch prediction, post-processing
│   └── GROMACS/        system preparation, MD (mdp/, simulation/), analysis, plotting
├── flowchart/          flowcharts of the pipeline stages
├── results/
│   ├── top27.csv
│   └── images/
└── configs/
    └── environment.yml
```

## Requirements

- Retrieval, parsing and analysis scripts: standard computer, Python 3.12, `configs/environment.yml`.
- Boltz-2 and GROMACS/MM-PBSA stages ran on Narval (Digital Research Alliance of Canada), one job per batch or compound:

| Stage             | GPU                                    | CPUs | RAM   |
| Boltz-2           | A100 MIG slice (`a100_3g.20gb`, 20 GB) | 1    | 16 GB |
| GROMACS + MM-PBSA | A100 MIG slice (`a100_3g.20gb`, 20 GB) | 8    | 16 GB |

- Boltz-2: batches of up to 500 compounds, about 3.3 min per compound (~211 GPU-hours in total).
- GROMACS: about 6.5 h per compound (~176 GPU-hours for 27).
- Software: GROMACS 2025.4, AmberTools 25.0, gmx_MMPBSA, RDKit 2024.09.6, Boltz-2.
- MD setup: amber99sb-ildn (protein), GAFF2 (ligand, via antechamber), SPC/E water, 2 fs time step, 10 ns production, v-rescale thermostat, Parrinello-Rahman barostat.

## Usage

Scripts in `pipeline/` are meant to be run in numerical stage order (see `flowchart/`). 
Paths and Slurm accounts are cluster-specific and must be edited depending on the system.

## Limitations

- Computational only; no experimental validation.
- MM-PBSA without an entropy term gives approximate rankings, not absolute binding free energies.
- Single 10 ns replica per compound; ligand RMSD is large for most candidates (mean 5.1 Å; only 1 of 27 below 3 Å), so poses are not demonstrably stable.
- Small sample (n = 27) drawn from the top of the Boltz-2 ranking; the restricted range weakens any correlation.
- No decoys or negative controls, and possible overlap between Boltz-2 training data and the compounds screened.
- Not tested on consumer hardware; the Boltz-2 and MD stages need a GPU.
- Scripts contain cluster-specific paths and are not yet packaged for other systems.

## TO-DO

- Organize the scripts contained in `pipeline/`.
- Streamline the pipeline into a reduced number of scripts.
- Score the known DrugBank actives with the same protocol to test enrichment.
- Integrate other tools such as DrugCLIP or BindFlow.
