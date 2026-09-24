# drug_discovery_pipeline

Automated in-silico screening pipeline for protein targets. Given a target
and known active compounds, it retrieves structurally similar molecules from
public databases, predicts complex structure and binding affinity with
Boltz-2, ranks candidates, and re-evaluates the top hits with GROMACS
molecular dynamics and MM-GBSA.

Developed for the MSc thesis "Targeting BACE1 Inhibitors: In Silico Drug
Discovery Approach" (Politecnico di Torino), supervised by Prof. Jacek A.
Tuszynski. Applied to BACE1, a target in Alzheimer's disease.

## Pipeline overview

1. **Retrieval:** First 29 compounds from DrugBank, then parsing from ChEMBL and PubChem
2. **Prediction:** Boltz-2 structure and affinity calculated for 3819 compounds
3. **Ranking:** top 27 by confidence, predicted affinity, structural diversity aggregated and ranked
4. **Validation:** GROMACS MD and MM-GBSA on the top candidates

## Results on BACE1
- 3819 compounds screened, 27 selected
- Correlation between Boltz-2 metrics and MM-GBSA: [r=-0.37 p=0.055 for aggregated score] [r=0.39 p=0.046 for Boltz-2 affinity score] Pearson correlation
- Ranked list: `results/top27.csv`; figures: `results/figures/`

## Repository structure
drug_discovery_pipeline/
    configs/
        environment.yml # Environment variables and python libraries required for the various steps 
    results/
        figures/
    pipeline/
        1.screening/
        2.simulation/
        3.analysis/

## Requirements
- Database parsing and analysis scripts run on a standard computer
  (Python 3, `configs/environment.yml`)
- Boltz-2 and GROMACS stages ran on Narval (Digital Research Alliance of Canada): 

    CPU: 1 to 8 core (1 core for Boltz-2, 8 for GROMACS)
    GPU: NVIDIA A100, 20 GB, MIG slice
    RAM: 16 GB

  Boltz-2 run on batches of max 500 compounds to reduce time consumption to an average of 3.32min/compound. Each batch was separated in the parsing stage and then selected using a general script which separated the jobs depending from the directory of input.

  GROMACS and MMPBSA simulation runs were run using the same resource allocation, with an average run of 6 and an half hour for each of the 27 compounds run. 

## Usage
Script contained into 'pipeline/' should be run in numerical order. 

## Limitations
- Results are computational only; no experimental validation
- MM-GBSA gives approximate rankings, not absolute binding free energies
- Sample size for the correlation analysis is small (n = 27)
- Full runs need GPU cluster resources
- Weak-to-moderate borderline correlation. 

## TO-DO

- Organize the scripts contained in 'pipeline/'
- Streamline the pipeline into a reduced number of scripts
- Run more tests in order to further validate the pipeline
- Integrate other tools like DrugClip or BindFlow
