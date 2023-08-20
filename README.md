# Kernelized Discriminant Analysis
This repository contains the simulation code for the paper titled "Kernelized Discriminant Analysis for Multivariate Categorical Response Regression."

## Running the Simulation
To execute the simulation:

1. Open the `Sim.sh` file.

2. Replace all instances of `XX` in the Sim.sh file with the desired model code: `A4`, `A6`, `B4`, or `B6`.

For instance, to run the simulation for Model A-4:
```
#!/bin/bash
#SBATCH --job-name=KernelDA
#SBATCH -o Results-A4/Rep_%a.Rout
#SBATCH --array=1-1200
#SBATCH --mail-user=xxx@xyz
#SBATCH --mail-type=END
#SBATCH --account=abc
#SBATCH --qos=abc
#SBATCH --mem=2gb
#SBATCH -t 96:00:00

mkdir Results-A4
ml R/4.2
R CMD BATCH --vanilla Main-A4.R Results-A4/Rep_${SLURM_ARRAY_TASK_ID}.Rout
```

3. Execute the modified `Sim.sh` file.

This process will launch 1200 jobs, corresponding to 100 replicates for each combination of:
+ $p$ values: 50, 100, 150
+ $\nu$ values: 0.8, 1.0, 1.2, 1.4

For every replicate, the script will:

+ Create a `Results-A4` directory to store the results.
+ Execute `Main-A4.R` to generate data based on the specified data-generating model, fit various models, and assess their performance. All the results are saved in RDS files under the `Results-A4` directory.

All necessary functions can be found in the Functions directory.

