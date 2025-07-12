# Kernelized Discriminant Analysis for Joint Modeling of Multivariate Categorical Response

This repository provides the simulation code for the paper titled *Kernelized Discriminant Analysis for Joint Modeling of Multivariate Categorical Response*.


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

ml R/4.2
R CMD BATCH --vanilla Main-A4.R Results-A4/Rep_${SLURM_ARRAY_TASK_ID}.Rout
```
3. Create a directory to store the results.

```
mkdir Results-A4
```

4. Execute the modified `Sim.sh` file.

This process will launch 1200 jobs, corresponding to 100 replicates for each combination of:
+ $p$ values: 50, 100, 150
+ $\nu$ values: 0.8, 1.0, 1.2, 1.4

For every replicate, the script will execute `Main-A4.R` to generate data based on the specified data-generating model, fit various models, and assess their performance. All the results are saved in RDS files under the `Results-A4` directory.

All necessary functions can be found in the Functions directory.


## Notes on Computing Environment

The provided workflow is configured for SLURM-based high-performance computing (HPC) environments. The `Sim.sh` script is a SLURM batch job script that manages parallel execution using job arrays.

If you wish to run the simulation on a different computing platform (e.g., a non-SLURM cluster, local multicore machine, or cloud environment), **you will need to modify**:

- The shell script `Sim.sh`: Replace SLURM-specific directives (e.g., `#SBATCH` flags) with those suitable for your job scheduling system (e.g., PBS, LSF, or a custom script for GNU parallel or background jobs).
- The R script `Main-XX.R`: Ensure that any job-specific environment variables (e.g., `SLURM_ARRAY_TASK_ID`) are replaced or emulated appropriately (e.g., using `commandArgs()` or manually setting replicate IDs).

Additionally, make sure all dependencies such as R packages and source files from the `Functions` directory are properly loaded in your environment.
