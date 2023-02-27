# GATK - Best Practices

This workflow is designed to perform the preprocessing steps for GATK

## Setup
### Code
From GitLab, fork this project filling out the project name, description, and ensure under namespace "SGT-Projects" is selected.  
Then, clone this repository into your project directory
```
git clone <url>
```

Checkout new analysis branch and push this to your project repository
```
git checkout -b analysis
git push -u origin analysis
```

Change your default branch to `analysis` in GitLab via setting -> repository -> default branch.

### Conda
Install conda to `conda/` directory
```
wget wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p conda
```

Activate the base environment
```
source conda/bin/activate base
```

Install mamba package manager
```
conda install mamba conda -y
```

Install Snakemake and activate environment
```
mamba env create -f snakemake.yaml
conda activate snakemake
```

To run the example data using SLURM, from the `.test/` directory, run
```
./snakemake_slurm.sh -d .test --use-conda -j500
```

## Running
Setup the correct genome and ensure the settings in `config/config.yaml` are correct, then prepare the workflow configuration files `fastqs.tsv`, `samples.tsv`, `analysis.yaml`.    
See documentation [here](config/README.md) for more information.  

It is often useful to see what Snakemake is planning on doing without running anything
```
snakemake -n
```

To run the workflow on SLURM one could use
```
./snakemake_slurm.sh --use-conda -j<threads>
```

This will instruct Snakemake to resolve dependencies using conda and up to `<threads>` number of parallel tasks (or number of threads per-task).  
When running locally use `snakemake`.  
See Snakemake documentation for more [command line options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options).


Once an analysis is complete, ensure all your changes are commited and pushed to your project's repository.
