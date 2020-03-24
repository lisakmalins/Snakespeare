# Snakespeare
Snakespeare is a simple, entertaining, and bioinformatics-free Snakemake workflow
designed for first-time workflow users and workflow developers.

The moving parts of Snakespeare are similar to many bioinformatics Snakemake pipelines:
- `Snakefile` – rules for all steps of workflow are written out here
- `config.yaml` – parameters users can customize are listed here
- `environment.yaml` – lists dependencies to be installed into conda virtual environment
- `data/` – all input and output files live in this directory
- `scripts/` – all Python and R scripts live in this directory

This is an ideal practice pipeline to become familiar with Snakemake before running other workflows or writing a workflow yourself.

## Snakespeare results
This workflow calculates and plots how much different characters speak in Shakespeare's tragedies _Romeo & Juliet_ and _Hamlet_.

### Interesting statistics
- **Hamlet** talks the most with over **1428 lines** of iambic pentameter.

- Hamlet's uncle **Claudius** talks the second most with over **500 lines** of iambic pentameter. It must run in the family.

- Besides the Chorus in _Romeo and Juliet_, the **Ghost of King Hamlet** is the most long-winded with an average monologue length of **6.3 lines**.

- **Friar Lawrence** is a close second with an average monologue length of **6.2 lines**.

- **Romeo talks slightly more** than Juliet; however, **Juliet's lines are wittier**.


## Usage
You can run Snakespeare either on your computer or on a server.

### STEP 1: Install miniconda and git

#### Installing conda and git on your computer
If you want to run Snakespeare on your computer, you can __download git and miniconda__ from the following links:
- Download git: https://git-scm.com/downloads
- Download miniconda: https://conda.io/en/latest/miniconda.html

After you have git and conda installed, continue to STEP 2.

#### Installing conda and git on a server
If you are using a work or lab server, ask your sysadmin if git and conda are installed already. If so, skip to STEP 2.


To __download miniconda__:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The installer will ask you some questions to complete installation. Review and accept the license, accept or change home location, and answer yes to placing it in your path.

To finish configuring miniconda:
```
source $HOME/.bashrc
```

To __install git__:
```
conda install git
```

### STEP 2: Clone the repository

In the terminal, navigate to where you want to download Snakespeare.

Copy and paste these commands to __clone this repository__ and then "change directory" into the folder.
```
git clone https://github.com/lisakmalins/Snakespeare.git
cd Snakespeare
```

### STEP 3: Build and activate the conda environment
When you __build the conda environment__, Conda obtains all the software listed in `environment.yaml`. You only need to do this step once.
```
# Recommended: prevent conda from crashing if home folder is not writable
conda config --add envs_dirs ./.conda/envs
conda config --add pkgs_dirs ./.conda/pkgs

# Build snakespeare conda environment
conda env create -f environment.yaml
```

Finally, you will need to __activate the environment__. The environment is named "snakespeare," and the software will only be accessible while the environment is active.
```
source activate snakespeare
```

### STEP 4: Run Snakespeare
__Run the snakemake workflow__ like this:
```
snakemake
```

That's it! The workflow should finish within a few seconds. The plots will appear in the folder `Snakespeare/data/plots`.
