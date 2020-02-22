# Snakespeare
Snakespeare is a simple, entertaining, and bioinformatics-free Snakemake workflow
designed for first-time workflow users and workflow developers.

Snakespeare has many of the same steps and features as typical bioinformatics Snakemake pipelines. This is an ideal "practice pipeline" to become familiar with Snakemake before running other workflows or writing a workflow yourself.

## Snakespeare Results
This workflow calculates and plots how much different characters speak in Shakespeare's tragedies _Romeo & Juliet_ and _Hamlet_.

### Interesting statistics
- **Hamlet** talks the most with over **1428 lines** of iambic pentameter.

- Hamlet's uncle **Claudius** talks the second most with over **500 lines** of iambic pentameter. It must run in the family.

- Besides the Chorus in _Romeo and Juliet_, the **Ghost of King Hamlet** is the most long-winded with an average monologue length of **6.3 lines**.

- **Friar Lawrence** is a close second with an average monologue length of **6.2 lines**.

- **Romeo talks slightly more** than Juliet; however, **Juliet's lines are wittier**.


## Usage

### Step 1: Clone the repository
In the terminal, navigate to where you want to download Snakespeare.

Copy and paste these commands to "clone" this repository and then "change directory" into the folder.
```
git clone https://github.com/lisakmalins/Snakespeare.git
cd Snakespeare
```

### Step 2: Install miniconda3
Snakespeare uses conda, a popular package manager used by many bioinformaticians. By creating a conda environment, you can easily obtain all the software to run a pipeline without needing to install it system-wide.

First, you will need to __install miniconda3__. If you are running Snakespeare on a server that already has miniconda installed, skip to building the conda environment.

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The installer will ask you some questions to complete installation. Review and accept the license, accept or change home location, and answer yes to placing it in your path.

Some final commands to finish configuring conda:

```
source $HOME/.bashrc
conda config --add envs_dirs ./.conda/envs
conda config --add pkgs_dirs ./.conda/pkgs
```

### Step 3: Build and activate the conda environment
When you __build the conda environment__, Conda obtains all the software listed in `environment.yaml`. You only need to do this step once.
```
conda env create -f environment.yaml
```

Finally, you will need to __activate the environment__. The environment is named "snakespeare," and the software will only be accessible while the environment is active.
```
source activate snakespeare
```

### Step 4: Run Snakespeare
Run the snakemake workflow like this:
```
snakemake
```

That's it! The plots will appear in the folder `data/plots`.
