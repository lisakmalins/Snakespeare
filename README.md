# Snakespeare

<br /><br />
<div align="center">
  <img src="images/snakespeare_logo.png" alt="Snakespeare logo" width="500px" />
</div>
<br /><br />

Snakespeare is a simple, entertaining, and bioinformatics-free Snakemake workflow designed for first-time workflow users and workflow developers.

The moving parts of Snakespeare are similar to many bioinformatics Snakemake pipelines:
- [`Snakefile`](Snakefile) – contains rules for all steps of workflow
- [`config.yaml`](config.yaml) – parameters users can customize are listed here
- [`environment.yaml`](environment.yaml) – lists software dependencies to be installed into conda virtual environment
- [`scripts/`](scripts) – all Python and R scripts live in this directory
- [`data/`](data) – all input and output files live in this directory


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

### STEP 1: Install miniconda and git
To run Snakespeare, you will need __git__ and __conda__.

<details>
<summary><b>Click here for instructions for Windows</b></summary>

<details>
<summary>Run Snakespeare via Anaconda prompt (recommended for new users)</summary>

#### Installing Miniconda3 for Windows
Head over to the Anaconda website and download a [Windows installer for Miniconda3](https://docs.conda.io/en/latest/miniconda.html#windows-installers).
> If you are not sure which to choose, pick the highest version of Python.
>
> You can check whether your system is 64-bit or 32-bit under __Settings__ > __About__ > __Device specifications__ > __System type__.

Run the installer and follow the instructions to complete your installation of Miniconda3.

#### Open Anaconda Prompt
Now click the Start menu and search for "__Anaconda prompt__." This is a version of the Windows "command prompt" terminal that includes miniconda.

#### Installing Git for Windows
In Anaconda prompt, copy and paste the following to install git:
```sh
conda install git
```

That's it! Continue to STEP 2.
</details>

<details>
<summary>Run Snakespeare via WSL (advanced users)</summary>

If you are already using Windows Subsystem for Linux (aka the Ubuntu app), see [command-line instructions](command_line_install.md) for how to install miniconda and git in your Linux terminal. Then continue to STEP 2.
</details>
</details>

<details>
<summary><b>Click here for instructions for Mac</b></summary>

#### Installing Git for Mac
On your Mac, open Terminal. Type `git` and press Enter.
- If a bunch of text appears (these are the usage instructions for git), congratulations, you already have git installed! Skip to **Installing Miniconda for Mac**.
- If you see `git: command not found`, then you will need to get git for Mac. The easiest method is to [install Xcode](https://apps.apple.com/us/app/xcode/id497799835), which is a suite of developer tools provided by Apple.
- After installing Xcode, open a *new* terminal window and try typing `git` again. You should see the usage instructions now.
> If you still see `git: command not found`, please [let me know](https://github.com/lisakmalins/Snakespeare/issues/new) so I can help.

#### Installing Miniconda for Mac
- To get Miniconda for Mac, download an installer from the [Anaconda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html).
- If you are not sure which to choose, download the __Python 3.9 Miniconda3 MacOSX 64-bit pkg__.
- Run the installer that just downloaded, and follow the instructions to complete your installation of Miniconda.

Done! Make sure you close any terminal windows that you have open, then continue to STEP 2.

</details>

<details>
<summary><b>Click here for instructions for Linux</b></summary>

<details>
<summary>Linux desktop users</summary>

- Head to the [git website](https://git-scm.com/download/linux) for instructions to install git with your distribution's package manager.
- Head to the [Anaconda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) for instructions to download and run a Miniconda installer.

After installing git and miniconda, close any terminal windows you have open and continue to STEP 2.
</details>


<details>
<summary>Linux server users</summary>

If you would like to run Snakespeare on a server, check with your supervisor or sysadmin to see if the server already has git and conda installed. If you do need to install software (and have permission to do so), please see the [command-line instructions](command_line_install.md). Then continue to STEP 2.
</details>
</details>

### STEP 2: Clone the repository
Open a new terminal window and navigate to where you want to download Snakespeare.

If you are not sure, I recommend you copy and paste these commands to __make a new directory called `GitHub_repos`__, then "change directory" into the folder:
```bash
mkdir GitHub_repos
cd GitHub_repos
```

Copy and paste these commands to __clone this repository__ and then "change directory" into the folder.
```bash
git clone https://github.com/lisakmalins/Snakespeare.git
cd Snakespeare
```

### STEP 3: Build and activate the conda environment
When you __build the conda environment__, Conda obtains all the software listed in `environment.yaml`. You only need to do this step once.

```bash
conda env create -f environment.yaml
```

Finally, you will need to __activate the environment__. The environment is named "snakespeare," and the software will only be accessible while the environment is active.
```bash
conda activate snakespeare
```
> Note: for older versions of Anaconda, you may need to use the command `source activate snakespeare` instead.

When you want to deactivate the environment later, you can do so with the command `conda deactivate`.

### STEP 4: Run Snakespeare
__Run the snakemake workflow__ like this:
```bash
snakemake
```

That's it! The workflow should finish within a few seconds. The plots will appear in the folder `Snakespeare/data/plots/`.
