# Snakespeare

<br /><br />
<div align="center">
  <img src="images/snakespeare_logo.png" alt="Snakespeare logo" width="500px" />
</div>
<br /><br />

Snakespeare is a simple, entertaining, and bioinformatics-free workflow designed for first-time Snakemake users and developers.

## Snakespeare learning goals
This tutorial is ideal for beginners, including folks who have never used the terminal before.

After working through this tutorial, you will learn:
- How to __access the terminal__ on your computer (Windows, Mac, or Linux)
- How to __clone a repository__ from GitHub
- How to __build and activate a conda environment__
- How to __run a Snakemake workflow__ and view the results

### Diving deeper
In addition, this repository is a simple demonstration of Snakemake for workflow developers.

The moving parts of Snakespeare are identical to the Snakemake pipelines I use every day for bioinformatics work:
- [`Snakefile`](Snakefile) – contains rules for all steps of workflow
- [`config.yaml`](config.yaml) – parameters users can customize are listed here
- [`environment.yaml`](environment.yaml) – lists software dependencies to be installed into conda virtual environment
- [`scripts/`](scripts) – all Python and R scripts live in this directory
- [`data/`](data) – all input and output files live in this directory

## Snakespeare results
This workflow calculates and plots how much characters speak in Shakespeare's tragedies _Hamlet_ and _Romeo & Juliet_.

![Example output from Snakespeare workflow.](images/all_statistics.png)

### Interesting statistics
- **Hamlet** talks the most with over **1428 lines** of iambic pentameter.

- Hamlet's uncle **Claudius** talks the second most with over **500 lines** of iambic pentameter. It must run in the family.

- Besides the Chorus in _Romeo and Juliet_, the **Ghost of King Hamlet** is the most long-winded with an average speech length of **6.3 lines**.

- **Friar Lawrence** is a close second with an average speech length of **6.2 lines**.

- **Romeo talks slightly more than Juliet** (however, Juliet's lines are wittier).


## Usage

### STEP 1: Install miniconda and git
To run Snakespeare, you will need two pieces of software: __git__ and __conda__.
- __git__ is a tool for downloading the code for Snakespeare from GitHub.
- __conda__ is a tool for accessing all software dependencies (including R, Python, and Snakemake).
> All software dependencies will be installed into a "virtual environment," so Snakespeare will not conflict with any Python or R software you have set up already.

<!------------------------------ Begin Windows instructions ------------------------------>
<details>
<summary><b>Click here for instructions for Windows</b></summary>
<table>

<!-- - - - - - - - - - - - - - - Windows + Anaconda Prompt - - - - - - - - - - - - - - -->
<tr><td><details>
<summary>Run Snakespeare via Anaconda prompt (easiest for beginning users)</summary>

#### Installing Miniconda3 + Anaconda Prompt for Windows
Head over to the Anaconda website and download a [Windows installer for Miniconda3](https://docs.conda.io/en/latest/miniconda.html#windows-installers).
> If you are not sure which to choose, pick the highest version of Python.
>
> You can check whether your system is 64-bit or 32-bit under __Settings__ > __System__ > __About__ > __Device specifications__ > __System type__.

Run the installer and follow the instructions to complete the installation. This software bundle includes Miniconda3 as well as Anaconda Prompt, which is a terminal app that you can use to run Snakespeare.

#### Open Anaconda Prompt
Now click the Start menu and search for "__Anaconda prompt__." This is a modified version of Windows Command Prompt (`cmd.exe`) that is pre-loaded with the conda executable.

#### Installing Git in Anaconda Prompt
In Anaconda prompt, copy and paste the following to install git:
```sh
conda install -y git
```

That's it! Continue to STEP 2.
</details></td></tr>

<tr><td><details>
<summary>Run Snakespeare via Git Bash (good for beginning users)</summary>

#### Installing Miniconda3 for Windows
Head over to the Anaconda website and download a [Windows installer for Miniconda3](https://docs.conda.io/en/latest/miniconda.html#windows-installers).
> If you are not sure which to choose, pick the highest version of Python.
>
> You can check whether your system is 64-bit or 32-bit under __Settings__ > __System__ > __About__ > __Device specifications__ > __System type__.

Run the installer and follow the instructions to complete your installation of Miniconda3.

#### Installing Git + Git Bash for Windows
Head to the [git website](https://git-scm.com/download/win) and download an installer for Windows. Run the installer and follow the instructions to complete the installation. This software bundle includes git as well as Git Bash, which is a terminal app that you can use to run Snakespeare.

__Important:__ While installing Git for Windows, be sure to check the box to add "Git Bash here" to the File Explorer context menu. You'll need it for the next step.

#### Enabling Conda in Git Bash
To enable Conda within Git Bash, you'll need to add the Conda startup script to your `~/.bashrc` file, which executes every time you open Git Bash.

From the Start menu, search for "Miniconda3" and click "Open File Location." Within that folder, navigate to `etc` and then `profile.d`. You should see a file called `conda.sh` in this folder. Right-click inside the window and select "Git Bash here" to open a terminal window in this folder.

Run the following command to to add the Conda startup script to your `~/.bashrc`:
```sh
echo ". '${PWD}'/conda.sh" >> ~/.bashrc
```

After that, close the terminal window.

Finally, let's double-check that conda is working in your new Git Bash terminal. From the Start menu, open __Git Bash__. Type `conda` and press Enter. If a bunch of text appears (these are the usage instructions for conda), congratulations, you're all set up! Continue to STEP 2.


</details></td></tr>

<!-- - - - - - - - - - - - - - - Windows + WSL - - - - - - - - - - - - - - -->
<tr><td><details>
<summary>Run Snakespeare via Windows Subsystem for Linux (advanced users)</summary>
If you are already using Windows Subsystem for Linux, follow the instructions below for how to install miniconda and git in your Ubuntu terminal.

#### Installing Git in WSL
Head to the [git website](https://git-scm.com/download/linux) and follow the installation instructions for Ubuntu.

#### Installing Miniconda in WSL
Head to the [Anaconda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) for instructions to download and run a Miniconda installer for Linux.

After installing git and miniconda, close any terminal windows you have open and continue to STEP 2.
</details></td></tr>

</table>
</details>

<!------------------------------ End Windows instructions ------------------------------>
<!------------------------------ Begin Mac instructions ------------------------------>
<details>
<summary><b>Click here for instructions for Mac</b></summary>
<table>

<tr><td>

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

</td></tr>
</table>
</details>

<!------------------------------ End Mac instructions ------------------------------>
<!------------------------------ Begin Linux instructions ------------------------------>
<details>
<summary><b>Click here for instructions for Linux</b></summary>
<table>

<!-- - - - - - - - - - - - - - - Linux Desktop - - - - - - - - - - - - - - -->
<tr><td><details>
<summary>Linux desktop users</summary>

#### Installing Git for Linux
Head to the [git website](https://git-scm.com/download/linux) for instructions to install git with your distribution's package manager.

#### Installing Miniconda for Linux
Head to the [Anaconda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) for instructions to download and run a Miniconda installer.

After installing git and miniconda, close any terminal windows you have open and continue to STEP 2.
</details></td></tr>

<!-- - - - - - - - - - - - - - - Linux Server - - - - - - - - - - - - - - -->
<tr><td><details>
<summary>Linux server users</summary>

If you would like to run Snakespeare on a work or lab server, check with your supervisor or sysadmin to see if git and conda are installed already. If so, continue to STEP 2.

Otherwise, if you need to install software (and have permission to do so), follow the instructions below.

#### Installing Miniconda on a Linux Server
To __install miniconda__ from the command line:
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The installer will ask you some questions to complete installation. Review and accept the license, accept or change home location, and answer yes to placing it in your path.

To finish configuring miniconda:
```bash
source $HOME/.bashrc
```

> Note: If your home folder is not writable on your server, conda will crash. If you experience this issue, run these commands to tell conda to store the environment in the current folder.
> ```bash
> conda config --add envs_dirs ./.conda/envs
> conda config --add pkgs_dirs ./.conda/pkgs
> ```

#### Installing Git on a Linux Server
To __install git__:
```bash
conda install git
```

</details></td></tr>
</table>
</details>
<!------------------------------ End Linux instructions ------------------------------>

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
