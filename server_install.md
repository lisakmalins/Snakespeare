# Server Installation
If you are using a work or lab server, ask your sysadmin if git and conda are installed already. If so, skip to STEP 2.

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

To __install git__:
```bash
conda install git
```

----
Note: On some server setups, conda may crash if your home folder is not writable. These commands will fix this by telling conda to store the environment in the current folder.

Navigate into the repository and run these lines to store the environment within the current folder.
```bash
conda config --add envs_dirs ./.conda/envs
conda config --add pkgs_dirs ./.conda/pkgs
```
