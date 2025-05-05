## AutoDock Vina Installation Steps

This guide provides a streamlined procedure to install and set up AutoDock Vina.

**1. Install Anaconda (if not already installed):**

   ```bash
   wget [https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh](https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh)
   bash Anaconda3-2024.02-1-Linux-x86_64.sh
   source ~/.bashrc
```



This downloads and installs Anaconda, a Python distribution that simplifies package management. Follow the prompts during the installation.


**2. Create a Conda Environment:**
```bash
conda create -n py3112 python=3.11.2
conda activate py3112
```

This creates an isolated environment named "py3112" with Python 3.11.2. Using environments helps avoid conflicts with other software. conda activate py3112 activates the newly created environment

**3. Install AutoDock Vina:**

```bash
# conda install -c conda-forge autodock-vina
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install vina
```

**4. Install Open Babel:**

```bash
sudo apt-get install openbabel
obabel -V
```
