
---

# 🧬 Molecular Docking Environment Setup (Linux, Python 3.11.2)

This guide provides step-by-step instructions to set up a molecular docking environment using Python 3.11.2. The setup includes:

* **Anaconda**: Package and environment management.
* **AutoDock Vina**: Docking engine.
* **AutoDockTools (Python 3 version)**: Tools for preparing `.pdbqt` files.
* **Open Babel**: Format conversion tool.
---

## 📦 1. Install Anaconda

Anaconda simplifies package management and deployment.

1. **Download the Anaconda installer:**

   ```bash
   wget https://repo.anaconda.com/archive/Anaconda3-2023.07-1-Linux-x86_64.sh
   ```

2. **Run the installer:**

   ```bash
   bash Anaconda3-2023.07-1-Linux-x86_64.sh
   ```

   Follow the prompts to complete the installation.

3. **Activate Anaconda:**

   ```bash
   source ~/.bashrc
   ```

---

## 🧪 2. Create and Activate a Conda Environment

Creating a separate environment helps manage dependencies effectively.

```bash
conda create -n vina-env python=3.11.2 -y
conda activate vina-env
```

---

## 🔧 3. Install AutoDock Vina

AutoDock Vina is a popular docking engine for molecular docking simulations.

1. **Add the `conda-forge` channel and set channel priority:**

   ```bash
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```

2. **Install AutoDock Vina:**

   ```bash
   conda install vina
   ```

   *Note*: If you encounter issues with the above command, you can alternatively install AutoDock Vina using `pip`:

   ```bash
   pip install vina
   ```

   Ensure that necessary dependencies like `numpy`, `swig`, and `boost` are installed in your environment.

3. **Verify the installation:**

   ```bash
   vina --help
   ```

   This should display the help information for AutoDock Vina.

---

## 🔬 4. Install AutoDockTools (Python 3 Version)

AutoDockTools provides scripts for preparing `.pdbqt` files required by AutoDock Vina.

Install the Python 3-compatible version from GitHub:

```bash
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

Verify the installation:

```bash
which prepare_ligand4.py
which prepare_receptor4.py
```

These scripts should now be accessible in your environment.

---

## 🔄 5. Install Open Babel

Open Babel is a chemical toolbox designed to speak the many languages of chemical data.

Install Open Babel using `conda`:

```bash
conda install -c conda-forge openbabel
```

Verify the installation:

```bash
obabel -V
```

---

## 🧬 6. Install Requirements.txt

To install rest of the lib.

```bash
pip install -r requirements.txt
```
