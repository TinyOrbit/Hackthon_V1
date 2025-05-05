# **Team Name: Flagbearers**

## 🧬 Overview

Our project tackles the challenge of discovering and prioritizing new drug targets by using an automated, multi-agent in silico modeling pipeline. The system starts from a disease name, identifies relevant protein targets, retrieves known ligands, modifies them, docks them into protein structures, and ranks the resulting protein-ligand interactions to identify promising therapeutic avenues.

## 🔍 Use Case Summary

### 🩺 **Input**

Accepts any user-provided **disease name** to initiate the in-silico discovery workflow.

### 🧠 **Pipeline Agents & Roles**

1. **Target Identification Agent**

   * Retrieves disease-associated protein targets from databases such as UniProt and RCSB PDB.

2. **Enrichment Agent**

   * Fetches ligands known to bind those targets from sources like PubChem and ChEMBL.

3. **Property Prioritization Agent(Dummy)**

   * Applies filters such as Lipinski’s Rule of Five, QED score, and molecular weight thresholds to prioritize drug-like compounds.

4. **Molecule Generation Agent**

   * Modifies ligands using RDKit to generate structural analogs for improved or novel binding potential.

5. **Structure Generation Agent**

   * Converts 2D SMILES or mol blocks into 3D structures using RDKit and prepares proteins and ligands for docking.

6. **PDBQT Conversion Agent**

   * Converts prepared 3D ligand and receptor files into PDBQT format using Open Babel for compatibility with AutoDock Vina.

7. **Docking Agent**

   * Performs molecular docking using AutoDock Vina and returns binding affinity scores.

## 🎯 Project Goal

To identify **novel or alternative protein targets** that demonstrate strong binding affinity with modified ligands, suggesting their potential relevance in treating the specified disease.

The **receptors** are ranked by their binding strength to modified ligands, aiming to:

* Repurpose ligands to new targets.
* Discover new protein targets.
* Improve original ligand efficacy.

## 🧬 Implementation Summary

* **Framework**: [CrewAI](https://www.crewai.com/) used to orchestrate modular agents in a sequential workflow.

* **Tools Used**:

  * **UniProt / RCSB**: Protein target retrieval
  * **PubChem / ChEMBL**: Ligand enrichment
  * **RDKit**: Molecule modification, 3D structure generation
  * **Open Babel**: Format conversion
  * **AutoDock Vina**: Docking and affinity scoring

* **Programming Language**: Python (with shell calls to external tools where required)

* **Input**: A disease name (e.g., "EGFR")

* **Output**: A ranked list of protein targets based on their docking scores with modified ligands

## 💡 Intended Users

* Drug discovery researchers seeking novel protein targets or repurposing opportunities
* Computational chemists designing or refining small molecule libraries

## 🏭  Use Case
### 🧬 Life Science
  🔬 In-Silico Modelling

## 📊 Example Scenarios

* **Virtual Screening**: Rank large databases of molecules to identify high-affinity binders.
* **Drug Repurposing**: Test known ligands against new protein targets.
* **Lead Optimization**: Compare ΔAffinity of modified ligands to improve interaction strength.

## 👥 Contributors

* **\[Rahul Suresh Shedge]** – \[AI Enginer/GenAI Developer]
* **\[To be completed: Contributor 2 Name]** – \[Role]
* **\[To be completed: Contributor 3 Name]** – \[Role]
* **\[To be completed: Contributor 4 Name]** – \[Role]

## 📸 Screenshots

*
*

## 📈 Future Enhancements

* **Drug-likeness & ADMET-aware Ranking**: Rank receptors not only by docking scores but also considering ligand ADMET and drug-likeness.
* **Receptor-Level ΔAffinity Analysis**: Highlight receptors that gained high affinity only after ligand modification.
* **Disease Relevance Scoring**: Combine DisGeNET, UniProt annotations with docking scores to prioritize actionable targets.
* **Automation of PDB Retrieval & Preprocessing**: Integrate steps to auto-download and prepare receptors.
* **User Interface**: Develop a GUI or web front end for easy interaction.
* **Parallel Docking Execution**: Improve runtime for large ligand libraries.

---

