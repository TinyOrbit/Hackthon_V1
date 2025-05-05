
import re
import os
import json
import random 
import requests
import tempfile
import subprocess
from rdkit import Chem
from typing import List,Type
from pydantic import BaseModel
from dotenv import load_dotenv
from rdkit.Chem import AllChem
from typing import List,Optional
from crewai.tools import BaseTool
from pydantic import BaseModel, Field
from rdkit.Chem import Descriptors, rdMolDescriptors
from crewai import Agent, Task, Crew, Process, LLM
from custom_tools import *
from custom_agents import *


identify_task = Task(
    description="Identify the protein target from the user's query and fetch UniProt data.",
    expected_output="Protein target information including UniProt ID, sequence, and potential PDB IDs.",
    agent=target_agent,
)

enrich_task = Task(
    description="Using the UniProt data, find inhibitors from PubChem, bioactivity from ChEMBL, and download structures from PDB. Use FileReadTool to read structures if needed.",
    expected_output="Enriched data with known inhibitors, bioactivity results, and 3D structure file paths.",
    agent=enrichment_agent,
    context=[identify_task]
)

prioritize_task = Task(
    description="Prioritize molecular properties using SMILES data collected from PubChem and ChEMBL sources, applying RDKit descriptors. If PubChem compounds are unavailable, fallback to ChEMBL bioactivities.",
    expected_output="A dictionary of prioritized drug-like properties including molecular weight range, logP range, H-bond donor/acceptor counts, and ring counts.",
    agent=property_agent,
    context=[enrich_task])

generation_task = Task(
    description=(
        "Using known ligand SMILES retrieved from PubChem or ChEMBL, create structurally related analogs "
        "by applying simple modifications such as replacing hydrogen with a methyl group."
    ),
    expected_output="A list of modified SMILES with corresponding original molecules. The output has to be list containing dictionries of modified and its orginal smile for every smile pass.",
    agent=molecule_generation_agent,
    context=[enrich_task,prioritize_task],
)

structure_task = Task(
    description="Take a list of orginal and modified smiles.Generate 3D structures for the modified ligands using RDKit, including energy minimization and conformer embedding.",
    expected_output="A list of 3D-optimized molecules with MolBlock data.",
    agent=structure_generation_agent,
    context=[generation_task]
)


pdbqt_conversion_task = Task(
    description=(
        "Convert ligands and receptor structures to .pdbqt format using AutoDockTools. "
        "Use ligands from the structure generation task and receptors from the enrichment task."
        "ligands is the list of dictionaries contains keys like original,modified and mol_block , Pass entire list as input."
        "receptors is the list of dictornaries contains keys like pubchem_compounds,chembl_bioactivities,pdb_structures, pass only pdb_structures as input."
    ),
    expected_output="A dictionary with 'ligands' and 'receptors' keys, each containing .pdbqt file paths.",
    agent=pdbqt_conversion_agent,
    context=[structure_task,enrich_task
    ]
)

perform_docking_task = Task(
    name="PerformDocking",
    description=( "Execute molecular docking simulations using AutoDock Vina."  ),
    agent=docking_agent,
    context=[pdbqt_conversion_task],
    expected_output="Dictionary containing docking results for each receptor-ligand pair."
)



# Update Task for Ranking Receptors
rank_receptors_task = Task(
    name="RankReceptors",
    description="Rank receptors based on their binding affinity with ligands extracted from docking output PDBQT files.",
    agent=receptor_analyst_agent,  # Updated agent to ReceptorAnalyst
    context=[perform_docking_task],  # Assuming perform_docking_task is defined elsewhere
    expected_output="A list of receptors ranked by their binding affinity with ligands."
)
