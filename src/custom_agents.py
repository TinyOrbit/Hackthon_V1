
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

# -------------------- Agents --------------------

target_agent = Agent(
    role="Protein Target Identifier",
    goal="Identify protein information from UniProt.",
    backstory="An expert bioinformatician skilled in protein databases.",
    tools=[UniProtTool()],
        #    , FileReadTool()],  # Add FileReadTool
    llm=llm
)

enrichment_agent = Agent(
    role="Protein Information Enricher",
    goal="Find inhibitors, bioactivity data, and 3D structures for the given protein.",
    backstory="Specializes in drug discovery and protein structure retrieval.",
    tools=[BioinformaticsEnrichmentTool()],
    #   FileReadTool()],  # Add FileReadTool
    llm=llm
)

property_agent = Agent(
    role="Molecular Property Prioritizer",
    goal="Evaluate and prioritize molecular properties to guide compound selection.",
    backstory="A chemoinformatics expert using drug-likeness heuristics to guide design decisions.",
    tools=[PropertyPrioritizationTool()],
    llm=llm
)


molecule_generation_agent = Agent(
    role="Ligand Designer",
    goal="Create new ligands by modifying existing molecules from PubChem or ChEMBL.",
    backstory="A molecular chemist with expertise in synthetic drug design and structure-activity relationships.",
    tools=[MoleculeGenerationTool()],
    llm=llm,
    # output_pydantic=
)


structure_generation_agent = Agent(
    role="3D Structure Generator",
    goal="Convert 2D molecules to 3D-optimized structures using cheminformatics tools. Use the result of generation_task which contains list of smiles original & modified.",
    backstory="A computational chemist skilled in conformer generation and molecular geometry optimization.",
    tools=[StructureGenerationTool()],
    llm=llm
)


pdbqt_conversion_agent = Agent(
    role="Molecule Docking Preparation Agent",
    goal="Convert ligand and receptor structures to .pdbqt format using AutoDockTools",
    backstory=(
        "An expert in molecular docking prep, this agent ensures all ligands and receptors "
        "are properly converted to the .pdbqt format for AutoDock simulations."
    ),
    verbose=True,
    llm = llm,
    tools=[PDBQTConversionTool()]
)

docking_agent = Agent(
    name="DockingAgent",
    role="Molecular Docking Specialist",
    goal="Perform molecular docking simulations to evaluate ligand binding affinities.",
    backstory=(
        "An expert in computational chemistry, proficient in using AutoDock Vina "
        "for simulating ligand-receptor interactions to aid in drug discovery."
    ),
    llm = llm,
    tools=[AutoDockVinaTool()]
)

receptor_analyst_agent = Agent(
    name="ReceptorAnalyst",
    role="Receptor Evaluation Specialist",
    goal="Evaluate and rank receptors based on their binding affinity with ligands to identify the most promising therapeutic targets.",
    backstory=(
        "An expert in evaluating receptor-ligand interactions, focused on identifying receptors with the highest binding affinities "
        "to guide drug discovery and development for disease treatment."
    ),
    llm=llm,
    tools=[ReceptorRankingTool()]  # Using ReceptorRankingTool to rank receptors
)
