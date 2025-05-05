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

class UniProtToolSchema(BaseModel):
    target_name: str

class UniProtTool(BaseTool):
    name: str = "UniProt Fetcher"
    description: str = "Fetch protein information from UniProt given a target name."
    args_schema = UniProtToolSchema

    def _run(self, target_name):
        print(f"\n[UniProtTool] Querying UniProt for: {target_name}")
        url = f"https://rest.uniprot.org/uniprotkb/search?query=({target_name})AND(organism_name:human)&format=json"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            if data and "results" in data and len(data["results"]) > 0:
                first_result = data["results"][0]
                uniprot_id = first_result.get("primaryAccession")
                sequence = first_result.get("sequence", {}).get("value", "")
                protein_name = first_result.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", target_name)
                alternative_names = [name.get("value") for name in first_result.get("protein", {}).get("alternativeName", []) if name.get("value")]
                pdb_entries = first_result.get("uniProtKBCrossReferences", [])
                pdb_ids = [entry.get("id") for entry in pdb_entries if entry.get("database") == "PDB"]

                output = {
                    "target_name": protein_name,
                    "alternative_names": alternative_names,
                    "uniprot_id": uniprot_id,
                    "sequence": sequence,
                    "potential_pdb_ids": pdb_ids
                }
                print(f"[UniProtTool] Output: {json.dumps(output, indent=2)}")
                return output
            else:
                print(f"[UniProtTool] No results found for {target_name}")
                return None
        except Exception as e:
            print(f"[UniProtTool] Error: {e}")
            return None

class BioinformaticsEnrichmentToolSchema(BaseModel):
    target_name: str
    uniprot_id: str
    pdb_ids: list[str] = []
    file_format: str = 'pdb'

class BioinformaticsEnrichmentTool(BaseTool):
    name: str = "Bioinformatics Enrichment Tool"
    description: str = "Fetch compound info from PubChem, bioactivity from ChEMBL, and 3D structure from PDB."
    args_schema = BioinformaticsEnrichmentToolSchema

    def _run(self, target_name, uniprot_id, pdb_ids=[], file_format='pdb'):
        result = {
            "pubchem_compounds": [],
            "chembl_bioactivities": [],
            "pdb_structures": []
        }

        # --- PubChem Part ---
        try:
            print(f"\n[PubChem] Querying PubChem for: {target_name}, {uniprot_id}")
            url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{target_name}/cids/JSON"
            response = requests.get(url_cid)
            response.raise_for_status()
            cid_data = response.json()

            cids = cid_data.get('IdentifierList', {}).get('CID', [])
            compounds = []
            for cid in cids:
                url_xrefs = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/JSON"
                xref_response = requests.get(url_xrefs)
                xref_response.raise_for_status()
                xref_data = xref_response.json()

                information_list = xref_data.get('InformationList', {}).get('Information', [])
                is_linked = any('UniProt' in info and uniprot_id in info['UniProt'] for info in information_list)

                if is_linked:
                    smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
                    smiles_response = requests.get(smiles_url)
                    smiles_response.raise_for_status()
                    smiles_data = smiles_response.json()
                    smiles = smiles_data.get('PropertyTable', {}).get('Properties', [{}])[0].get('CanonicalSMILES')
                    if smiles:
                        compounds.append({
                            "smiles": smiles,
                            "cid": cid,
                            "name": f"Compound_{cid}"
                        })
            result["pubchem_compounds"] = compounds
            print(f"[PubChem] Found {len(compounds)} compounds.")
        except Exception as e:
            print(f"[PubChem] Error: {e}")

        # --- ChEMBL Part ---
        try:
            print(f"\n[ChEMBL] Querying ChEMBL for UniProt ID: {uniprot_id}")
            target_url = f"https://www.ebi.ac.uk/chembl/api/data/target?target_components.accession={uniprot_id}&format=json"
            target_response = requests.get(target_url)
            target_response.raise_for_status()
            target_data = target_response.json()
            targets = target_data.get("targets", [])
            if targets:
                target_chembl_id = targets[0].get("target_chembl_id")
                print(f"[ChEMBL] Found target_chembl_id: {target_chembl_id}")
                activity_url = f"https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id={target_chembl_id}&format=json"
                activity_response = requests.get(activity_url)
                activity_response.raise_for_status()
                activity_data = activity_response.json()
                activities = activity_data.get("activities", [])

                chembl_output = []
                for activity in activities:
                    molecule_chembl_id = activity.get("molecule_chembl_id")
                    if not molecule_chembl_id:
                        continue
                    molecule_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{molecule_chembl_id}?format=json"
                    molecule_response = requests.get(molecule_url)
                    if molecule_response.status_code != 200:
                        continue
                    molecule_data = molecule_response.json()
                    smiles = molecule_data.get("molecule_structures", {}).get("canonical_smiles")
                    if smiles:
                        chembl_output.append({
                            "chembl_id": molecule_chembl_id,
                            "smiles": smiles,
                            "target_chembl_id": activity.get("target_chembl_id"),
                            "activity": activity.get("standard_type"),
                            "value": activity.get("standard_value"),
                            "unit": activity.get("standard_units")
                        })
                result["chembl_bioactivities"] = chembl_output
                print(f"[ChEMBL] Found {len(chembl_output)} bioactivities.")
            else:
                print(f"[ChEMBL] No target found for UniProt ID: {uniprot_id}")
        except Exception as e:
            print(f"[ChEMBL] Error: {e}")

        import os
        import glob

        # Create the folder if it doesn't exist
        os.makedirs('pdb_files', exist_ok=True)

        # Delete all files in the folder if it exists
        if os.path.exists('pdb_files'):
            for f in glob.glob('./pdb_files/*'):
                if os.path.isfile(f):
                    os.remove(f)

        for pdb_id in pdb_ids:
            # pass
            pdb_id = pdb_id.upper()
            print(f"\n[PDB] Retrieving structure for PDB ID: {pdb_id}")
            if not re.match(r'^[A-Z0-9]{4}$', pdb_id):
                print(f"[PDB] Invalid PDB ID: {pdb_id}")
                continue

            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            pdb_file_path = f'pdb_files/{pdb_id}.pdb'
            try:
                response = requests.get(url, stream=True)
                response.raise_for_status()
                with open(pdb_file_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                result["pdb_structures"].append({
                    "pdb_id": pdb_id,
                    "file_path": pdb_file_path
                })
                print(f"[PDB] Saved structure for {pdb_id} at {pdb_file_path}.")
            except requests.exceptions.HTTPError as http_err:
                print(f"[PDB] HTTP error for {pdb_id}: {http_err}")
            except Exception as err:
                print(f"[PDB] Other error for {pdb_id}: {err}")

        return result



class PropertyPrioritizationSchema(BaseModel):
    pubchem_compounds: list = []
    chembl_bioactivities: list = []

class PropertyPrioritizationTool(BaseTool):
    name: str = "Property Prioritization Tool"
    description: str = "Analyze molecular properties and prioritize based on drug-likeness criteria."
    args_schema : type = PropertyPrioritizationSchema

    def _run(self, pubchem_compounds=None, chembl_bioactivities=None):
        print("\n[PropertyPrioritization] Evaluating molecular properties...")

        properties = {
            "molecular_weight_range": [200, 500],
            "logp_range": [0.5, 4.5],
            "hydrogen_bond_donors_range": [0, 5],
            "hydrogen_bond_acceptors_range": [2, 8],
            "num_rings": None
        }

        # Try PubChem compounds first, then ChEMBL
        smiles = None
        if pubchem_compounds:
            for c in pubchem_compounds:
                smiles = c.get("smiles")
                if smiles:
                    break
        if not smiles and chembl_bioactivities:
            for b in chembl_bioactivities:
                smiles = b.get("smiles")
                if smiles:
                    break

        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                properties["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
                # print(properties["num_rings"])

        print(f"[PropertyPrioritization] Output: {json.dumps(properties, indent=2)}")
        return {"prioritized_properties": properties}


class MoleculeStructure(BaseModel):
    original: str = Field(..., description="Original SMILES string")
    modified: str = Field(..., description="Modified SMILES string")

class MoleculeStructureList(BaseModel):
    molecules: List[MoleculeStructure]

class MoleculeStructurewithMOL(BaseModel):
    original: str = Field(..., description="Original SMILES string")
    modified: str = Field(..., description="Modified SMILES string")
    mol_block: Optional[str] = None 

class MoleculeStructureListwithMOL(BaseModel):
    molecules: List[MoleculeStructurewithMOL]

class MoleculeGenerationSchema(BaseModel):
    smiles_list: list

class MoleculeGenerationTool(BaseTool):
    name: str = "Molecule Generation Tool"
    description: str = "Modify known ligands using RDKit to generate novel analogs."
    args_schema: type = MoleculeGenerationSchema

    def _run(self, smiles_list):
        modified_molecules = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            # Add explicit hydrogens before modification
            mol_with_H = Chem.AddHs(mol)
            hydrogen_idxs = [atom.GetIdx() for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 1]

            if not hydrogen_idxs:
                continue  # No replaceable hydrogens

            # Replace one H with a methyl group (-CH3)
            replace_idx = random.choice(hydrogen_idxs)
            editable = Chem.RWMol(mol_with_H)
            editable.ReplaceAtom(replace_idx, Chem.Atom("C"))

            # Sanitize and get SMILES of modified molecule
            modified_mol = editable.GetMol()
            try:
                Chem.SanitizeMol(modified_mol)
                modified_smiles = Chem.MolToSmiles(modified_mol)
                modified_molecules.append({
                    "original": smiles,
                    "modified": modified_smiles
                })
            except:
                continue  # Skip if invalid

        return modified_molecules
        # return MoleculeStructureList(molecules=modified_molecules)


from rdkit.Chem import AllChem

# Define the schema for molecules
class MoleculeStructure(BaseModel):
    original: str  # original SMILES
    modified: str  # modified SMILES

# Define the schema for the tool's input
class StructureGenerationSchema(BaseModel):
    molecules_list: List[MoleculeStructure]  # List of MoleculeStructure

# Tool class for 3D structure generation
class StructureGenerationTool(BaseTool):
    name: str = "Structure Generation Tool"
    description: str = "Generate 3D structures from modified SMILES using RDKit."
    args_schema: type = StructureGenerationSchema

    # Run method to generate 3D structures
    def _run(self, molecules_list: List[Dict[str, str]]):
        molecules_3d = []
        # Iterate through the list of molecule data
        for mol_data in molecules_list:
            print("Processing molecule:")
            print(mol_data)

            smiles = mol_data.get("modified")
            if not smiles:
                continue

            # Convert SMILES to RDKit molecule object
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f"Invalid SMILES: {smiles}")
                continue

            # Add hydrogens and perform 3D embedding
            mol = Chem.AddHs(mol)
            try:
                success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                if success == 0:
                    AllChem.UFFOptimizeMolecule(mol)
                    mol_block = Chem.MolToMolBlock(mol)
                    molecules_3d.append({
                        "original": mol_data.get("original"),
                        "modified": smiles,
                        "mol_block": mol_block
                    })
                else:
                    print(f"3D embedding failed for {smiles}")
            except Exception as e:
                print(f"[3D Generation Error] for {smiles}: {e}")

        return molecules_3d









# class MoleculeStructure(BaseModel):
#     original: str = Field(..., description="Original SMILES string")
#     modified: str = Field(..., description="Modified SMILES string")

# class StructureGenerationSchema(BaseModel):
#     molecules_list: List[MoleculeStructure]  # list of {"original": ..., "modified": ...}

# class StructureGenerationTool(BaseTool):
#     name: str = "Structure Generation Tool"
#     description: str = "Generate 3D structures from modified SMILES using RDKit."
#     args_schema: type = StructureGenerationSchema

#     # def _run(self, molecules_list: List[MoleculeStructure]):
#     def _run(self, molecules_list):
#         molecules_3d = []
#         # print(">>>"*20)
#         # print(molecules_list)
#         # print(">>>"*20)
#         for mol_data in molecules_list:
#             print("##"*30)
#             print(mol_data)
#             smiles = mol_data.get("modified") 

#             print("**"*30)
#             print(smiles)
#             print("**"*20)
#             mol = Chem.MolFromSmiles(smiles)
#             if not mol:
#                 continue
#             mol = Chem.AddHs(mol)
#             try:
#                 success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
#                 if success == 0:
#                     AllChem.UFFOptimizeMolecule(mol)
#                     mol_block = Chem.MolToMolBlock(mol)
#                     molecules_3d.append({
#                         "original": mol_data.get("original"),
#                         "modified" : smiles,
#                         "mol_block" : mol_block  
#                     })
#             except Exception as e:
#                 print(f"[3D Generation Error] for {smiles}: {e}")
#         # return MoleculeStructureListwithMOL(molecules_3d)
#         return molecules_3d


class ReceptorInput(BaseModel):
    pdb_id: str
    file_path: str

class LigandInput(BaseModel):
    original: str
    modified: str
    mol_block: str

class PDBQTConversionInput(BaseModel):
    pdb_structures: List[ReceptorInput] = Field(..., description="List of receptor PDB structures.")
    ligands: List[LigandInput] = Field(..., description="List of ligands with mol_block data.")

class PDBQTConversionTool(BaseTool):
    name:str = "PDBQT Conversion Tool"
    description: str = "Converts receptor PDB files and ligand mol_blocks into .pdbqt format using MGLTools."
    args_schema: type = PDBQTConversionInput

    def _run(self, pdb_structures: List[dict], ligands: List[dict]) -> dict:
        receptors = pdb_structures
        ligands = ligands
        print("***"*20)
        print("Recep :",receptors)
        print("##"*20)
        print("Ligands : ",ligands)

        output = {"receptors": [], "ligands": []}
        receptor_dir = "receptor_pdbqt"
        ligand_dir = "ligand_pdbqt"
        os.makedirs(receptor_dir, exist_ok=True)
        os.makedirs(ligand_dir, exist_ok=True)

        # Convert Receptors
        for rec in receptors:
            pdb_id = rec["pdb_id"]
            pdb_path = rec["file_path"]
            # pdb_id = rec.pdb_id
            # pdb_path = rec.file_path
            pdbqt_path = os.path.join(receptor_dir, f"{pdb_id}.pdbqt")
            try:
                subprocess.run([
                    "prepare_receptor4",
                    "-r", pdb_path,
                    "-o", pdbqt_path,
                    "-A", "hydrogens"
                ], check=True)

                output["receptors"].append({
                    "pdb_id": pdb_id,
                    "pdb_file": pdb_path,
                    "pdbqt_file": pdbqt_path
                })

            except Exception as e:
                print(e)
        # Convert Ligands
        for i, ligand in enumerate(ligands):
            ligand_id = ligand["original"] or f"lig_{i}"
            mol_block = ligand["mol_block"]
            # ligand_id = ligand.original or f"lig_{i}"
            # mol_block = ligand.mol_block
            mol_path = os.path.join(ligand_dir, f"{ligand_id}.mol")
            pdb_path = mol_path.replace(".mol", ".pdb")
            pdbqt_path = mol_path.replace(".mol", ".pdbqt")

            try:
                # Save mol_block and convert to PDB using Open Babel
                with open(mol_path, "w") as f:
                    f.write(mol_block)

                subprocess.run(["obabel", mol_path, "-O", pdb_path], check=True)

                subprocess.run([
                    "prepare_ligand4",
                    "-l", pdb_path,
                    "-o", pdbqt_path,
                    "-A", "hydrogens"
                ], check=True)

                output["ligands"].append({
                    "ligand_id": ligand_id,
                    "mol_file": mol_path,
                    "pdb_file": pdb_path,
                    "pdbqt_file": pdbqt_path
                })
            except Exception as e:
                print(e)

        return output



class DockingInput(BaseModel):
    receptors: List[str] = Field(..., description="Paths to receptor .pdbqt files.")
    ligands: List[str] = Field(..., description="Paths to ligand .pdbqt files.")
    center_x: float = Field(10.0, description="X-coordinate of the docking box center.")
    center_y: float = Field(12.5, description="Y-coordinate of the docking box center.")
    center_z: float = Field(15.0, description="Z-coordinate of the docking box center.")
    size_x: float = Field(20.0, description="Size of the docking box along the X-axis.")
    size_y: float = Field(20.0, description="Size of the docking box along the Y-axis.")
    size_z: float = Field(20.0, description="Size of the docking box along the Z-axis.")
    exhaustiveness: int = Field(8, description="Exhaustiveness of the global search.")
    num_modes: int = Field(9, description="Maximum number of binding modes to generate.")


class AutoDockVinaTool(BaseTool):
    name :str = "AutoDock Vina Tool"
    description : str = "Performs molecular docking using AutoDock Vina."
    args_schema: Type[BaseModel] = DockingInput

    def _run(
        self,
        receptors: List[str],
        ligands: List[str],
        center_x: float,
        center_y: float,
        center_z: float,
        size_x: float,
        size_y: float,
        size_z: float,
        exhaustiveness: int,
        num_modes: int
    ) -> dict:
        output = {"docking_results": []}
        output_dir = "docking_results"
        os.makedirs(output_dir, exist_ok=True)

        print("Receptors : : : :",receptors)
        print("ligands : : : :",ligands)

        for receptor_path in receptors:
            receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
            for ligand_path in ligands:
                ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
                result_prefix = f"{receptor_name}_{ligand_name}"
                out_pdbqt = os.path.join(output_dir, f"{result_prefix}_out.pdbqt")
                log_file = os.path.join(output_dir, f"{result_prefix}_log.txt")

                try:
                    subprocess.run([
                        "vina",
                        "--receptor", receptor_path,
                        "--ligand", ligand_path,
                        "--center_x", str(center_x),
                        "--center_y", str(center_y),
                        "--center_z", str(center_z),
                        "--size_x", str(size_x),
                        "--size_y", str(size_y),
                        "--size_z", str(size_z),
                        "--exhaustiveness", str(exhaustiveness),
                        "--num_modes", str(num_modes),
                        "--out", out_pdbqt,
                        # "--log", log_file
                    ], check=True)

                    output["docking_results"].append({
                        "receptor": receptor_path,
                        "ligand": ligand_path,
                        "output_pdbqt": out_pdbqt,
                        # "log": log_file
                    })
                except Exception as e:
                    print(e)

        return output


# Update ReceptorRankingInput
class ReceptorRankingInput(BaseModel):
    docking_results: List[dict] = Field(..., description="List of docking results with receptor and ligand information.")

# Update ReceptorRankingTool
class ReceptorRankingTool(BaseTool):
    name: str = "Receptor Ranking Tool"
    description: str = "Ranks receptors based on their binding affinity with ligands extracted from docking results."
    args_schema: Type[BaseModel] = ReceptorRankingInput

    def _run(self, docking_results: List[dict]) -> List[dict]:
        ranked_receptors = []
        for result in docking_results:
            receptor_path = result.get('receptor')
            ligand_path = result.get('ligand')
            output_pdbqt = result.get('output_pdbqt')
            receptor_id = os.path.basename(receptor_path).split('.')[0] if receptor_path else 'unknown_receptor'
            affinity = None

            if output_pdbqt and os.path.exists(output_pdbqt):
                try:
                    with open(output_pdbqt, 'r') as f:
                        for line in f:
                            if line.startswith('REMARK VINA RESULT:'):
                                parts = line.strip().split()
                                if len(parts) >= 4:
                                    affinity = float(parts[3])
                                    break
                except Exception as e:
                    print(f"Error reading output PDBQT file {output_pdbqt}: {e}")

            if affinity is not None:
                ranked_receptors.append({
                    'receptor_id': receptor_id,
                    'binding_affinity': affinity
                })

        # Sort receptors by binding affinity (more negative is better)
        ranked_receptors.sort(key=lambda x: x['binding_affinity'])
        return ranked_receptors
