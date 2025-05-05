

---

### 🧠 Use Case:

1. 🩺 **Input**: A *disease*.
2. 🧬 **Agents**:

   * Identify *target proteins* (from RCSB, UniProt, etc.).
   * Retrieve *ligands* (from PubChem, ChEMBL) known to interact with those targets.
3. 🧪 **Modify ligands** and dock them into **existing receptor structures** (PDBs).
4. 🎯 **Goal**: Discover *new or better protein-ligand interactions* that may suggest **alternative or novel targets** for treatment.

---

### ✅ So What Matters for Ranking?

#### **Your goal is to discover new useful protein targets for a disease.**

That means:

> ✅ **You are ranking *receptors* based on their binding affinity with a modified ligand.**

The assumption is: if a **new protein** (receptor) binds **strongly** to a modified or repurposed ligand, it might be a **novel therapeutic target** — either primary or off-target with benefit.

---

### 🧬 Modified Ligands in Your Pipeline:

* You **start with known ligands** and **modify them** (e.g., via RDKit).
* These modifications aim to improve binding or find new protein matches.
* So the **ligand design** is your way of probing *new protein spaces*.

Thus:

> ✅ **You rank receptors to discover promising targets** — possibly new ones — for the modified ligands.

---

### 🔬 Affinity Interpretation:

* If **modified ligands** bind **stronger** (i.e. more negative affinity) than original ones:

  * The **target receptor is promising**.
  * The **modification is useful**.
* If a **receptor never had a known interaction with the original ligand**, and now shows strong binding:

  * It might be a **new therapeutic target**.

---

### 📌 Final Summary:

| Component            | Role in Your Pipeline                            |
| -------------------- | ------------------------------------------------ |
| **Receptors**        | Potential targets to be **ranked**               |
| **Modified Ligands** | Probes to identify useful receptors              |
| **Affinity Scores**  | Used to **rank receptors** for each ligand       |
| **Goal**             | Find **new proteins** useful in treating disease |


