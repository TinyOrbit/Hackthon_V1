from custom_tools import *
from custom_agents import *
from crewai import Agent, Task, Crew, Process, LLM

# -------------------- Crew --------------------

load_dotenv()

llm = LLM(
    model="azure/o3-mini",
    api_key=os.environ["AZURE_API_KEY"],
    api_base=os.environ["AZURE_API_BASE"],
    api_version=os.environ["AZURE_API_VERSION"]
)



bioinfo_crew = Crew(
    agents=[target_agent, enrichment_agent, property_agent,molecule_generation_agent,
            structure_generation_agent,pdbqt_conversion_agent,docking_agent,receptor_analyst_agent

            ],
    tasks=[identify_task, enrich_task, prioritize_task
           ,generation_task,structure_task,pdbqt_conversion_task,perform_docking_task,rank_receptors_task
           ],
    process=Process.sequential,
    verbose=True
)

# -------------------- Run Crew --------------------

def run_bioinformatics_pipeline(query):
    print("\n=== Running Bioinformatics Crew ===")
    results = bioinfo_crew.kickoff(inputs={"user_input": query})
    print("\n=== Final Results ===")
    print(results)
    return results

if __name__ == "__main__":
    run_bioinformatics_pipeline("EGFR")
