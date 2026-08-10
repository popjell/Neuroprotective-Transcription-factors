import pandas as pd
import networkx as nx
import obonet
import sys

globals()
# ==========================
# Load your GO enrichment CSV
# ==========================
if "input_df" in globals():
    traversal_df = input_df
    print("Running via R integration.")
    del input_df
else:
    input_file = "/home/jathavp/lab/Neuroprotective-Transcription-factors/all_BP_terms.csv"
    traversal_df = pd.read_csv(input_file)
    print("running_manually")


# ==========================
# Load GO graph
# ==========================

print("Loading GO graph...")

graph = obonet.read_obo(
    "https://current.geneontology.org/ontology/go-basic.obo"
)

# GO edges are child -> parent
# Reverse so we can travel parent -> children
go_graph = graph.reverse()

print("Loaded", len(go_graph.nodes), "GO terms")


# ==========================
# Define broad categories
# ==========================

anchors = {

    "metabolism": [

        "GO:0008152",

    ],


    "exclude": [

        "GO:0001816", # cytokine production, technically a metabolic process, but we want to exclude it,
        "GO:0010467"

    ],


    "cell_death": [

        "GO:0008219",   # cell death
        "GO:0001906",   # cell killing
        "GO:0031341"    # regulation of cell killing

    ]
}


# ==========================
# Validate anchors
# ==========================

print("\nChecking anchors:")

for category, ids in anchors.items():

    for go_id in ids:

        if go_id not in go_graph:
            print("Missing:", category, go_id)

        else:
            print(
                "OK:",
                category,
                go_id,
                go_graph.nodes[go_id].get("name")
            )


# ==========================
# Classify GO terms
# ==========================

def classify_term(go_id):

    if go_id not in go_graph:
        return ["unknown"]

    matches = []

    for category, anchor_list in anchors.items():

        for anchor in anchor_list:

            if anchor not in go_graph:
                continue

            if (
                go_id == anchor
                or go_id in nx.descendants(go_graph, anchor)
            ):
                matches.append(category)
                break

    return matches if matches else ["other"]



# ==========================
# Apply classification
# ==========================

traversal_df["raw_category"] = traversal_df["term_id"].apply(
    classify_term
)



traversal_df["category"] = traversal_df["raw_category"].apply(lambda x: ",".join(x))
output_file = "all_bp_terms_graphtraversal.csv"

traversal_df.to_csv(
    output_file,
    index=False
)


print("\nFinished!")
print("Saved:", output_file)



print(
    traversal_df["category"].value_counts()
)