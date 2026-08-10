import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity


# ==============================
# Load GO terms
# ==============================

csv_path = "/home/jathavp/lab/Neuroprotective-Transcription-factors/all_BP_terms.csv"

df = pd.read_csv(csv_path)

go_terms = df["term_name"].astype(str).tolist()


# ==============================
# Load model
# ==============================

model = SentenceTransformer("all-MiniLM-L6-v2")


# ==============================
# Biological concept anchors
# ==============================
references = {

    "metabolism": [
        "cellular metabolic processes involving biochemical reactions that synthesize, modify, and degrade metabolites",
        "metabolism of carbohydrates, lipids, amino acids, nucleotides, and other cellular molecules",
        "biosynthetic and catabolic pathways that control cellular metabolites and chemical transformations",
        "enzymatic reactions involved in nutrient utilization, energy production, and molecular turnover",
        "regulation of cellular metabolic pathways and metabolic homeostasis",
        "small molecule metabolism and biochemical processes required for cellular function"
    ],


    "lipid_metabolism": [
        "lipid biosynthesis, fatty acid metabolism, cholesterol metabolism, and lipid degradation",
        "cellular regulation of lipid synthesis, storage, transport, and utilization",
        "fatty acid oxidation and lipid-derived energy metabolism",
        "phospholipid metabolism and regulation of lipid homeostasis"
    ],


    "mitochondrial_energy": [
        "mitochondrial respiration, oxidative phosphorylation, and ATP production",
        "electron transport chain activity and mitochondrial energy generation",
        "cellular respiration and mitochondrial metabolic pathways",
        "mitochondrial regulation of energy balance and metabolic activity"
    ],


    "redox_metabolism": [
        "reactive oxygen species metabolism and regulation of cellular redox reactions",
        "oxidoreductase activity, antioxidant metabolism, and reactive molecule processing",
        "superoxide, hydrogen peroxide, and other reactive oxygen species metabolic pathways",
        "glutathione metabolism and cellular redox homeostasis",
        "detoxification and chemical modification of reactive oxygen species"
    ],


    "cell_death": [
        "programmed cell death pathways including apoptosis, necroptosis, and pyroptosis",
        "cellular processes causing loss of viability, elimination, or destruction of cells",
        "regulation of apoptotic signaling pathways and activation of cell death machinery",
        "caspase activation and molecular pathways controlling apoptosis",
        "intrinsic and extrinsic pathways of programmed cell death",
        "clearance of dead and dying cells"
    ],


    "cell_killing": [
        "immune-mediated killing of target cells by cytotoxic lymphocytes",
        "natural killer cell mediated cytotoxicity and T cell mediated killing",
        "leukocyte mediated cytotoxic mechanisms that destroy target cells",
        "complement-mediated cell lysis and immune effector killing"
    ],


    "immune": [
        "immune system activation, inflammation, and defense responses",
        "innate and adaptive immune responses involving leukocytes and immune signaling",
        "cytokine production, chemokine signaling, and inflammatory pathways",
        "immune cell activation, migration, differentiation, and host defense"
    ],


    "transcription": [
        "regulation of gene expression through transcription factors and DNA binding proteins",
        "RNA polymerase activity and transcriptional regulation",
        "control of transcriptional programs and regulation of messenger RNA production",
        "chromatin regulation and gene regulatory networks"
    ],


    "translation": [
        "ribosome assembly, translation, and protein synthesis",
        "production of proteins from messenger RNA",
        "ribosomal function and regulation of translational machinery",
        "peptide synthesis and translation initiation"
    ],


    "signaling": [
        "cellular signaling pathways involving receptors, kinases, second messengers, and signal transduction",
        "regulation of cellular responses through intracellular communication pathways",
        "protein phosphorylation cascades and signaling networks"
    ],


    "transport": [
        "transport of molecules across membranes and between cellular compartments",
        "vesicle trafficking, intracellular localization, and molecular transport",
        "ion transport and regulation of membrane movement"
    ],


    "development_structure": [
        "cell differentiation, tissue development, and organismal development",
        "cell adhesion, cytoskeleton organization, and maintenance of cellular structure",
        "developmental programs controlling cell identity and tissue organization"
    ],


    "neuronal": [
        "neuronal development, axon growth, neurite extension, and synaptic processes",
        "nervous system development and neuronal signaling",
        "neuronal regeneration and maintenance of neural structures"
    ]
}

# ==============================
# Embed GO terms
# ==============================

print("Embedding GO terms...")

term_embeddings = model.encode(
    go_terms,
    show_progress_bar=True
)


# ==============================
# Create category embeddings
# ==============================

print("Embedding biological categories...")

category_embeddings = {}

for category, descriptions in references.items():

    embeddings = model.encode(descriptions)

    # average descriptions into one category vector
    category_embeddings[category] = embeddings.mean(axis=0)


category_embeddings = np.vstack(
    list(category_embeddings.values())
)


categories = list(references.keys())


# ==============================
# Calculate similarity
# ==============================

print("Calculating similarity...")

scores = cosine_similarity(
    term_embeddings,
    category_embeddings
)


score_df = pd.DataFrame(
    scores,
    columns=categories
)


# ==============================
# Combine classifieds
# ==============================

classified = pd.concat(
    [
        df.reset_index(drop=True),
        score_df
    ],
    axis=1
)


# ==============================
# Classification
# ==============================

classified["best_category"] = score_df.idxmax(axis=1)


# confidence = distance between best and second best

classified["confidence"] = score_df.apply(
    lambda row: row.nlargest(2).iloc[0] - row.nlargest(2).iloc[1],
    axis=1
)


# ==============================
# Sort by enrichment
# ==============================

if "negative_log10_of_adjusted_p_value" in classified.columns:
    classified = classified.sort_values(
        "negative_log10_of_adjusted_p_value",
        ascending=False
    )


# ==============================
# Save output
# ==============================

output_path = (
    "/home/jathavp/lab/Neuroprotective-Transcription-factors/"
    "all_BP_terms_semantic_categories.csv"
)

classified.to_csv(
    output_path,
    index=False
)


print("\nSaved:")
print(output_path)

print("\nTop results:")
print(
    classified[
        [
            "term_name",
            "best_category",
            "confidence"
        ]
    ].head(30)
)