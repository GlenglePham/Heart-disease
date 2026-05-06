# simulate_1000genotypes.py

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Create output directories
os.makedirs("data", exist_ok=True)
os.makedirs("pic", exist_ok=True)

# Load processed SNP data
df = pd.read_csv("data/cad_snps_processed.csv")

# Number of simulated individuals
N = 1000

# Create genotype dataframe
genotype_df = pd.DataFrame()

# Create individual IDs
genotype_df["person_id"] = [
    f"P{i+1}" for i in range(N)
]

# --------------------------------------------------
# Simulate genotypes using Hardy-Weinberg equilibrium
# --------------------------------------------------

for _, row in df.iterrows():

    rsid = row["rsID"]
    raf = row["RAF"]

    # Allele frequencies
    q = raf
    p = 1 - q

    # Genotype probabilities
    probs = [
        p**2,       # genotype 0
        2 * p * q,  # genotype 1
        q**2        # genotype 2
    ]

    # Generate synthetic genotypes
    genotypes = np.random.choice(
        [0, 1, 2],
        size=N,
        p=probs
    )

    genotype_df[rsid] = genotypes

# --------------------------------------------------
# Calculate PRS
# --------------------------------------------------

# Initialize PRS column
genotype_df["PRS"] = 0.0

# Calculate PRS for each individual
for _, row in df.iterrows():

    rsid = row["rsID"]
    beta = row["beta_final"]

    genotype_df["PRS"] += (
        genotype_df[rsid] * beta
    )

# Print simulated data
print("\nSimulated genotype matrix with PRS:\n")
print(genotype_df.head())

# Print PRS summary statistics
print("\nPRS summary:\n")
print(genotype_df["PRS"].describe())

# # Save simulated genotype dataset
# genotype_df.to_csv(
#     "data/simulated_genotypes.csv",
#     index=False
# )

# print("\nSaved: data/simulated_genotypes.csv")

# --------------------------------------------------
# PRS distribution visualization
# --------------------------------------------------

plt.figure(figsize=(10, 6))

plt.hist(
    genotype_df["PRS"],
    bins=30
)

plt.xlabel("Polygenic Risk Score (PRS)")
plt.ylabel("Number of Individuals")
plt.title("Distribution of CAD Polygenic Risk Scores")

plt.grid(True)

plt.savefig(
    "pic/prs_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nSaved: pic/prs_distribution.png")

# --------------------------------------------------
# Top / Bottom PRS analysis
# --------------------------------------------------

# Define percentile thresholds
top_5_threshold = genotype_df["PRS"].quantile(0.95)
bottom_5_threshold = genotype_df["PRS"].quantile(0.05)

# Select high-PRS individuals
top_group = genotype_df[
    genotype_df["PRS"] >= top_5_threshold
]

# Select low-PRS individuals
bottom_group = genotype_df[
    genotype_df["PRS"] <= bottom_5_threshold
]

print("\nTop 5% PRS individuals:\n")
print(
    top_group[
        ["person_id", "PRS"]
    ].head()
)

print("\nBottom 5% PRS individuals:\n")
print(
    bottom_group[
        ["person_id", "PRS"]
    ].head()
)

print("\nTop 5% threshold:", top_5_threshold)
print("Bottom 5% threshold:", bottom_5_threshold)

# --------------------------------------------------
# Top vs Bottom genotype enrichment analysis
# --------------------------------------------------

comparison_list = []

for _, row in df.iterrows():

    rsid = row["rsID"]

    # Mean genotype dosage
    top_mean = top_group[rsid].mean()
    bottom_mean = bottom_group[rsid].mean()

    # Difference between groups
    difference = top_mean - bottom_mean

    comparison_list.append({
        "rsID": rsid,
        "gene": row["gene"],
        "top_mean_genotype": top_mean,
        "bottom_mean_genotype": bottom_mean,
        "difference": difference
    })

# Create dataframe
comparison_df = pd.DataFrame(comparison_list)

# Sort by enrichment difference
comparison_df = comparison_df.sort_values(
    by="difference",
    ascending=False
)

# Print results
print("\nTop vs Bottom PRS genotype comparison:\n")
print(comparison_df.head(10))

# --------------------------------------------------
# Visualization
# --------------------------------------------------

plt.figure(figsize=(12, 6))

plt.bar(
    comparison_df["rsID"].head(10),
    comparison_df["difference"].head(10)
)

plt.xlabel("SNP")
plt.ylabel("Mean Genotype Difference")

plt.title(
    "Genotype Enrichment in Top 5% PRS Individuals"
)

plt.xticks(rotation=45)

plt.grid(True)

plt.savefig(
    "pic/top_vs_bottom_prs_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(
    "\nSaved: pic/top_vs_bottom_prs_comparison.png"
)