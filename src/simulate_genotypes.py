# simulate_genotypes.py

import pandas as pd
import numpy as np

# Load processed data
df = pd.read_csv("data/cad_snps_processed.csv")

# Number of individuals
N = 1000

# Create empty genotype dataframe
genotype_df = pd.DataFrame()

# Add individual IDs
genotype_df["person_id"] = [f"P{i+1}" for i in range(N)]

# Simulate genotypes for each SNP
for _, row in df.iterrows():
    rsid = row["rsID"]
    raf = row["RAF"]
    # Hardy-Weinberg frequencies
    q = raf
    p = 1 - q
    probs = [p**2, 2*p*q, q**2]
    # p**2,       : genotype 0
    # 2 * p * q,  : genotype 1
    # q**2        : genotype 2

    # Generate synthetic genotypes
    genotypes = np.random.choice([0, 1, 2], size=N, p=probs)
    genotype_df[rsid] = genotypes

# Check
print(genotype_df.head())


# Initialize PRS column
genotype_df["PRS"] = 0.0

# Calculate PR
for _, row in df.iterrows():
    rsid = row["rsID"]
    beta = row["beta_final"]
    genotype_df["PRS"] += (
        genotype_df[rsid] * beta
    )

# Check output
print("\nSimulated genotype matrix with PRS:\n")
print(genotype_df.head())

# Print PRS summary statistics
print("\nPRS summary:\n")
print(genotype_df["PRS"].describe())

# # Save output
# genotype_df.to_csv(
#     "data/simulated_genotypes.csv",
#     index=False
# )

# print("\nSaved: data/simulated_genotypes.csv")
# -----------------------------
# PRS histogram visualization
# -----------------------------

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))

plt.hist(
    genotype_df["PRS"],
    bins=30
)

plt.xlabel("Polygenic Risk Score (PRS)")
plt.ylabel("Number of Individuals")
plt.title("Distribution of CAD Polygenic Risk Scores")

plt.grid(True)

# Save figure
plt.savefig(
    "pic/prs_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nSaved: pic/prs_distribution.png")

# -----------------------------
# Top / bottom percentile analysis
# -----------------------------

# Calculate percentile thresholds
top_5_threshold = genotype_df["PRS"].quantile(0.95)
bottom_5_threshold = genotype_df["PRS"].quantile(0.05)

# Extract high-risk individuals
top_5 = genotype_df[
    genotype_df["PRS"] >= top_5_threshold
]

# Extract low-risk individuals
bottom_5 = genotype_df[
    genotype_df["PRS"] <= bottom_5_threshold
]

print("\nTop 5% PRS individuals:\n")
print(
    top_5[
        ["person_id", "PRS"]
    ].head()
)

print("\nBottom 5% PRS individuals:\n")
print(
    bottom_5[
        ["person_id", "PRS"]
    ].head()
)

print("\nTop 5% threshold:", top_5_threshold)
print("Bottom 5% threshold:", bottom_5_threshold)

# -----------------------------
# SNP contribution analysis
# -----------------------------

contrib_list = []

for _, row in df.iterrows():

    rsid = row["rsID"]
    beta = row["beta_final"]
    raf = row["RAF"]

    # Approximate contribution to PRS variance
    contribution = (
        2 * raf * (1 - raf) * (beta ** 2)
    )

    contrib_list.append({
        "rsID": rsid,
        "gene": row["gene"],
        "RAF": raf,
        "beta": beta,
        "contribution": contribution
    })

# Create contribution dataframe
contrib_df = pd.DataFrame(contrib_list)

# Sort SNPs by contribution
contrib_df = contrib_df.sort_values(
    by="contribution",
    ascending=False
)

print("\nTop SNP contributions:\n")
print(contrib_df.head(10))

# -----------------------------
# SNP contribution visualization
# -----------------------------

plt.figure(figsize=(12, 6))

plt.bar(
    contrib_df["rsID"].head(10),
    contrib_df["contribution"].head(10)
)

plt.xlabel("SNP")
plt.ylabel("Contribution to PRS Variance")
plt.title("Top SNP Contributions to PRS")

plt.xticks(rotation=45)

plt.grid(True)

# Save figure
plt.savefig(
    "pic/snp_contributions.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nSaved: pic/snp_contributions.png")