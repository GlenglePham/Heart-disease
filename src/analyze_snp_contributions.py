import os

import pandas as pd
import matplotlib.pyplot as plt

# Create output directory
os.makedirs("pic", exist_ok=True)

# Load processed SNP data
df = pd.read_csv("data/cad_snps_processed.csv")

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

# Print results
print("\nTop SNP contributions:\n")
print(contrib_df.head(10))

# -----------------------------
# Visualization
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

plt.savefig(
    "pic/snp_contributions.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nSaved: pic/snp_contributions.png")