# process_data.py

import pandas as pd
import numpy as np

df = pd.read_csv("data/SNP_Dataset.csv")

print(df.head())

# # -----------------------------
# # OR -> beta
# # -----------------------------
df["pvalue"] = (df["pvalue"].astype(float))
df["beta_final"] = np.nan

# # 1. use beta if it exists
beta_mask = df["beta"].notna()

df.loc[beta_mask, "beta_final"] = df.loc[beta_mask, "beta"]

# # 2. if OR exists, use log(OR)
or_mask = df["OR"].notna()

df.loc[or_mask, "beta_final"] = np.log(
    df.loc[or_mask, "OR"]
)

# Check
print("\nFinal PRS weights:\n")

print(
    df[
        [
            "rsID",
            "risk_allele",
            "OR",
            "beta",
            "beta_final",
            "RAF"
        ]
    ]
)
# print(df.dtypes)
# # save
df.to_csv("data/cad_snps_processed.csv", index=False)
print("\nSaved: data/cad_snps_processed.csv")