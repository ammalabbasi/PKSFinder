#!/usr/bin/env python3

import pandas as pd
import sys
import os

# ----------------- Input arguments ----------------- #
counts_files = sys.argv[1:-1]   # List of counts.txt files
output_file = sys.argv[-1]      # Output combined table

# ----------------- Collect and Merge ----------------- #
merged_df = None

for file_path in counts_files:
    sample_name = os.path.basename(file_path).replace(".counts.txt", "")
    
    # Skip header lines starting with #
    df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
    
    # Extract relevant columns
    if merged_df is None:
        # First file: keep all annotation columns
        merged_df = df.iloc[:, :6].copy()  # Geneid, Chr, Start, End, Strand, Length
        merged_df[sample_name] = df.iloc[:, 6]
    else:
        merged_df[sample_name] = df.iloc[:, 6]

# ----------------- Save Output ----------------- #
merged_df.to_csv(output_file, sep="\t", index=False)
print(f"Merged counts table written to: {output_file}")
