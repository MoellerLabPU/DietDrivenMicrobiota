import pandas as pd
import glob
import os
import sys

# Check if a folder path was provided
if len(sys.argv) < 2:
    print("Usage: python script.py /path/to/ratio_files")
    sys.exit(1)

input_folder = sys.argv[1]

# Make sure the folder exists
if not os.path.isdir(input_folder):
    print(f"Error: folder '{input_folder}' does not exist.")
    sys.exit(1)

# Create output folder inside the input folder
output_folder = os.path.join(input_folder, "filtered_output_new")
os.makedirs(output_folder, exist_ok=True)

# Loop through all .ratio.txt files in the folder
for file in glob.glob(os.path.join(input_folder, "*.ratio.txt")):
    df = pd.read_csv(file, sep="\t")
    
    # Apply filtering conditions
    filtered = df[
       (
            (df["Pe_ratio"] >= 0.01) &
            (df["Pe_R"] > 5) &
            (df["Span_R"] > 3)
       )       
       |
       (
            (df["Pe_F"] > 5) &
            (df["Span_F"] > 3)
       )
    ]
    
    # Save only if there are matching rows
    if not filtered.empty:
        out_file = os.path.join(output_folder, os.path.basename(file).replace(".ratio.txt", "_filtered.txt"))
        filtered.to_csv(out_file, sep="\t", index=False)
        print(f"Filtered rows saved to {out_file}")
    else:
        print(f"No matching rows in {file}")
