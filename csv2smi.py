import sys
import pandas as pd
if len(sys.argv) == 1:
    print("usage: csv2smi.py input.csv [output.smi]")
    exit()
file_out: str = "output.smi"
if len(sys.argv) >= 2:
    file_in = sys.argv[1]
if len(sys.argv) >= 3:
    file_out = sys.argv[2]

f = pd.read_csv(file_in)
s = [s for s in f["smiles"]]
with open(file_out, "w") as f:
    for ss in s:
        f.write(f"{ss:s}\n")
