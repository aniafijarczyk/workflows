import pandas as pd
import glob

coords = glob.glob("nucmer_*_qr.coords")
fname = coords[1]

D = []
for fname in coords:
    df = pd.read_csv(fname, sep="\t", skiprows = [0,1,2,3], header=None, names=["S1","E1","S2","E2","Len1","Len2","%ID",
                                                                                       "LenR","LenQ","CovR","CovQ","Ref","ContigName"])
    name = fname.split("/")[-1].replace("nucmer_","").replace("_qr.coords","")
    #print(name)
    if df.empty == True:
        print("empty")
    else:
        df["chrom"] = df["ContigName"]
        dh = df.groupby(["Ref"]).agg({'CovR':'sum','%ID':'mean'}).reset_index().rename(columns = {'CovR':'frac_mtdna','%ID':'%ID_mtdna'})
        dh['species'] = name
        D.append(dh)
dD = pd.concat(D)
dD.head()
dD.to_csv('getMtdnaInfo.tab', sep="\t", index=False, header=True, na_rep="NA")
