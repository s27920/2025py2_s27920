from Bio import Entrez as E, SeqIO as S
from io import StringIO as T
import csv
import matplotlib.pyplot as P
G=input
N=None
O=open
R=range
L=len
X=Exception
I=int
Q=exit
PR=print
TX="taxid"
AC="accessions"
SL="Sequence Length"
NC="nucleotide"
LN=" len [optional]: "
if __name__=="__main__":
    E.email=G("email:")
    E.api_key=G("API key:")
    E.tool='BioScriptEx10'
    t=G(f"Enter {TX}:").strip()
    try:
        mni=G(f"Min{LN}").strip()
        mxi=G(f"Max{LN}").strip()
        mnl=I(mni) if mni else N
        mxl=I(mxi) if mxi else N
        if mnl and mxl and mnl>mxl:
            mnl,mxl=mxl,mnl
    except:
        PR("Invalid filter");Q()
    PR(f"Searching for records with {TX}: {t}")
    on=E.read(E.efetch(db="taxonomy",id=t,retmode="xml"))[0]["ScientificName"]
    PR(f"Organism: {on} ({TX}: {t})")
    sr=E.read(E.esearch(db=NC,term=f"txid{t}[Organism] AND {mnl if mnl is not N else 0}:{mxl if mxl is not N else 1e12}[SLEN]",usehistory="y"))
    c=I(sr["Count"])
    if c==0:
        PR("No records found matching the criteria")
        Q()
    PR(f"Found {c} records")
    dt,bs = [],500
    for s in R(0,c,bs):
        dt.append(E.efetch(db=NC, rettype="gb", retmode="text", retstart=s, retmax=min(bs,c-s), webenv=sr["WebEnv"], query_key=sr["QueryKey"]).read())
    rt = ''.join(dt)
    rc = list(S.parse(T(rt), "genbank"))
    with O(f"{TX}_{t}_records.gb", "w") as f:
        f.write(rt)
    with O(f"{TX}_{t}_report.csv", "w", newline='') as f:
        csv.writer(f).writerows([[r.annotations.get(AC, [""])[0], L(r.seq), r.description] for r in rc])
    sr=sorted(rc,key=lambda x: L(x.seq),reverse=True)
    ac=[rec.annotations.get(AC,[''])[0] for rec in sr]
    P.figure(figsize=(12, 6))
    P.plot(R(L(ac)),[L(rec.seq) for rec in sr],marker='o',linestyle='-',color='b')
    P.xticks(R(L(ac)),ac,rotation=45,ha='right')
    P.ylabel(SL)
    P.title(f'{SL} Dist')
    P.savefig(f"{TX}_{t}_lengths.png")
    P.close()