from genetools.gene_data import GeneData
import pandas as pd
from datetime import datetime
from pathlib import Path

# genes to compare exons of
query = GeneData("CFTR", species = "Human")
target = GeneData("CFTR", species = "Mouse")

# list to hold results
results = []

# get exon lists
qExons = query.exons
tExons = target.exons

# go through each exon pair, align with respect to AA and cDNA, save in results dict
for idx,_ in enumerate(tExons):

    # get this set of exons
    q = qExons[idx]
    t = tExons[idx]

    # cDNA sequences for exons
    cdna1 = query.get_exon_cdna(q)
    cdna2 = target.get_exon_cdna(t)

    # AA sequences for exons
    aa1 = q['aa_seq']
    aa2 = t['aa_seq']

    # generate alignments
    a = GeneData.align_aa(aa1,aa2)
    c = GeneData.align_cdna(cdna1,cdna2)

    # get alingments as human readable strings
    strA = str(a) if len(a) > 0 else 'ssa Error'
    strC = str(c) if len(c) > 0 else 'ssa Error'

    # splice site analysis for acceptor
    og3,og5 = target.splice_site_analysis(idx)

    # swich donor with acceptor and analyze new splice site
    check,new = target.swap_exons(idx,cdna1)
    new3,new5 = GeneData.swapped_ssa(check,new)

    # calculate differnces between new and old exons
    try:
        old3 = float(og3)
    except ValueError:
        old3 = 'N/A'
    try:
        old5 = float(og5)
    except ValueError:
        old5 = 'N/A'
    try:
        swap3 = float(new3)
    except ValueError:
        swap3 = 'N/A'
    try:
        swap5 = float(new5)
    except ValueError:
        swap5 = 'N/A'
    try:    
        diff3 = abs(swap3-old3)
    except TypeError:
        diff3 = 'N/A'
    try:
        diff5 = abs(swap5-old5)
    except TypeError:
        diff5 = 'N/A'

    # tell if new value is better or worse than old
    # 3'ss
    if old3 == 'N/A':
        val3 = 'N/A'
    elif old3 > swap3:
        val3 = 'Better'
    elif old3 == swap3:
        val3 = 'Same'
    else:
        val3 = 'Worse'
    #5'ss
    if old5 == 'N/A':
        val5 = 'N/A'
    if old5 > swap5:
        val5 = 'Better'
    elif old5 == swap5:
        val5 = 'Same'
    else:
        val5 = 'Worse'

    # save data in a dictionary
    data = {
        'exon' : t['rank'],
        'alignAA' : strA,
        'aligncDNA' : strC,
        'aaScore' : a.score,
        'cDNAScore' : c.score,
        'og3ssa' : old3,
        'og5ssa' : old5,
        'new3ssa' : swap3,
        'new5ssa' : swap5,
        'diff3' : diff3,
        'diff3_significance' : val3,
        'diff5' : diff5,
        'diff5_significance' : val5
    }
    
    # append dictionary to results list
    results.append(data)
    print(f"Exon {t['rank']} saved")

# convert dict to pandas dataframe
df = pd.DataFrame(results)

# define outputs folder (relative to this script)
outputs_dir = Path(__file__).parent / "outputs"
outputs_dir.mkdir(parents=True, exist_ok=True)  # create folder if missing

# get timestamp for analysis
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# build full file path
output_file = outputs_dir / f"{timestamp}_alignment_results.xlsx"

# write results to excel
df.to_excel(output_file, index=False)