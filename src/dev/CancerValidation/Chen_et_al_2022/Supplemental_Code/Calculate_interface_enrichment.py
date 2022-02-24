# Function to calculate the enrichment of given mutations on protein interaction interfaces
# Input: 
# 1) a list of protein residue numbers where mutations happened
# 2) confidence level of interface prediction 
# Returns:
# numbers and percentages of mutated residues on protein interfaces,
# numbers and percentages of all residues on protein interfaces (background rate),
# and enrichment of mutations on protein interaction interfaces (fold change and p-value by a binomial exact test)

def iresEnrichment(list_res, eclair_level):
    x1,n1 = 0,0
    x2,n2 = 0,0
    if eclair_level == "pdb": eclair = prot_ECLAIRiresPDB_perprot
    elif eclair_level == "pi": eclair = prot_ECLAIRiresPI_perprot
    elif eclair_level == "hq": eclair = prot_ECLAIRiresHQ_perprot
    elif eclair_level == "all": eclair = prot_ECLAIRiresALL_perprot
    for prot,res in list_res:
        if len(eclair[prot]) == 0: continue
        n1+=1
        if res in eclair[prot]: x1+=1
        if prot not in prots: 
            x2 += len(eclair[prot])
            n2 += uniprot2len[prot]
        prots.add(prot)
        
    r1 = float(x1)/n1
    r2 = float(x2)/n2
    fc = r1/r2
    return [x1,n1,r1],[x2,n2,r2],fc,scipy.stats.binom_test(x1,n1, p=float(x2)/n2)
    