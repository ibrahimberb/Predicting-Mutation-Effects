def unzip_res_range(res_range):
    '''Converts ranges in the form: [2-210] or [3-45,47A,47B,51-67] into lists of strings including all numbers in these ranges in order'''
    import re
    res_ranges = res_range.strip()[1:-1].split(',')
    index_list = []
    for r in res_ranges:
        if re.match('.+-.+', r):
            a, b = r.split('-')
            index_list += [str(n) for n in range(int(a), int(b)+1)]
        else:
            index_list.append(r)
    return index_list

def setValues(nested_list):
    '''Converts nested lists to a list of sorted unique values'''
    l = []
    for li in nested_list: l+=li
    return sorted(list(set(l)))




# Define a set of high-quality binary human protein-protein interactions (the HINT database)
hint = set([tuple(sorted(line.strip("\n").split("\t")[:2])) for line in open("HomoSapiens_binary_hq.txt").readlines()[1:]])

# Map interface predictions from Interactome INSIDER onto HINT protein interactions
# on a per-interaction level
prot_ECLAIRiresALL_perint={}
for line in open("H_sapiens_interfacesALL.txt"):
    p1,p2,src,ires1,ires2=line.strip().split("\t")
    if tuple(sorted([p1,p2])) not in hint: continue
    if p1 not in prot_ECLAIRiresALL_perint: prot_ECLAIRiresALL_perint[p1] = {}
    if p2 not in prot_ECLAIRiresALL_perint: prot_ECLAIRiresALL_perint[p2] = {}
    prot_ECLAIRiresALL_perint[p1][p2] = unzip_res_range(ires1)
    prot_ECLAIRiresALL_perint[p2][p1] = unzip_res_range(ires2)

    
# Dissect interface predictions based on prediction source and confidence 
# PDB: predictions based on protein PDB structures; 
# PI: predictions based on protein PDB structures or homology models
# HQ: high quality predictions
# PO: prediction only, excluding PI
eclair_model = open("interaction_model_index.PI.tsv").readlines()
prot_ECLAIRresPI_perint={}
prot_ECLAIRresPDB_perint = {}
prot_ECLAIRiresHQ_perint={}
prot_ECLAIRiresHQPO_perint={} 
for line in open("H_sapiens_interfacesHQ.txt"):
    p1,p2,src,ires1,ires2=line.strip().split("\t")
    if tuple(sorted([p1,p2])) not in hint: continue
    ires1 = unzip_res_range(ires1)
    ires2 = unzip_res_range(ires2)
    if p1 not in prot_ECLAIRiresHQ_perint: prot_ECLAIRiresHQ_perint[p1] = {}
    if p2 not in prot_ECLAIRiresHQ_perint: prot_ECLAIRiresHQ_perint[p2] = {}
    prot_ECLAIRiresHQ_perint[p1][p2] = ires1
    prot_ECLAIRiresHQ_perint[p2][p1] = ires2
    if src == "PDB":
        if p1 not in prot_ECLAIRiresPDB_perint: prot_ECLAIRiresPDB_perint[p1] = {}
        if p2 not in prot_ECLAIRiresPDB_perint: prot_ECLAIRiresPDB_perint[p2] = {}
        prot_ECLAIRiresPDB_perint[p1][p2] = ires1
        prot_ECLAIRiresPDB_perint[p2][p1] = ires2
    if src in ["PDB","I3D"]:
        if p1 not in prot_ECLAIRiresPI_perint: prot_ECLAIRiresPI_perint[p1] = {}
        if p2 not in prot_ECLAIRiresPI_perint: prot_ECLAIRiresPI_perint[p2] = {}
        prot_ECLAIRiresPI_perint[p1][p2] = ires1
        prot_ECLAIRiresPI_perint[p2][p1] = ires2
    else:
        if p1 not in prot_ECLAIRiresHQPO_perint: prot_ECLAIRiresHQPO_perint[p1] = {}
        if p2 not in prot_ECLAIRiresHQPO_perint: prot_ECLAIRiresHQPO_perint[p2] = {}
        prot_ECLAIRiresHQPO_perint[p1][p2] = ires1
        prot_ECLAIRiresHQPO_perint[p2][p1] = ires2

        
# Collapse per-interaction level interfaces to per-protein level
prot_ECLAIRiresALL_perprot={}
prots = prot_ECLAIRiresALL_perint.keys()
for prot in prots:
    prot_ECLAIRiresALL_perprot[prot] = [i for i in setValues(prot_ECLAIRiresALL_perint[prot].values()) if i !=""]
prot_ECLAIRiresHQ_perprot={}
prots = prot_ECLAIRiresHQ_perint.keys()
for prot in prots:
    prot_ECLAIRiresHQ_perprot[prot] = [i for i in setValues(prot_ECLAIRiresHQ_perint[prot].values()) if i !=""]
prot_ECLAIRiresPI_perprot={}
prots = prot_ECLAIRiresPI_perint.keys()
for prot in prots:
    prot_ECLAIRiresPI_perprot[prot] = [i for i in setValues(prot_ECLAIRiresPI_perint[prot].values()) if i !=""]
prot_ECLAIRiresPDB_perprot={}
prots = prot_ECLAIRiresPDB_perint.keys()
for prot in prots:
    prot_ECLAIRiresPDB_perprot[prot] = [i for i in setValues(prot_ECLAIRiresPDB_perint[prot].values()) if i !=""]


# Functions to obtain interaction partners associated (or not) with a particular protein residue
# Input: protein_id, protein_residue_number, confidence_level_of_interface_prediction 
# Returns: a list of protein interaction partners that are associated with the query protein residue and a list of partners that are not associated
def getIresPartners(prot, res, eclair_level):
    if eclair_level == "pdb": eclair = prot_ECLAIRiresPDB_perint
    elif eclair_level == "pi": eclair = prot_ECLAIRiresPI_perint
    elif eclair_level == "hq": eclair = prot_ECLAIRiresHQ_perint
    elif eclair_level == "all": eclair = prot_ECLAIRiresALL_perint
    if prot not in eclair: return []
    else: return [interactor for interactor in eclair[prot] if res in eclair[prot][interactor]]
    
def getModeledPartners(prot, res, eclair_level):
    if eclair_level == "pdb": 
        eclair = prot_ECLAIRresPDB_perint
        if prot not in eclair: return []
        return [interactor for interactor in eclair[prot] if res in eclair[prot][interactor]]
    elif eclair_level == "pi": 
        eclair = prot_ECLAIRresPI_perint
        if prot not in eclair: return []
        return [interactor for interactor in eclair[prot] if res in eclair[prot][interactor]]
    elif eclair_level == "hq": 
        eclair = prot_ECLAIRiresHQPO_perint
        po = list(eclair[prot].keys()) if prot in eclair else []
        return po+getModeledPartners(prot, res, "pi")
    elif eclair_level == "all": 
        eclair = prot_ECLAIRiresALLPO_perint
        po = list(eclair[prot].keys()) if prot in eclair else []
        return po+getModeledPartners(prot, res, "pi")

def getIresPartners2(prot, res, eclair_level):
    ires_par = getIresPartners(prot, res, eclair_level)
    all_par = getModeledPartners(prot, res, eclair_level)
    nires_par = [i for i in all_par if i not in set(ires_par)]
    return ires_par,nires_par


