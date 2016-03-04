#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install rdkit rdkit-devel libfreetype6-dev
# sudo pip install biopython

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def iterMol2(file):
    handle = open(file, 'r')
    if handle:
        curMol2 = []
        for line in handle:
            if '@<TRIPOS>MOLECULE' in line:
                if len(curMol2):
                    yield curMol2
                    curMol2 = []
            curMol2.append(line)
        if len(curMol2):
            yield curMol2
        handle.close()

def clusterMolecules(file):
    molecules=[]
    simil = []
    for file in iterMol2(file):
        molecules.append(Chem.MolFromMol2Block(''.join(file)))
    print(len(molecules))
    vects = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in molecules]
    for mol1 in range(len(vects)):
        simil.append(DataStructs.BulkTanimotoSimilarity(vects[mol1], vects[:mol1+1]))
    dm =_DistanceMatrix(['mol' + str(x+1) for x in range(len(simil))], simil)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return tree
#Bio.Phylo.draw(tree)
def ClusterFps(fps,cutoff=0.2):
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs
