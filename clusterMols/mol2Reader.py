#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

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

        simil.append(DataStructs.BulkTanimotoSimilarity(vects[mol1], vects[:mol1]))
    return simil
    #dm=_DistanceMatrix(['mol' + str(x) for x in range(99)], mols[1:])
    # tree = constructor.nj(dm)
    #constructor = DistanceTreeConstructor()
    #from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    #from Bio.Phylo.TreeConstruction import _DistanceMatrix

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
