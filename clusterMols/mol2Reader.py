#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install rdkit rdkit-devel libfreetype6-dev
# sudo pip install biopython

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo import draw

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
def readMolecules(file):
    molecules=[]
    names=[]
    for lines in iterMol2(file):
        molecules.append(Chem.MolFromMol2Block(''.join(lines)))
        names.append(lines[1])
    print(len(molecules))
    return molecules, names

def clusterMolecules(molecules, names):

    simil = []
    vects = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in molecules]
    for mol1 in range(len(vects)):
        simil.append([1-x for x in DataStructs.BulkTanimotoSimilarity(vects[mol1], vects[:mol1+1])])
    dm =_DistanceMatrix(names, simil)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return tree
#Bio.Phylo.draw(tree)

def getTrainingCompounds():
    from chemspipy import ChemSpider
    cs = ChemSpider('2228d430-a955-416b-b920-14547d28df9e')
    inhibitorsChsIds=[331, 650, 1353, 1906, 2000, 2066, 2121, 2157, 2483, 2562, 2925, 3065, 3097, 3192, 3225, 3254, 3277, 3544, 3584, 3693, 3694, 3897, 3904, 4339, 4393, 4480, 4617, 5308, 8711, 21386, 28714, 29843, 32983, 33051, 34911, 65084, 137720, 147913, 171462, 388601, 392692, 392977, 393569, 394812, 394942, 581047, 1265915, 2871682, 3512137, 4444105, 4444112, 4445953, 4510145, 5365258, 8081394, 8098374, 8443182, 8530675, 9569918, 10442653, 10442740, 12951165, 20572535, 21106585, 23122887, 23122889, 23122978, 25057753]
    notInhibitorsChsIds=[682, 733, 864, 937, 971, 1116, 1512, 2562, 2925, 2971, 3065, 3350, 4339, 4911, 5611, 5653, 5764, 5768, 6170, 6257, 6312, 7742, 10610, 73505, 83361, 96749, 110209, 111188, 120261, 133236, 216840, 391555, 392800, 4576521, 8082544, 8257952, 10442740, 20572534, 23122865, 23122927, 23123076, 112728, 10442445, 1710, 5304]
    names = []
    molecules = []
    for chsId in inhibitorsChsIds:
        name = 'yes_' + str(chsId)
        print(name)
        comp = cs.get_compound(chsId)
        smiles=comp.smiles.encode('ascii','ignore')
        molecules.append(Chem.MolFromSmiles(smiles))
        names.append(name)
    for chsId in notInhibitorsChsIds:
        name = 'not_' + str(chsId)
        print(name)
        comp = cs.get_compound(chsId)
        smiles=comp.smiles.encode('ascii','ignore')
        molecules.append(Chem.MolFromSmiles(smiles))
        names.append(name)
    return molecules, names




def operate(inFile, outFile):
    molecules, names = readMolecules(inFile)
    tree = clusterMolecules(molecules, names)
    #draw(tree)
    Phylo.write(tree, outFile, 'newick')
def operate_tr(outFile):
    molecules, names = getTrainingCompounds()
    tree = clusterMolecules(molecules, names)
    #draw(tree)
    Phylo.write(tree, outFile, 'newick')



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
