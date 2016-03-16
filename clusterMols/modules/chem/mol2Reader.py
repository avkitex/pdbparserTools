#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install python-rdkit librdkit1 rdkit-data libfreetype6-dev
# sudo pip install biopython chemspipy

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo import draw

from ..common.f import iterMol2

def readMoleculesRd(file):
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
    inhibitorsChsIds=[331, 1353, 1906, 2000, 2066, 2121, 2157, 2562, 2925, 3065, 3097, 3192, 3225, 3254, 3544, 3584, 3693, 3694, 3897, 3904, 4339, 4393, 4480, 4617, 4911, 5308, 8711, 133236, 28714, 29843, 32983, 33051, 34911, 65084, 110209, 137720, 147913, 171462, 388601, 392692, 392977, 393569, 394812, 394942, 581047, 8082544, 1265915, 2871682, 3512137, 4444105, 4444112, 4445953, 4510145, 5365258, 8081394, 8098374, 8443182, 8530675, 9569918, 10442653, 10442740, 12951165, 20572535, 21106585, 23122887, 23122889, 23122978, 25057753]
    notInhibitorsChsIds=[650, 682, 733, 864, 937, 971, 1116, 1512, 2971, 3350, 5611, 5653, 5764, 5768, 6170, 6257, 6312, 7742, 10610, 73505, 83361, 96749, 111188, 120261, 216840, 391555, 392800, 4576521, 8257952, 20572534, 23122865, 23122927, 23123076, 112728, 10442445, 1710, 5304]
    names = []
    molecules = []
    for chsId in inhibitorsChsIds:
        comp = cs.get_compound(chsId)
        name = 'yes_' + str(chsId) + '_' + comp.common_name.encode('ascii','ignore')
        #.replace('(', '_').replace(')', '_').replace('[', '_').replace(']', '_').replace(',', '_').replace(' ', '_').replace(';', '_')[:25]
        print(name)
        smiles=comp.smiles.encode('ascii','ignore')
        molecules.append(Chem.MolFromSmiles(smiles))
        names.append(name)
    for chsId in notInhibitorsChsIds:
        comp = cs.get_compound(chsId)
        name = 'not_' + str(chsId) + '_' + comp.common_name.encode('ascii','ignore')
        #.replace('(', '_').replace(')', '_').replace('[', '_').replace(']', '_').replace(',', '_').replace(' ', '_').replace(';', '_')[:25]
        print(name)
        smiles=comp.smiles.encode('ascii','ignore')
        molecules.append(Chem.MolFromSmiles(smiles))
        names.append(name)
    return molecules, names


def similarityMatrixToTree(inFile):
    handle = open(inFile)
    amount = int(handle.next().strip())
    names = handle.next().strip().split()
    simil = []
    for line in handle:
        simil.append([float(x) for x in line.strip().split()])
    dm =_DistanceMatrix(names, simil)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    draw(tree)
    return tree
def operate(inFile, outFile):
    molecules, names = readMoleculesRd(inFile)
    tree = clusterMolecules(molecules, names)
    #draw(tree)
    #Phylo.write(tree, outFile, 'newick')
def operate_tr(outFile):
    molecules, names = getTrainingCompounds()
    tree = clusterMolecules(molecules, names)
    #draw(tree)
    #Phylo.write(tree, outFile, 'newick')
    Phylo.draw_ascii(tree)



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
#operate_tr('out.dnd')