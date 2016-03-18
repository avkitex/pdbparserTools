import string

from chemspipy import ChemSpider

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

def leftOnlyLettersDigits(sn):
    res = ''
    for let in sn:
        if let in string.ascii_letters + string.digits + ['_', '|', '-']:
            res += let
    return let

def getChemspiderCompounds(token, list, pref, delim = '_', longNames = True, onlyLettersDigits = False):
    cs = ChemSpider(token)
    names = []
    molecules = []
    for chsId in list:
        comp = cs.get_compound(chsId)
        name = pref + delim + str(chsId)
        if longNames:
            name += delim
            sn = comp.common_name.encode('ascii','ignore')
            if onlyLettersDigits:
                sn = leftOnlyLettersDigits(sn)
            name += sn
        #.replace('(', '_').replace(')', '_').replace('[', '_').replace(']', '_').replace(',', '_').replace(' ', '_').replace(';', '_')[:25]
        print(name)
        smiles=comp.smiles.encode('ascii','ignore')
        molecules.append(Chem.MolFromSmiles(smiles))
        names.append(name)
    return molecules, names

def getTrainingCompounds(inhibitorsChsIds=[], notInhibitorsChsIds=[], longNames = True, onlyLettersDigits = False, token = '2228d430-a955-416b-b920-14547d28df9e'):

    moleculesY, namesY = getChemspiderCompounds(token, inhibitorsChsIds, 'yes', '_', longNames, onlyLettersDigits)

    moleculesN, namesN = getChemspiderCompounds(token, notInhibitorsChsIds, 'not', '_', longNames, onlyLettersDigits)
    return moleculesY + moleculesN, namesY + namesN


# def similarityMatrixToTree(inFile):
#     handle = open(inFile)
#     amount = int(handle.next().strip())
#     names = handle.next().strip().split()
#     simil = []
#     for line in handle:
#         simil.append([float(x) for x in line.strip().split()])
#     dm =_DistanceMatrix(names, simil)
#     constructor = DistanceTreeConstructor()
#     tree = constructor.upgma(dm)
#     draw(tree)
#     return tree
# def operate(inFile, outFile):
#     molecules, names = readMoleculesRd(inFile)
#     tree = clusterMolecules(molecules, names)
#     #draw(tree)
#     #Phylo.write(tree, outFile, 'newick')
# def operate_tr(outFile):
#     molecules, names = getTrainingCompounds()
#     tree = clusterMolecules(molecules, names)
#     #draw(tree)
#     #Phylo.write(tree, outFile, 'newick')
#     Phylo.draw_ascii(tree)
#
#
#
# def ClusterFps(fps,cutoff=0.2):
#     from rdkit.ML.Cluster import Butina
#
#     # first generate the distance matrix:
#     dists = []
#     nfps = len(fps)
#     for i in range(1,nfps):
#         sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
#         dists.extend([1-x for x in sims])
#
#     # now cluster the data:
#     cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
#     return cs
# #operate_tr('out.dnd')
