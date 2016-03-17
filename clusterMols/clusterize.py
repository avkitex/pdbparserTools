#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install python-rdkit librdkit1 rdkit-data libfreetype6-dev python-networkx python-PyGraphviz
# sudo pip install biopython chemspipy pylab



from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from rdkit.DataStructs.cDataStructs import SparseBitVect
from Bio import Phylo
from Bio.Phylo import draw
import pylab

from modules.positional.kitsite import *
from modules.chem.mol2Reader import *

def getFormatedTime():
	return str(datetime.utcnow().strftime("[%H:%M:%S]"))

def appendBitVectors(bitVector1, bitVector2):
    res = SparseBitVect(bitVector1.GetNumBits() + bitVector2.GetNumBits())
    for i in range(bitVector1.GetNumBits()):
        res[i] = bitVector1[i]
    for i in range(bitVector2.GetNumBits()):
        res[i + bitVector1.GetNumBits()] = bitVector2[i]
    return res
def appendChemBoxBitVectors(chemVectors, chemNames, boxVectors, boxNames):
    newNames = []
    newVectors = []
    for cname in range(len(chemNames)):
        shortName = chemNames[cname][4:chemNames[cname][4:].find('_') + 4]
        #print(shortName)
        for bname in range(len(boxNames)):
            if shortName in boxNames[bname]:
                #print(boxNames[bname], chemNames[cname])
                newNames.append(chemNames[cname])
                newVectors.append(appendBitVectors(chemVectors[cname], boxVectors[bname]))
                break
    return newVectors, newNames
def clusterBitVectors(vectors, names):
	simil = []
	for mol1 in range(len(bothVectors)):
	    simil.append([1-x for x in DataStructs.BulkTanimotoSimilarity(bothVectors[mol1], bothVectors[:mol1+1])])
	dm =_DistanceMatrix(bothNames, simil)
	constructor = DistanceTreeConstructor()
	tree = constructor.upgma(dm)
	#Phylo.draw_ascii(tree)
	Phylo.draw_graphviz(tree)
	pylab.show()
	return tree
print(getFormatedTime() + " Downloading ideal training molecules for chemical clustering")

chemMolecules, chemNames = getTrainingCompounds()
print(getFormatedTime() + " Obtaining bit vectors representation")
chemVectors = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in chemMolecules]



############################ PARAMS #####################################
protinFile='3NTB_D.mol2'

asCenterX = -26.9
asCenterY = 21.9
asCenterZ = -76.9

bondLenBoxExtend = 0
bondLenClustering = 3.5

gridSizeX = 30
gridSizeY = 30
gridSizeZ = 30

stepSize = 0.5
minCavSize = 5
topAtomsPersent = 20

cornerPenalty = 10
############################ PARAMS #####################################

############################################# params adopting ############################
gridSize = point3D(gridSizeX + bondLenBoxExtend, gridSizeY + bondLenBoxExtend, gridSizeZ + bondLenBoxExtend)
centerCoords = point3D(asCenterX, asCenterY, asCenterZ)
box = boxParams(centerCoords, gridSize)
############################################# params adopting ############################

print(getFormatedTime() + " Getting protein grid box and filtering active site atoms")

filteredBox=filterBoxAtoms('3NTB_D.mol2', box, stepSize, topAtomsPersent)
print(getFormatedTime() + " Getting bit vector representation")
boxVectors, boxNames = getMoleculesContactsAsBitVect('../../training/training_bp.mol2', filteredBox)

print(getFormatedTime() + " Composing bit vectors")
bothVectors, bothNames = appendChemBoxBitVectors(chemVectors, chemNames, boxVectors, boxNames)

print(getFormatedTime() + " Final tree")
#clusterBitVectors(bothVectors, bothNames)
