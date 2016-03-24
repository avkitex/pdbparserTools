#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install -y python-rdkit librdkit1 rdkit-data libfreetype6-dev python-networkx python-PyGraphviz
# sudo pip install biopython chemspipy

#deprecated pylab

# sudo apt-get install -y emboss embassy-phylip
# fneighbor -matrixtype l -datafile bigDM.dm -treetype u -outfile out

from __future__ import print_function
import os.path, argparse, gc
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import SparseBitVect

from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo import draw

#import pylab

from modules.positional.kitsite import *
from modules.chem.mol2Reader import *

parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-tf', '--trainingMol2', metavar='GlobConfig', type=str, help='Full path to multiMol2 training compounds docked file.', required=True)
parser.add_argument('-dm', '--distanceMatrixOutput', metavar='GlobConfig', type=str, help='Full path to distance matrix output file.', required=True)
#parser.add_argument('-inhl', '--inhibitorsList', metavar='GlobConfig', type=str, help='inhibitors list.', required=True)
#parser.add_argument('-ninhl', '--notInhibitorsList', metavar='GlobConfig', type=str, help='Not inhibitors list.', required=True)
args = parser.parse_args()


############################ PARAMS #####################################
if not os.path.isfile(args.trainingMol2):
	print('Docked training compounds mol2 file is not exists')
#if not os.path.isfile(args.inhibitorsList):
#	print('Docked training compounds mol2 file is not exists')
#if not os.path.isfile(args.notInhibitorsList):
#	print('Docked training compounds mol2 file is not exists')

inhibitorsChsIds=[331, 1353, 1906, 2000, 2066, 2121, 2157, 2562, 2925, 3065, 3097, 3192, 3225,\
3254, 3544, 3584, 3693, 3694, 3897, 3904, 4339, 4393, 4480, 4617, 4911, 5304, 5308, 8711, 133236,\
28714, 29843, 32983, 33051, 34911, 65084, 110209, 137720, 147913, 171462, 388601, 392692, 392977,\
393569, 394812, 394942, 581047, 8082544, 1265915, 2871682, 3512137, 4444105, 4444112, 4445953, \
4510145, 5365258, 8081394, 8098374, 8443182, 8530675, 9569918, 10442653, 10442740, 20572535, \
21106585, 23122887, 23122889, 23122978, 25057753] # 12951165 no results case of bor
notInhibitorsChsIds=[650, 682, 733, 864, 971, 1116, 1512, 1710, 2971, 3350, 5611, 5653, 5764, 5768, 6170, 6257, 6312, 7742, 10610, 73505, 83361, 96749, 111188, 120261, 216840, 391555,\
392800, 4576521, 8257952, 20572534, 23122865, 23122927, 23123076, 112728, 10442445]


def getFormatedTime():
	return str(datetime.utcnow().strftime("[%H:%M:%S]"))

def appendBitVectors(bitVector1, bitVector2):
	res = SparseBitVect(bitVector1.GetNumBits() + bitVector2.GetNumBits())
	for i in range(bitVector1.GetNumBits()):
		res[i] = bitVector1[i]
	for i in range(bitVector2.GetNumBits()):
		res[i + bitVector1.GetNumBits()] = bitVector2[i]
	return res
def appendChemBoxBitVectors(chemNames, chemVectors, boxNames, boxVectors):
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
def getSimilarityFromBitVectors(bitVectors):
	simil = []
	count = 0
	n = 1
	for mol1 in range(len(bitVectors)):
		simil.append([1-x for x in DataStructs.BulkTanimotoSimilarity(bitVectors[mol1], bitVectors[:mol1+1])])
		if count >= 100:
			print(n * 100)
			n += 1
			count = 0
		count += 1
	return simil

def getDistanceMatrix(names, bitVectors):
	return _DistanceMatrix(names, getSimilarityFromBitVectors(bitVectors))

def distanceMatrixToTree(distanceMatrix, method = True):
	constructor = DistanceTreeConstructor()
	if method:
		return constructor.upgma(distanceMatrix)
	else:
		return constructor.nj(distanceMatrix)
def drawTree(tree):
	Phylo.draw_ascii(tree)
	#Phylo.draw_graphviz(tree)
	#pylab.show()
def genTreeBitVectors(names, vectors):
	dm = getDistanceMatrix(names, vectors)
	tree = distanceMatrixToTree(dm, True)
	drawTree(tree)
	return tree

def getChemBitVectorsArray(molecules):
	return [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in molecules]
def getChemThainingCompondsAsVectors():
	print(getFormatedTime() + " Downloading ideal training molecules for chemical clustering")
	chemNames, chemMolecules = getTrainingCompounds(inhibitorsChsIds, notInhibitorsChsIds)
	print(getFormatedTime() + " Obtaining bit vectors representation")
	return chemNames, getChemBitVectorsArray(chemMolecules)

def getDistanceVectors():
	print(getFormatedTime() + " Getting protein grid box and filtering active site atoms")
	filteredBox=filterBoxAtoms(args.proteinMol2, box, stepSize, minCavSize, topAtomsPersent, outBoxFile)
	print(getFormatedTime() + " Getting bit vector representation")
	return getMoleculesContactsAsBitVect(args.trainingMol2, filteredBox, bondLenClustering)


def composeBitVectorsToTree(namesC, vectorsC, namesD, vectorsD):
	print(getFormatedTime() + " Chem tree")
	genTreeBitVectors(namesC, vectorsC)
	print(getFormatedTime() + " Dist tree")
	genTreeBitVectors(namesD, vectorsD)

	print(getFormatedTime() + " Composing bit vectors")
	bothVectors, bothNames = appendChemBoxBitVectors(namesC, vectorsC, namesD, vectorsD)

	print(getFormatedTime() + " Final tree")
	genTreeBitVectors(bothNames, bothVectors)

def combineSimilarity(similarity1, similarity2, coeff):
	simil = []
	for row in range(len(similarity1)):
		newRow = []
		for item in range(len(similarity1[row])):
			newRow.append((similarity1[row][item] + similarity2[row][item] * coeff)/(1 + coeff))
		simil.append(newRow)
	return simil
def getTreeCombineSimilarityFromBitVectors(namesC, vectorsC, namesD, vectorsD):
	similC = getSimilarityFromBitVectors(vectorsC)
	print('chemTree')
	drawTree(distanceMatrixToTree(_DistanceMatrix(namesC, similC)))
	similD = getSimilarityFromBitVectors(vectorsD)
	print('PosTree')
	drawTree(distanceMatrixToTree(_DistanceMatrix(namesD, similD)))
	similboth = combineSimilarity(similC, similD, ratioCoeff)
	dm = _DistanceMatrix(namesC, similboth)
	tree = distanceMatrixToTree(dm, True)
	drawTree(tree)
	return tree
def outLowerTriangularDistanceMatrix(ofile, names, simil):
	fh = open(ofile, 'w')
	print('OutDistanceMatrix sm')
	print(len(namesD), file=fh)
	for i in range(len(namesD)):
		if i % 100 == 0:
			print(i)
		print (namesD[i], file=fh, end = '')
		for dist in similD[i][:-1]:
			print('\t', "%.4f" % dist, sep='', end = '', file=fh)
		print(file=fh)
	fh.close()

def genDistanceMatrixFileFewCompounds():
	namesD, vectorsD = getDistanceVectors()
	similD = getSimilarityFromBitVectors(vectorsD)
	del vectorsD
	gc.collect()
	outLowerTriangularDistanceMatrix(args.distanceMatrixOutput, namesD, similD)

def genDistanceMatrixFileManyCompounds(ofile, names, vectors):
	fh = open(ofile, 'w')
	print('OutDistanceMatrix obo')
	print(len(names), file=fh)

	for cnum1 in range(len(names)):
		if cnum1 % 100 == 0:
			print(cnum1)
		print (names[cnum1], file=fh, end = '')
		bsimil = [1-x for x in DataStructs.BulkTanimotoSimilarity(vectors[cnum1], vectors[:cnum1])]
		for sim in bsimil:
			print('\t', "%.4f" % sim, sep='', end = '', file=fh)
		print(file=fh)
	fh.close()

############################ CHEM bitVectors #####################################
namesC, vectorsC = getChemThainingCompondsAsVectors()
#namesC, vectorsC = getChemMoleculesAsBitVectorsOneByOne(args.trainingMol2)

genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, namesC, vectorsC)
############################################# Distance bitVectors ############################

############################################# params ############################
#gridSize = point3D(gridSizeX + bondLenBoxExtend, gridSizeY + bondLenBoxExtend, gridSizeZ + bondLenBoxExtend)
#centerCoords = point3D(asCenterX, asCenterY, asCenterZ)
#box = boxParams(centerCoords, gridSize)

############################################# Few compounds ############################
#genDistanceMatrixFileFewCompounds()
############################################# Few compounds ############################


#namesD, vectorsD = getDistanceVectors()
#genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, namesD, vectorsD)



#dmD = _DistanceMatrix(namesD, similD)
#print('Dm')

#del namesD
#del similD
#gc.collect()
#print('Now')
#drawTree(distanceMatrixToTree(dmD))



############################################# BOTH ############################

#composeBitVectorsToTree(namesC, vectorsC, namesD, vectorsD)
#getTreeCombineSimilarityFromBitVectors(namesC, vectorsC, namesD, vectorsD)
