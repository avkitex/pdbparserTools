#!/usr/bin/python

# sudo apt-get update

# sudo apt-get install -y python python-dev python-rdkit librdkit1 rdkit-data libfreetype6-dev python-networkx python-PyGraphviz
# sudo pip install biopython chemspipy

#deprecated pylab

# sudo apt-get install -y emboss embassy-phylip
# fneighbor -matrixtype l -datafile bigDM.dm -treetype n -outfile out

from __future__ import print_function

from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import SparseBitVect

from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

from modules.positional.contactsFunctions import *
from modules.chem.chemFunctions import *
from modules.common.f import errorMsg, logMsg

defaults = {}

defaults['bondLenClustering'] = 4.5

defaults['stepSize'] = 0.5
defaults['minCavSize'] = 4
defaults['topAtomsPersent'] = 20


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

def getChemThainingCompondsAsVectors(inhibitorsArr, notInhibitorsArr, longNames = False, vectorSize = 1024):
	logMsg("Downloading ideal training molecules for chemical clustering")
	chemNames, chemMolecules = getTrainingCompounds(inhibitorsArr, notInhibitorsArr, longNames)#onlyLettersDigits = False
	logMsg("Obtaining bit vectors representation")
	return chemNames, getChemBitVectorsArrayFromMolecules(chemMolecules, vectorSize)

def getProteinActiveSiteAtoms(proteinFile, box, topAtomsPersent = defaults['topAtomsPersent'], stepSize = defaults['stepSize'], minCavSize = defaults['minCavSize'], outBoxFile = ''):
	logMsg("Getting protein grid box and filtering active site atoms")
	return filterBoxAtoms(proteinFile, box, stepSize, minCavSize, topAtomsPersent, outBoxFile)
	
def getProteinContactsAsBitVectors(filteredBox, moleculesFile, bondLenClustering = defaults['bondLenClustering'], limit=-1):
	logMsg("Getting contacts with protein bit vector representation")
	return getMoleculesContactsAsBitVect(moleculesFile, filteredBox, bondLenClustering, limit = limit)


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

def genDistanceMatrixFileManyCompounds(ofile, names, vectors):
	fh = open(ofile, 'w')
	logMsg('Generating and outputting distance matrix')
	print(len(names), file=fh)

	for cnum1 in range(len(names)):
		if cnum1 % 100 == 0:
			print(cnum1)
		print (names[cnum1], file=fh, end = '')
		bsimil = [1 - x for x in DataStructs.BulkTanimotoSimilarity(vectors[cnum1], vectors[:cnum1])]
		for sim in bsimil:
			print('\t', "%.4f" % sim, sep='', end = '', file=fh)
		print(file=fh)
	fh.close()

def getTrainigToBaseSimilarityMatrix(trainMolsDict, namesBase, vectorsBase):
	result = {}
	for tname in trainMolsDict:
		result[tname] = dict(zip(namesBase, DataStructs.BulkTanimotoSimilarity(trainMolsDict[tname], vectorsBase)))
	return result
def getCompoundToSetSimilarity(vector, molsDict):
	res = {}
	for i in molsDict:
		res[i] =  DataStructs.TanimotoSimilarity(vector, molsDict[i])
	return res
			
			
			
			
	