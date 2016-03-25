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
#from Bio.Phylo import draw

from modules.positional.kitsite import *
from modules.chem.mol2Reader import *

parser = argparse.ArgumentParser(prog='genTestDistDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-tf', '--inputMol2', metavar='GlobConfig', type=str, help='Full path to multiMol2 training compounds docked file.', required=True)
parser.add_argument('-pr', '--proteinMol2', metavar='GlobConfig', type=str, help='Full path to file with protein in mol2 format.', required=True)
parser.add_argument('-dm', '--distanceMatrixOutput', metavar='GlobConfig', type=str, help='Full path to distance matrix output file.', required=True)
#parser.add_argument('-bp', '--boxParams', metavar='GlobConfig', type=double, help='Box params file', required=True)
args = parser.parse_args()


if not os.path.isfile(args.inputMol2):
	print('Docked training compounds mol2 file is not exists')
if not os.path.isfile(args.proteinMol2):
	print('Protein mol2 file is not exists')

############################ PARAMS #####################################
asCenterX = -26.9
asCenterY = 21.9
asCenterZ = -76.9

gridSizeX = 30
gridSizeY = 30
gridSizeZ = 30

bondLenBoxExtend = 0
bondLenClustering = 4.5

stepSize = 0.5
minCavSize = 4
topAtomsPersent = 20

############################ PARAMS #####################################

def getFormatedTime():
	return str(datetime.utcnow().strftime("[%H:%M:%S]"))

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
def getDistanceVectors(proteinFile, box, inputMol2):
	print(getFormatedTime() + " Getting protein grid box and filtering active site atoms")
	filteredBox = filterBoxAtoms(proteinFile, box, stepSize, minCavSize, topAtomsPersent)
	print(getFormatedTime() + " Getting bit vector representation")
	return getMoleculesContactsAsBitVect(inputMol2, filteredBox, bondLenClustering)


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

gridSize = point3D(gridSizeX + bondLenBoxExtend, gridSizeY + bondLenBoxExtend, gridSizeZ + bondLenBoxExtend)
centerCoords = point3D(asCenterX, asCenterY, asCenterZ)
box = boxParams(centerCoords, gridSize)


names, vectors = getDistanceVectors(args.proteinMol2, box, args.inputMol2)
#drawTree(distanceMatrixToTree(getDistanceMatrix(names, vectors)))
genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, names, vectors)
