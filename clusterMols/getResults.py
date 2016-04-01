#!/usr/bin/python
from __future__ import print_function
import os.path, argparse, sys
from copy import deepcopy

from clusterize import getChemMoleculesAsBitVectorsOneByOne, getChemThainingCompondsAsVectors, getTrainigToBaseSimilarityMatrix, getProteinContactsAsBitVectors, getCompoundToSetSimilarity, errorMsg
from clusterize import point3D, boxParams
from clusterize import drawTree, distanceMatrixToTree, getDistanceMatrix
from htmlGenerator import createHtmlReport

parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-bm', '--baseMol2', metavar='GlobConfig', type=str, help='Full path to database multiMol2', required=True)
parser.add_argument('-td', '--trainingDockedMol2', metavar='GlobConfig', type=str, help='Full path to docked training multiMol2', required=True)
parser.add_argument('-bd', '--baseDockedMol2', metavar='GlobConfig', type=str, help='Full path to database docked multiMol2', required=True)
parser.add_argument('-bdr', '--baseDockingResult', metavar='GlobConfig', type=str, help='Full path to database docking result tsv file', required=True)
parser.add_argument('-pr', '--proteinMol2', metavar='GlobConfig', type=str, help='Full path to file with protein in mol2 format.', required=True)
#parser.add_argument('-bp', '--boxParams', metavar='GlobConfig', type=double, help='Box params file', required=True)


#parser.add_argument('-inhl', '--inhibitorsList', metavar='GlobConfig', type=str, help='inhibitors list.', required=True)
#parser.add_argument('-ninhl', '--notInhibitorsList', metavar='GlobConfig', type=str, help='Not inhibitors list.', required=True)
args = parser.parse_args()

averageIsolatedDist = 'average'

if not os.path.isfile(args.baseMol2):
	print('Base compounds mol2 file does not exists')
	sys.exit(1)
if not os.path.isfile(args.trainingDockedMol2):
	print('Docked training compounds mol2 file does not exists')
	sys.exit(1)
if not os.path.isfile(args.baseDockedMol2):
	print('Base docked compounds mol2 file does not exists')
	sys.exit(1)
if not os.path.isfile(args.baseDockingResult):
	print('Base docking results file does not exists')
	sys.exit(1)
if not os.path.isfile(args.proteinMol2):
	print('Protein mol2 file is not exists')
	sys.exit(1)
	
#if not os.path.isfile(args.inhibitorsList):
#	print('Docked training compounds mol2 file is not exists')
#if not os.path.isfile(args.notInhibitorsList):
#	print('Docked training compounds mol2 file is not exists')

############################ PARAMS #####################################
inhibitorsChsIds=[331, 1353, 1906, 2000, 2066, 2121, 2157, 2562, 2925, 3065, 3097, 3192, 3225,\
3254, 3544, 3584, 3693, 3694, 3897, 3904, 4339, 4393, 4480, 4617, 4911, 5304, 5308, 8711, 133236,\
28714, 29843, 32983, 33051, 34911, 65084, 110209, 137720, 147913, 171462, 388601, 392692, 392977,\
393569, 394812, 394942, 581047, 8082544, 1265915, 2871682, 3512137, 4444105, 4444112, 4445953, \
4510145, 5365258, 8081394, 8098374, 8443182, 8530675, 9569918, 10442653, 10442740, 20572535, \
21106585, 23122887, 23122889, 23122978, 25057753] # 12951165 no results case of bor
notInhibitorsChsIds=[650, 682, 733, 864, 971, 1116, 1512, 1710, 2971, 3350, 5611, 5653, 5764, 5768, 6170, 6257, 6312, 7742, 10610, 73505, 83361, 96749, 111188, 120261, 216840, 391555,\
392800, 4576521, 8257952, 20572534, 23122865, 23122927, 23123076, 112728, 10442445]

asCenterX = -26.9
asCenterY = 21.9
asCenterZ = -76.9

gridSizeX = 30
gridSizeY = 30
gridSizeZ = 30
############################ PARAMS #####################################

def averageDistances(distances, subsetTraining):
	if subsetTraining[0] not in distances:
		print('There is no ' + subsetTraining[0] + ' in distances matrix')
		return []
	out = deepcopy(distances[subsetTraining[0]])
	for tname in subsetTraining[1:]:
		if tname not in distances:
			print('There is no ' + tname + ' in distances matrix')
			continue
		for compName in distances[tname]:
			out[compName] += distances[tname][compName]
	for compName in out:
		out[compName] /= len(subsetTraining)
	return out
def isolatedDistances(distances, subsetTraining, amount):
	amount = max(amount, 10)
	res = set()
	for tname in subsetTraining:
		print(tname, '******')
		if tname not in distances:
			print('There is no ' + tname + ' in distances matrix')
			continue
		out = deepcopy(distances[tname])
		for i in sorted(out.keys(), key=out.get, reverse=True)[:amount]:
			res.add(i)
			print(i, out[i])
	return list(res)
def getTopSimilarCompounds(distances, subsetTraining, topPer = 10.):
	if averageIsolatedDist == 'average':
		averageDist = averageDistances(distances, subsetTraining)
		return sorted(averageDist.keys(), key=averageDist.get, reverse=True)[:int(topPer/100. * len(averageDist))]
	else:
		return isolatedDistances(distances, subsetTraining, int(topPer/100. * len(distances[subsetTraining[0]]) / len(subsetTraining)))
def filterVectorsByNamesList(names, vectors, list):
	comp = dict(zip(names, vectors))
	res = []
	for name in list:
		res.append(comp[name])
	return res
def dfs(curClade, listNotTerminals):
	if curClade.is_terminal():
		return
	else:
		for subClade in curClade:
			dfs(subClade, listNotTerminals)
		listNotTerminals.append(curClade)
def findYesClades(tree, threshold):
	allNotTerminals = []
	dfs(tree.clade, allNotTerminals)
	for subTree in allNotTerminals:
		pers = 0
		terminals = subTree.get_terminals()
		names = []
		for terminal in terminals:
			if terminal.name[:3].lower() == 'yes':
				pers += 1
			names.append(terminal.name)
		pers = pers * 1. / len(terminals)
		#print(subTree.name, pers, names)
		if pers >= threshold / 100.:
			yield names
def chemTreeWalkGenerator(namesTrain, vectorsTrain, namesBase, vectorsBase, minYesThreshold= 75, topPercent = 2):
	chemTrainMolsDict = dict(zip(namesTrain, vectorsTrain))
	chemTrainMolsTree = distanceMatrixToTree(getDistanceMatrix(namesTrain, vectorsTrain))
	#drawTree(chemTrainMolsTree)
	#print(chemTrainMolsTree)
	############################ Base bitVectors #####################################

	distances = getTrainigToBaseSimilarityMatrix(chemTrainMolsDict, namesBase, vectorsBase)
	#maxSimilarNames = set()
	for subset in findYesClades(chemTrainMolsTree, minYesThreshold):
		#print(subset)
		topBaseLikeTrainigSetNames = getTopSimilarCompounds(distances, subset, topPer = topPercent)
		yield subset, topBaseLikeTrainigSetNames
		#for name in topBaseLikeTrainigSetNames:
		#	maxSimilarNames.add(name)
		#tree2Names = subset + topBaseLikeTrainigSetNames
		
		#tree2Vectors = filterVectorsByNamesList(namesTrain, vectorsTrain, subset) + filterVectorsByNamesList(namesBase, vectorsBase, topBaseLikeTrainigSetNames)
		#tree2 = distanceMatrixToTree(getDistanceMatrix(tree2Names, tree2Vectors))
		#drawTree(tree2)
	#print(list(maxSimilarNames))

def contactsTreeWalkGenerator(namesTrain, vectorsTrain, namesBase, vectorsBase, minYesThreshold= 60, topPercent = 5):
	contactsTrainMolsDict = dict(zip(namesTrain, vectorsTrain))
	contactsTrainMolsTree = distanceMatrixToTree(getDistanceMatrix(namesTrain, vectorsTrain))
	#drawTree(contactsTrainMolsTree)
	distances = getTrainigToBaseSimilarityMatrix(contactsTrainMolsDict, namesBase, vectorsBase)

	for subset in findYesClades(contactsTrainMolsTree, minYesThreshold):
		#print(subset)
		topBaseLikeTrainigSetNames = getTopSimilarCompounds(distances, subset, topPer = topPercent)
		yield subset, topBaseLikeTrainigSetNames

def readDockingResults(resFile):
	fHandle = open(resFile)
	result= {}
	try:
		header = fHandle.next().strip().split()
		if "Ligand" not in header or "1_Energy" not in header:
			raise
	except:
		errorMsg("Docking result file is in bad format")
		return {}
	for line in fHandle:
		entry = dict(zip(header, line.strip().split()))
		result[entry["Ligand"]] = float(entry["1_Energy"].replace(',','.'))
	fHandle.close()
	return result
def getBestSimilarToMolset(vector, bolsDict):
	resDict = getCompoundToSetSimilarity(vector, bolsDict)
	return [{"id":x[4:], "type": "inhibitor" if x[:3].lower() == 'yes' else "notinhibitor", "similarity":resDict[x]} for x in sorted(resDict.keys(), key = resDict.get, reverse=True)]
############################ main #####################################

baseDockingResults = readDockingResults(args.baseDockingResult)

namesTrainChem, vectorsTrainChem = getChemThainingCompondsAsVectors(inhibitorsChsIds, notInhibitorsChsIds, False)
chemTrainMolsDict = dict(zip(namesTrainChem, vectorsTrainChem))
namesBaseChem, vectorsBaseChem = getChemMoleculesAsBitVectorsOneByOne(args.baseMol2)
chemBaseMolsDict = dict(zip(namesBaseChem, vectorsBaseChem))


box = boxParams(point3D(asCenterX, asCenterY, asCenterZ), point3D(gridSizeX, gridSizeY, gridSizeZ))

namesTrainContacts, vectorsTrainContacts = getProteinContactsAsBitVectors(args.proteinMol2, box, args.trainingDockedMol2)
contactsTrainMolsDict = dict(zip(namesTrainContacts, vectorsTrainContacts))
namesBaseContacts, vectorsBaseContacts = getProteinContactsAsBitVectors(args.proteinMol2, box, args.baseDockedMol2)
contactsBaseMolsDict = dict(zip(namesBaseContacts, vectorsBaseContacts))

finalResult = set()

for contactsTrainSet, contactsBaseSet in contactsTreeWalkGenerator(namesTrainContacts, vectorsTrainContacts, namesBaseContacts, vectorsBaseContacts, topPercent = 1):
	print("Analyzing group:")
	print(contactsTrainSet)
	#print('Subset1')
	#print(contactsBaseSet)
	bestCompounds = set()
	
	trainSubset = []
	baseSubset = []
	for x in contactsTrainSet:
		if x in chemTrainMolsDict:
			trainSubset.append(chemTrainMolsDict[x])
	for x in contactsBaseSet:
		if x in chemBaseMolsDict:
			baseSubset.append(chemBaseMolsDict[x])
	for chemTrainSet, chemBaseSet in chemTreeWalkGenerator(contactsTrainSet, trainSubset, contactsBaseSet, baseSubset, topPercent = 30):
		for cName in chemBaseSet:
			bestCompounds.add(cName)
	#print('Best compounds')
	#print(bestCompounds)
	bestCompoundsList = list(bestCompounds)
	bestCompoundsList.sort(key = lambda x: baseDockingResults[x])
	for cName in bestCompoundsList[:20]:
		#print('  ', cName, baseDockingResults[cName])
		finalResult.add(cName)
	#print("*******************")

finalResultList = list(finalResult)
finalResultList.sort(key = lambda x: baseDockingResults[x])
results = []
for cName in finalResultList:
	print('  ', cName, baseDockingResults[cName])
	results.append({"name":cName, "energy":baseDockingResults[cName], "similarChem":getBestSimilarToMolset(chemBaseMolsDict[cName], chemTrainMolsDict), "similarContacts":getBestSimilarToMolset(contactsBaseMolsDict[cName], contactsTrainMolsDict)})

createHtmlReport('out.html', results)




