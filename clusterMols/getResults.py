from __future__ import print_function
import os.path, argparse, sys

from clusterize import getChemMoleculesAsBitVectorsOneByOne, getChemThainingCompondsAsVectors, genDistanceMatrixFileManyCompounds, getDistTrainigToBaseDistMatrix
from clusterize import drawTree, distanceMatrixToTree, getDistanceMatrix

parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-bm', '--baseMol2', metavar='GlobConfig', type=str, help='Full path to database multiMol2 compounds file.', required=True)

#parser.add_argument('-inhl', '--inhibitorsList', metavar='GlobConfig', type=str, help='inhibitors list.', required=True)
#parser.add_argument('-ninhl', '--notInhibitorsList', metavar='GlobConfig', type=str, help='Not inhibitors list.', required=True)
args = parser.parse_args()


if not os.path.isfile(args.baseMol2):
	print('Docked training compounds mol2 file does not exists')
	sys.exit(1)
	
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

def averageDistances(distances, list):
	if list[0] not in distances:
		print('There is no ' + list[0] + ' in distances matrix')
		return []
	out = distances[list[0]]
	for tname in list[1:]:
		if tname not in distances:
			print('There is no ' + tname + ' in distances matrix')
			continue
		for bcomp in range(len(distances[tname])):
			out[bcomp] += distances[tname][bcomp]
	for i in out:
		i /= len(list)
	return out

def getTopSimilarCompounds(distances, namesBase, subsetTraining, topPer = 10):
	averageDist = zip(namesBase, averageDistances(distances, subsetTraining))
	averageDist.sort(key=lambda x: -x[1])
	res = []
	for cname in averageDist[:int(topPer/100. * len(averageDist))]:
		#print(cname)
		res.append(cname[0])
	return res
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
############################ Training bitVectors #####################################
namesTrain, vectorsTrain = getChemThainingCompondsAsVectors(inhibitorsChsIds, notInhibitorsChsIds, False)
chemTrainMolsDict = dict(zip(namesTrain, vectorsTrain))
trainChemTree = distanceMatrixToTree(getDistanceMatrix(namesTrain, vectorsTrain))
#drawTree(trainChemTree)
#print(trainChemTree)
############################ Base bitVectors #####################################

namesBase, vectorsBase = getChemMoleculesAsBitVectorsOneByOne(args.baseMol2, limit=5000)
distances = getDistTrainigToBaseDistMatrix(namesTrain, vectorsTrain, namesBase, vectorsBase)
maxSimilarNames = set()
for subset1 in findYesClades(trainChemTree, 75):
	#print(subset1)
	topBaseLikeTrainigSetNames = getTopSimilarCompounds(distances, namesBase, subset1, topPer = 2)
	for name in topBaseLikeTrainigSetNames:
		maxSimilarNames.add(name)
	#tree2Names = subset1 + topBaseLikeTrainigSetNames
	
	#tree2Vectors = filterVectorsByNamesList(namesTrain, vectorsTrain, subset1) + filterVectorsByNamesList(namesBase, vectorsBase, topBaseLikeTrainigSetNames)
	#tree2 = distanceMatrixToTree(getDistanceMatrix(tree2Names, tree2Vectors))
	#drawTree(tree2)
print(list(maxSimilarNames))
