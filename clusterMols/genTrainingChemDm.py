import os.path, argparse, sys

from clusterize import getChemThainingCompondsAsVectors, genDistanceMatrixFileManyCompounds
from clusterize import drawTree, distanceMatrixToTree

parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-dm', '--distanceMatrixOutput', metavar='GlobConfig', type=str, help='Full path to distance matrix output file.', required=True)
#parser.add_argument('-inhl', '--inhibitorsList', metavar='GlobConfig', type=str, help='inhibitors list.', required=True)
#parser.add_argument('-ninhl', '--notInhibitorsList', metavar='GlobConfig', type=str, help='Not inhibitors list.', required=True)
args = parser.parse_args()


############################ PARAMS #####################################
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


############################ CHEM bitVectors #####################################
names, vectors = getChemThainingCompondsAsVectors(inhibitorsChsIds, notInhibitorsChsIds, False)
#drawTree(distanceMatrixToTree(getDistanceMatrix(namesC, vectorsC)))
genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, names, vectors)
