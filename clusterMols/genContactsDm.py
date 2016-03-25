from __future__ import print_function
import os.path, argparse, sys
from clusterize import point3D, boxParams

from clusterize import getDistanceVectors, genDistanceMatrixFileManyCompounds
#from clusterize import drawTree, distanceMatrixToTree


parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-im', '--inputMol2', metavar='GlobConfig', type=str, help='Full path to multiMol2 training compounds docked file.', required=True)
parser.add_argument('-pr', '--proteinMol2', metavar='GlobConfig', type=str, help='Full path to file with protein in mol2 format.', required=True)
parser.add_argument('-dm', '--distanceMatrixOutput', metavar='GlobConfig', type=str, help='Full path to distance matrix output file.', required=True)
#parser.add_argument('-bp', '--boxParams', metavar='GlobConfig', type=double, help='Box params file', required=True)
args = parser.parse_args()

if not os.path.isfile(args.inputMol2):
	print('Docked training compounds mol2 file is not exists')
	sys.exit(1)
if not os.path.isfile(args.proteinMol2):
	print('Protein mol2 file is not exists')
	sys.exit(1)

############################ PARAMS #####################################
asCenterX = -26.9
asCenterY = 21.9
asCenterZ = -76.9

gridSizeX = 30
gridSizeY = 30
gridSizeZ = 30
############################ PARAMS #####################################

box = boxParams(point3D(asCenterX, asCenterY, asCenterZ), point3D(gridSizeX, gridSizeY, gridSizeZ))

names, vectors = getDistanceVectors(args.proteinMol2, box, args.inputMol2)
#drawTree(distanceMatrixToTree(getDistanceMatrix(names, vectors)))
genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, names, vectors)
