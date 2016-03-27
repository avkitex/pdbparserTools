from __future__ import print_function
import os.path, argparse, sys

from clusterize import getChemMoleculesAsBitVectorsOneByOne, genDistanceMatrixFileManyCompounds
#from clusterize import drawTree, distanceMatrixToTree


parser = argparse.ArgumentParser(prog='genTrainingChemDm.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-im', '--inputMol2', metavar='GlobConfig', type=str, help='Full path to multiMol2 training compounds docked file.', required=True)
parser.add_argument('-dm', '--distanceMatrixOutput', metavar='GlobConfig', type=str, help='Full path to distance matrix output file.', required=True)

args = parser.parse_args()

if not os.path.isfile(args.inputMol2):
	print('Docked training compounds mol2 file is not exists')
	sys.exit(1)

names, vectors = getChemMoleculesAsBitVectorsOneByOne(args.inputMol2)
#drawTree(distanceMatrixToTree(getDistanceMatrix(names, vectors)))
genDistanceMatrixFileManyCompounds(args.distanceMatrixOutput, names, vectors)

