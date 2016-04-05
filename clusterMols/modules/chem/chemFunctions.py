import string

from chemspipy import ChemSpider

from rdkit import Chem
from rdkit.Chem import AllChem

from ..common.f import iterMol2

def getChemBitVectorsArrayFromMolecules(molecules, vectorSize = 1024):
	return [AllChem.GetMorganFingerprintAsBitVect(x,2,vectorSize) for x in molecules]
	
def getChemMoleculesAsBitVectorsOneByOne(file, vectorSize = 1024, limit = -1):
	vectors = []
	names = []
	readMolecules = 0
	for lines in iterMol2(file):
		if readMolecules % 100 == 0:
			print(readMolecules)
		name=lines[1].strip()
		vectors.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromMol2Block(''.join(lines)),2,vectorSize))
		names.append(name)
		readMolecules += 1
		if limit > 0 and readMolecules >= limit:
			break
	print(len(names))
	return names, vectors

def readMoleculesRd(file):
	molecules=[]
	names=[]
	for lines in iterMol2(file):
		molecules.append(Chem.MolFromMol2Block(''.join(lines)))
		names.append(lines[1])
	print(len(molecules))
	return molecules, names

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
	return namesY + namesN, moleculesY + moleculesN

def getLongNames(molsDict, pref = 4, onlyLettersDigits = False, token = '2228d430-a955-416b-b920-14547d28df9e'):
	cs = ChemSpider(token)
	names = {}
	for mol in molsDict:
		comp = cs.get_compound(mol[pref:])
		nName = comp.common_name.encode('ascii','ignore')
		if onlyLettersDigits:
			nName = leftOnlyLettersDigits(nName)
		names[mol] = nName
	return names