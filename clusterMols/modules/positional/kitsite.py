#!/usr/bin/python
from __future__ import print_function

import math
from ..common.f import iterMol2
from rdkit.DataStructs.cDataStructs import SparseBitVect


############################ CONSTATNTS #################################
defaultTopAtomsPersent = 17
minCavSize = 3.5
defaultStepSize = 0.5
bondLenBoxExtend = 0
bondLenClustering = 3.5
cornerPenalty = 12


redius = {}
redius['H'] = 1.2
redius['C'] = 1.7
redius['N'] = 1.55
redius['O'] = 1.52
redius['F'] = 1.47
redius['P'] = 1.8
redius['S'] = 1.8
redius['Cl'] = 1.75
redius['Br'] = 1.85
redius['I'] = 1.98
############################ CONSTATNTS #################################

class point3D:
	x = 0.0
	y = 0.0
	z = 0.0
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
class boxParams():
	center = None
	size = None
	atoms = []
	def __init__(self, center, size):
		self.size = size
		self.center = center
		self.minx = self.center.x - self.size.x / 2
		self.maxx = self.center.x + self.size.x / 2
		self.miny = self.center.y - self.size.y / 2
		self.maxy = self.center.y + self.size.y / 2
		self.minz = self.center.z - self.size.z / 2
		self.maxz = self.center.z + self.size.z / 2
		self.atoms = []
	def atomIn(self, atom):
		return self.minx <= atom.x <= self.maxx and \
			   self.miny <= atom.y <= self.maxy and \
			   self.minz <= atom.z <= self.maxz
	def extractAtoms(self, atoms):
		self.atoms = []
		for atom in atoms:
			if self.atomIn(atom):
				self.atoms.append(atom)
	def outBox(self, outFile):
		outH = open(outFile, 'w')
		print(len(self.atoms), file=outH)
		print('ProteinBox', file=outH)
		for atom in self.atoms:
			atom.aprint(outH)
		outH.close()
	def getActiveSiteAtomsX(self, step = defaultStepSize, minCavSize = minCavSize):
		self.atoms.sort(key=lambda atom: atom.x)
		ycoord = self.miny
		while ycoord < self.maxy:
			prevAtom = -1
			zcoord = self.minz
			while zcoord < self.maxz:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineX(ycoord, zcoord):
						if prevAtom >= 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				self.atoms[prevAtom].score -= cornerPenalty
				zcoord += step

			prevAtom = -1
			xcoord = self.maxx
			while xcoord > self.minx:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineZ(xcoord, ycoord):
						if prevAtom >= 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				self.atoms[prevAtom].score -= cornerPenalty
				xcoord -= step

			ycoord += step
	def getActiveSiteAtomsY(self, step = defaultStepSize, cavSize = minCavSize):
		self.atoms.sort(key=lambda atom: atom.y)
		zcoord = self.minz
		while zcoord < self.maxz:
			prevAtom = -1
			xcoord = self.minx
			while xcoord < self.maxx:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineY(xcoord, zcoord):
						if prevAtom > 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							#print(dist)
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				xcoord += step
				self.atoms[prevAtom].score -= cornerPenalty

			prevAtom = -1
			ycoord = self.maxy
			while ycoord > self.miny:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineX(ycoord, zcoord):
						if prevAtom > 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							#print(dist)
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				ycoord -= step
				self.atoms[prevAtom].score -= cornerPenalty

			zcoord += step
	def getActiveSiteAtomsZ(self, step = defaultStepSize, cavSize = minCavSize):
		self.atoms.sort(key=lambda atom: atom.z)
		xcoord = self.minx
		while xcoord < self.maxx:
			prevAtom = -1
			ycoord = self.miny
			while ycoord < self.maxy:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineZ(xcoord, ycoord):
						if prevAtom > 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							#print(dist)
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				ycoord += step
				self.atoms[prevAtom].score -= cornerPenalty

			prevAtom = -1
			zcoord = self.maxz
			while zcoord > self.minz:
				for atomNum in range(len(self.atoms)):
					if self.atoms[atomNum].crossLineY(xcoord, zcoord):
						if prevAtom > 0:
							dist = self.atoms[atomNum].dist(self.atoms[prevAtom])
							#print(dist)
							if dist >= minCavSize:
								self.atoms[atomNum].score += 1
								self.atoms[prevAtom].score += 1
						else:
							self.atoms[atomNum].score -= cornerPenalty
						prevAtom = atomNum
				zcoord -= step
				self.atoms[prevAtom].score -= cornerPenalty

			xcoord += step
	def getActiveSiteAtoms(self, step = defaultStepSize, cavSize = minCavSize, topAtomsPersent = defaultTopAtomsPersent):
		self.getActiveSiteAtomsX(step = step)
		self.getActiveSiteAtomsY(step = step)
		self.getActiveSiteAtomsZ(step = step)
		self.atoms.sort(key=lambda atom: atom.score, reverse=True)
		atomsAmount=topAtomsPersent*len(self.atoms)//100
		#for atom in self.atoms[:atomsAmount]:
		#	print(atom.atom, atom.serial, atom.score, sep = ':')
		return sorted(self.atoms[:atomsAmount], key=lambda atom: atom.serial)
		# for atom in self.atoms:
			# print(atom.atom, atom.serial, atom.score, sep = ':', end = ',')
			# if atom.score > treshold:
				# pickedAtoms.append(atom)
		# return sorted(pickedAtoms, key=lambda atom: atom.serial)


def twoPointsDist(x1, y1, z1, x2, y2, z2):
	return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

class atom(point3D):
	serial = 0
	name = ''
	ctype = ''
	resnum = ''
	resname = ''
	charge = ''
	radius = 0
	atom = ''
	score = 0
	mol2String = ''
	def __init__(self, serial, name, x, y, z, ctype, resnum, resname, charge, mol2String = ''):
		self.serial = serial
		self.name = name
		self.x = x
		self.y = y
		self.z = z
		self.ctype = ctype
		point = ctype.find('.')
		if point != -1:
			self.atom = ctype[:point]
		else:
			self.atom = ctype
		if self.atom not in redius:
			self.radius = redius['C']
		else:
			self.radius = redius[self.atom]
		self.resnum = resnum
		self.resname = resname
		self.charge = charge
		self.score = 0
		self.mol2String = mol2String
	def crossLineX(self, y, z):
		return self.y - self.radius <= y <= self.y + self.radius and \
			   self.z - self.radius <= z <= self.z + self.radius
	def crossLineY(self, x, z):
		return self.x - self.radius <= x <= self.x + self.radius and \
			   self.z - self.radius <= z <= self.z + self.radius
	def crossLineZ(self, x, y):
		return self.x - self.radius <= x <= self.x + self.radius and \
			   self.y - self.radius <= y <= self.y + self.radius
	def dist(self, atom):
		return twoPointsDist(self.x, self.y, self.z, atom.x, atom.y, atom.z) - self.radius - atom.radius
	def aprint(self, fileH = None):
		if file:
			print(self.atom, self.x, self.y, self.z, sep='\t', file=fileH)
		else:
			print(self.atom, self.x, self.y, self.z, sep='\t')
class singleMolecule():
	atoms = []
	def __init__(self, atoms):
		self.atoms = atoms
		atoms.sort(key=lambda atom: atom.x)
		self.minx = atoms[0].x
		self.maxx = atoms[-1].x
		atoms.sort(key=lambda atom: atom.y)
		self.miny = atoms[0].y
		self.maxy = atoms[-1].y
		atoms.sort(key=lambda atom: atom.z)
		self.minz = atoms[0].z
		self.maxz = atoms[-1].z
	def minDist(self, sAtom):
		dist = 100000
		for atom in self.atoms:
			d = twoPointsDist(atom.x, atom.y, atom.z, sAtom.x, sAtom.y, sAtom.z)
			dist = min(dist, d)
		return dist
def atomFromMol2String(mol2String):
	items = mol2String.strip().split()
	try:
		return atom(int(items[0]), items[1], float(items[2]), float(items[3]), float(items[4]), items[5], int(items[6]), items[7], float(items[8]), mol2String)
	except Exception as e:
		print('Problems with atom. ', str(e))
		return None
def readMol2ProteinMolecule(protinFile):
	protH = open(protinFile)
	proteinAtoms = []
	curStr = protH.next()
	while len(curStr) and "@<TRIPOS>ATOM" not in curStr: #skipping info section
		curStr = protH.next()
	curStr = protH.next()
	#	  1  N		-29.9490	7.8020  -41.7450 N.4	 3  ALA3		0.0000
	#	  2  HT1	  -29.3460	8.0040  -40.9730 H	   3  ALA3		0.0000
	while len(curStr) and "@<TRIPOS>BOND" not in curStr: # reading all atoms till to bonds section
		atom = atomFromMol2String(curStr)
		if atom:
			proteinAtoms.append(atom)
		curStr = protH.next()
	protH.close()
	return proteinAtoms


############################################# main #######################################
def filterBoxAtoms(protInFile, box, step, topAtomsPersent = defaultTopAtomsPersent):
	pAtoms = readMol2ProteinMolecule(protInFile)
	print('Total protein atoms:', len(pAtoms))
	box.extractAtoms(pAtoms)
	if len(box.atoms) < 10:
		print('There are too little atoms in box! Check active site center or box size!')
	print('Box atoms:', len(box.atoms))
	#box.outBox(boxOutFile)
	activeSiteAtoms = box.getActiveSiteAtoms(topAtomsPersent = topAtomsPersent)
	print('Active site atoms: ' + str(len(activeSiteAtoms)) +  ', (' + str(topAtomsPersent) + '%)')
	return activeSiteAtoms


	#outfH = open(outFile, 'w')
	#print(len(pickedX), file=outfH)
	#print('Name', file=outfH)
	#for atom in pickedX:
	#	print(atom.atom, atom.x, atom.y, atom.z, sep='\t', file=outfH)
	#outfH.close()

def getContactsAsBitVec(molecule, filteredBoxAtoms, contact = bondLenClustering):
	bv = SparseBitVect(len(filteredBoxAtoms))
	for atomNum in range(len(filteredBoxAtoms)):
		if molecule.minx - contact  <= filteredBoxAtoms[atomNum].x <= molecule.maxx + contact and \
		   molecule.miny - contact <= filteredBoxAtoms[atomNum].y <= molecule.maxy + contact and \
		   molecule.minz - contact <= filteredBoxAtoms[atomNum].z <= molecule.maxz + contact:
			if molecule.minDist(filteredBoxAtoms[atomNum]) < contact:
				bv.SetBit(atomNum)
	return bv


def getMoleculesContactsAsBitVect(file, filteredBoxAtoms):
	filteredBoxAtoms.sort(key=lambda atom: (atom.x, atom.y, atom.z))
	contactsBV = []
	names = []
	for lines in iterMol2(file):
		name = lines[1].strip()
		atoms = []
		curline = 0
		while curline < len(lines) and "@<TRIPOS>ATOM" not in lines[curline]: #skipping info section
			curline += 1
		curline += 1
		while curline < len(lines) and "@<TRIPOS>BOND" not in lines[curline]: # reading all atoms till to bonds section
			atom = atomFromMol2String(lines[curline])
			if atom:
				atoms.append(atom)
			curline += 1
		contactsBV.append(getContactsAsBitVec(singleMolecule(atoms), filteredBoxAtoms))
		names.append(name)

	return contactsBV, names
