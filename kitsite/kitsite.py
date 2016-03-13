#!/usr/bin/python
from __future__ import print_function

import math

############################ PARAMS #####################################
protinFile='3NTB_D.mol2'
outFile='res.xyz'
boxOutFile = 'box.xyz'

asCenterX = -26.9
asCenterY = 21.9
asCenterZ = -76.9

maxBondLen = 3.5

gridSizeX = 30
gridSizeY = 30
gridSizeZ = 30

stepSize = 0.5
minCavSize = 5
topAtoms = 250

cornerPenalty = 10
############################ PARAMS #####################################


############################ CONSTATNTS #################################
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
    def getActiveSiteAtomsX(self, step = stepSize, cavSize = minCavSize):
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
    def getActiveSiteAtomsY(self, step = stepSize, cavSize = minCavSize):
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
    def getActiveSiteAtomsZ(self, step = stepSize, cavSize = minCavSize):
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
    def getActiveSiteAtoms(self, step = stepSize, cavSize = minCavSize, topAtoms = topAtoms):
        self.getActiveSiteAtomsX()
        self.getActiveSiteAtomsY()
        self.getActiveSiteAtomsZ()
        self.atoms.sort(key=lambda atom: atom.score, reverse=True)
        for atom in self.atoms[:topAtoms]:
            print(atom.atom, atom.serial, atom.score, sep = ':')
        return sorted(self.atoms[:topAtoms], key=lambda atom: atom.serial)
        # for atom in self.atoms:
            # print(atom.atom, atom.serial, atom.score, sep = ':', end = ',')
            # if atom.score > treshold:
                # pickedAtoms.append(atom)
        # return sorted(pickedAtoms, key=lambda atom: atom.serial)




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
        return math.sqrt((self.x - atom.x)**2 + (self.y - atom.y)**2 + (self.z - atom.z)**2) - self.radius - atom.radius
    def aprint(self, fileH = None):
        if file:
            print(self.atom, self.x, self.y, self.z, sep='\t', file=fileH)
        else:
            print(self.atom, self.x, self.y, self.z, sep='\t')

def readProtein(protinFile):
    protH = open(protinFile)
    proteinAtoms = []
    curStr = protH.next()
    while len(curStr) and "@<TRIPOS>ATOM" not in curStr: #skipping info section
        curStr = protH.next()
    curStr = protH.next()
    #      1  N        -29.9490    7.8020  -41.7450 N.4     3  ALA3        0.0000
    #      2  HT1      -29.3460    8.0040  -40.9730 H       3  ALA3        0.0000
    while len(curStr) and "@<TRIPOS>BOND" not in curStr: # reading all atoms till to bonds section
        items = curStr.strip().split()
        try:
            proteinAtoms.append(atom(int(items[0]), items[1], float(items[2]), float(items[3]), float(items[4]), items[5], int(items[6]), items[7], float(items[8]), curStr))
        except Exception as e:
            print('Problems with atom. ', str(e))
        curStr = protH.next()
    protH.close()
    return proteinAtoms




############################################# params adopting ############################
gridSize = point3D(gridSizeX + maxBondLen, gridSizeY + maxBondLen, gridSizeZ + maxBondLen)
centerCoords = point3D(asCenterX, asCenterY, asCenterZ)
box = boxParams(centerCoords, gridSize)
############################################# params adopting ############################
############################################# main #######################################

pAtoms = readProtein(protinFile)
print(len(pAtoms))
box.extractAtoms(pAtoms)
if len(box.atoms) < 10:
    print('There are too little atoms in box! Check active site center or box size!')
print(len(box.atoms))
box.outBox(boxOutFile)
pickedX = box.getActiveSiteAtoms()
print(len(pickedX))
outfH = open(outFile, 'w')
print(len(pickedX), file=outfH)
print('Name', file=outfH)
for atom in pickedX:
    print(atom.atom, atom.x, atom.y, atom.z, sep='\t', file=outfH)
outfH.close()
