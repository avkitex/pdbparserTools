#!/usr/bin/python


class point3D:
    x = 0.0
    y = 0.0
    z = 0.0
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
############################ PARAMS #####################################
protH = open('3NTB_D.mol2')
centerCoords = point3D(26.9, 21.9, -76.9)
gridSize = point3D(30, 30, 30)
maxBondLen = 3.5
############################ PARAMS #####################################

class atom(point3D):
    serial = 0
    name = ''
    ctype = ''
    resnum = ''
    resname = ''
    charge = ''
    def __init__(self, serial, name, x, y, z, ctype, resnum, resname, charge):
        self.serial = serial
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.ctype = ctype
        self.resnum = resnum
        self.resname = resname
        self.charge = charge

class residue():#FIXME do i realy need this?
    atoms = []
    def __init__(self):
        self.atoms = []

proteinAtoms = []
curStr = ' '
while len(curStr) and "@<TRIPOS>ATOM" not in curStr:
    curStr = protH.next()
curStr = protH.next()
while len(curStr) and "@<TRIPOS>BOND" not in curStr:
#      1  N        -29.9490    7.8020  -41.7450 N.4     3  ALA3        0.0000
#      2  HT1      -29.3460    8.0040  -40.9730 H       3  ALA3        0.0000
    items = curStr.strip().split()
    try:
        proteinAtoms.append(atom(int(items[0]), items[1], float(items[2]), float(items[3]), float(items[4]), items[5], int(items[6]), items[7], float(items[8])))
    except Exception as e:
        print('Problems with atom. ', str(e))
    curStr = protH.next()
protH.close()
print(len(proteinAtoms))
