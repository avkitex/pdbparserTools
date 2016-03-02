#!/usr/bin/python

from rdkit import Chem

def iterMol2(file):
    handle = open(file, 'r')
    if handle:
        curMol2 = []
        for line in handle:
            if '@<TRIPOS>MOLECULE' in line:
                if len(curMol2):
                    yield curMol2
                    curMol2 = []
            curMol2.append(line)
        if len(curMol2):
            yield curMol2
        handle.close()

def clusterMolecules(file):
    molecules=[]
    for file in iterMol2(file):
        molecules.append(Chem.MolFromMol2Block(''.join(file)))
    print(len(molecules))
