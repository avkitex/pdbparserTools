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
