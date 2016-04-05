from __future__ import print_function
from datetime import datetime
import json

def getFormatedTime():
	return str(datetime.utcnow().strftime("[%H:%M:%S]"))
def logMsg(msg):
	print(getFormatedTime(), msg)
def errorMsg(msg):
	print(getFormatedTime(), "ERROR:", msg)


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

def dumpJson(file, data):
	with open(file, 'w') as out:
		json.dump(data, out)
	return 0
def loadJson(file):
	try:
		with open(file, 'r') as inp:
			data = json.load(inp)
		return data
	except Exception as e:
		errorMsg("Not able to read json file. " + str(e))
		return None