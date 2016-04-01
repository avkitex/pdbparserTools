from __future__ import print_function
from Cheetah.Template import Template
#{"ZINC12515747": { "energy":-11.6, "similar":331},"ZINC05288595": { "energy":-11.5, "similar":3321},"ZINC20580778": { "energy":-11.3, "similar":3361},"ZINC09190718": { "energy":-11.2, "similar":4331},"ZINC10157289": { "energy":-11.1, "similar":54731}}
#"similarChem":[]
#"similarContacts":[]
def createHtmlReport(outFile, results, maxsimil=2, maxOutComp = 100):
	templateDef=''.join(open('template.html').readlines())
	nameSpace = {'maxSimilar':1, "results": results[:maxOutComp]}
	t = Template(templateDef, searchList=[nameSpace])
	oHandle=open(outFile, 'w')
	print(t, file=oHandle)
	oHandle.close()