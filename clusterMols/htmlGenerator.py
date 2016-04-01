from __future__ import print_function
from Cheetah.Template import Template

templateDef=''.join(open('template.html').readlines())
nameSpace = {'title': 'Example', 'maxSimilar':1, "results": {"ZINC12515747": { "score":-11.6, "similar":331},"ZINC05288595": { "score":-11.5, "similar":3321},"ZINC20580778": { "score":-11.3, "similar":3361},"ZINC09190718": { "score":-11.2, "similar":4331},"ZINC10157289": { "score":-11.1, "similar":54731}}}
t = Template(templateDef, searchList=[nameSpace])
print(t)