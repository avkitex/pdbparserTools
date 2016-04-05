import sys, argparse, os
from modules.common.f import loadJson
from htmlGenerator import createHtmlReport

parser = argparse.ArgumentParser(prog='htmlReportRegenerator.py', usage='%(prog)s [options]', description='Regeneration of html report',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-inj', '--inputJson', type=str, help='Input saved fson info', required=True)
parser.add_argument('-htmlrep', '--htmlreport', type=str, default='report.html',  help='Html out file.', required=True)
parser.add_argument('-tp', '--template', type=str, default='template.html', help='Path to template')
parser.add_argument('-mtop', '--maxtop', type=int, default=150,  help='Max report top compounds.')
parser.add_argument('-msim', '--maxsimilar', type=int, default=3,  help='Max report similar compounds.')

args = parser.parse_args()

if not os.path.isfile(args.template):
	print('Not able to locate template file')
	sys.exit(1)
if not os.path.isfile(args.inputJson):
	print('Not able to locate input json file')
	sys.exit(1)
	
results=loadJson(args.inputJson)
	
createHtmlReport(args.htmlreport, results, maxsimil=min(max(args.maxsimilar, 3), 10), maxOutComp = args.maxtop)