import os
import sys
import pandas as pd
import numpy as np

helper = "call by ~ <metadata> <inName> (<ratio>) (<outname>)\n"
helper += "\tget the 30-mer counting table from intermediate files from running mergeKMCTable.exe\n"
helper += "filter 30-mer that appears more than 60% in either Response or Nonresponse group\n"
helper += "output the OTU-like tbale as TSV table"

if len(sys.argv) == 3:
    dname = sys.argv[1]
    inName = sys.argv[2]
    ratio = 0.6
    outname = "K30"
elif len(sys.argv) == 4:
    dname = sys.argv[1]
    inName = sys.argv[2]
    ratio = int(sys.argv[3]) / 100.0
    outname = "K30"
elif len(sys.argv) == 5:
    dname = sys.argv[1]
    inName = sys.argv[2]
    ratio = int(sys.argv[3]) / 100.0
    outname = sys.argv[4]
else:
    print(helper)
    exit()

sample_group = {}
CB_MAPPING = {"PR":1, "CR":1, "PD":0, "DD":0, "SD":1}
with open(dname) as fin:
    for line in fin:
        if not line[0] == '#':
            cont = line.strip().split(',')
            sample_group[cont[0]] = cont[5]


TOTAL_YES = np.sum(sample_group)
TOTAL_NO = len(sample_group) - TOTAL_YES
print(TOTAL_YES, TOTAL_NO)

def getSampleName(name):
    return name[name.rfind('/')+1:]

def takeLine(line, ratio=0.6):
    cont = line.split()
    cnt_yes = 0
    cnt_no = 0
    for i, v in enumerate(cont[1:]):
        if i%2==0:
            if CB_MAPPING[sample_group[getSampleName(v)]]:
                cnt_yes += 1
            else:
                cnt_no += 1
    # if cnt_yes >= ratio*TOTAL_YES or cnt_no >= ratio*TOTAL_NO:
    #     return True
    # return False
    print(cnt_yes, cnt_no)
    return True

def getLine(line):
    cont = line.split()
    tmp = {"sample":cont[0]}
    for i in range(len(cont[1:])//2):
        tmp[getSampleName(cont[2*i+1])] = cont[2*i+2]
    return tmp

data = []
with open(inName, 'r') as fin:
    for line in fin:
        if takeLine(line, ratio):
            data.append(getLine(line))

pd.DataFrame(data).fillna(0).to_csv(f"{outname}.count.tsv", sep='\t', index=None)
