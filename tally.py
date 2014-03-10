import glob
import sys
import os
# to use:
# nohup python tally.py [file.sam] [chromosome]
arg=sys.argv
FILE = arg[1]


CHR=arg[2]


def readHG():
	f = open('mm10.hg','r')
	hg = f.read().split('\n')[2:-2]
	f.close()	
	chromosomes = {}
	for line in hg:
		parts =  line.split('\t')
		CHR = parts[1][3:]
		start = int(parts[2])
		end = int(parts[3])
		type = parts[-2]
		#print(CHR + '\t' + str(start) + '\t' + str(end) + '\t' + type)
		if not(CHR in chromosomes): 
			chromosomes.update({CHR:{'telomere':[],'centromere':[]}})
		chromosomes[CHR][type].append(start)
		chromosomes[CHR][type].append(end)
	for CHR in chromosomes:
		chromosomes[CHR]['telomere']=sorted(chromosomes[CHR]['telomere'])[1:-1]
	return chromosomes

def findClosest(what,CHR,where,hg):
	locations = hg[CHR][what]
	distances = []
	if locations == []:
		return 'none'
	for d in locations:
		distances.append(abs(d-where))
	dist=str(sorted(distances)[0])
	return dist


def getGTFline(line):
    details = line.split('\t')
    name=details[-1].split('gene_name "')[-1]
    name = name[:name.find('"')]
    chromosome =details[0]
    TYPE =details[1]
    dictionary={'chr':chromosome,'type':TYPE,'start':details[3],'end':details[4],'strand':details[6],'name':name}
    return dictionary

def makeGTFdict(FILE,chromosome, allowedTypes = ['protein_coding'], exons = False):
    hg = readHG()
    f= open(FILE,'r')
    bothDicts={'+':{},'-':{},'both':{}}
    n=0
    while True:
        line = f.readline()
        if line == '':
            break
        line = getGTFline(line)
        if line['chr']==chromosome and line['type'] in allowedTypes:
            if exons:
                if line['name'] in bothDicts[line['strand']]:
                    bothDicts[line['strand']][line['name']].append((int(line['start']),int(line['end'])))
                else:
                    bothDicts[line['strand']].update({line['name']:[(int(line['start']),int(line['end']))]})
            else:
                if line['name'] in bothDicts[line['strand']]:
                    if int(line['start']) < bothDicts[line['strand']][line['name']][0]:
                        bothDicts[line['strand']][line['name']][0]=int(line['start'])
                    if int(line['end']) > bothDicts[line['strand']][line['name']][1]:
                        bothDicts[line['strand']][line['name']][1]=int(line['end'])
                else:
                    bothDicts[line['strand']].update({line['name']:[int(line['start']),int(line['end'])]})
        n+=1
    f.close()
    x=bothDicts['+'].copy()
    for gene in x:
        x[gene].append('+')
    y=bothDicts['-'].copy()
    for gene in y:
        y[gene].append('-')
    bothDicts['both'].update(x)
    bothDicts['both'].update(y)
    for i in bothDicts['both']:
        mid = sum(bothDicts['both'][i][:-1])/2
        bothDicts['both'][i].append(findClosest('telomere',chromosome,mid,hg))
	bothDicts['both'][i].append(findClosest('centromere',chromosome,mid,hg))
    return bothDicts

def bin(n):
    digs = []
    s = ''
    if n<0:
        s = '-'
        n = -n
    while True:
        digs.append(str(n%2))
        n /= 2
        if not n: break
    if s: digs.append(s)
    digs.reverse()
    return ''.join(digs)

def readSAMline(read):
    line = read.split('\t')
    CHR =line[2]
    start = int(line[3])
    end = start + len(line[9])
    flag = int(line[1])
    flag = bin(flag)
    flag = str(flag)[-5]
    if flag == '1': strand = '-'
    else: strand = '+'
    return {'chr' : CHR,'start':start,'end':end,'strand':strand}

def Tally(FILE,chromosome,GFTdict):
    f= open(FILE,'r')
    while f.readline()[0] == '@':
	pass
    tally=dict.fromkeys(GFTdict['both'].keys())
    for gene in tally:
        tally[gene] = {'sense':0,'anti':0,'telomere':GFTdict['both'][gene][3],'centromere':GFTdict['both'][gene][4]}
    n=0
    while True:
        line = f.readline()
        if line == '':
            break
        line = readSAMline(line)
        if line['chr']==chromosome:
            for gene in GFTdict['both']:
                if line['start']>=GFTdict['both'][gene][0]:
                    if line['end']<=GFTdict['both'][gene][1]:
                        if not(gene in tally):
                            tally.update({gene:{'sense':0,'anti':0,'telomere':GFTdict['both'][gene][3],'centromere':GFTdict['both'][gene][4]}})
                        if gene in tally:
                            if line['strand']==GFTdict['both'][gene][2]:
                                tally[gene]['sense']+=1
                            else:
                                tally[gene]['anti']+=1
                        #print(gene+ str(tally[gene]))


            if n%50000==0:print(n)
        if n%1000000==0:print(n)
        n+=1
    print('scanned all bed')
    f.close()
    return tally


gtf = glob.glob('*.gtf')[0]

print('reading gtf')

biggusDicus = makeGTFdict(gtf,CHR)
#biggusDicus[strand][gene][exon]
#biggusDicus['-']['Hcrtr1'][0]

#gene =biggusDicus['-']['Cdkn2a']

print('scanning sam file')

results = Tally(FILE,CHR,biggusDicus)
print('saving table')
table=[]
for gene in results:

    table.append('\t'.join([gene,str(results[gene]['sense']),str(results[gene]['anti']),str(results[gene]['telomere']),str(results[gene]['centromere'])]))
table = '\n'.join(table)
header = 'gene\tsense\tantisense\ttelomere\tcentromere\n'
f=open(FILE[:-3]+'_'+CHR+'_'+'.csv','w')
f.write(header+table)
f.close()





