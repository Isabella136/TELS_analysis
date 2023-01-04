import requests
from bs4 import BeautifulSoup

aclame = {}
iceberg = {}
plasmid = {}
all_mges = open('../unique_mges.csv', 'r')
for line in all_mges:
    gene, type = line.split(',')
    if type[-1] == '\n':
        type = type[:-1]
    if gene == 'GENE':
        continue
    if 'gene:' in gene:
        aclame.update({gene:[type]})
    elif 'ICEberg' in gene:
        iceberg.update({gene:[type]})
    else:
        plasmid.update({gene:[type]})
all_mges.close()

for gene in iceberg:
    header = gene.split('|')
    iceberg[gene].append(header[2])
    URL = 'https://www.ncbi.nlm.nih.gov/nuccore/' + header[-2]
    r = requests.get(URL)
    soup = BeautifulSoup(r.content, 'xml')
    xml = soup.prettify().split('\n')
    i = 0
    for line in xml:
        i+=1
        if '<div class="rprtheader">' in xml[i-3]:
            iceberg[gene].append(line)
            break
    print(gene)

for gene in plasmid:
    accession = gene[gene.rfind('_')+1:]
    URL = 'https://www.ncbi.nlm.nih.gov/nuccore/' + accession
    r = requests.get(URL)
    soup = BeautifulSoup(r.content, 'xml')
    xml = soup.prettify().split('\n')
    i = 0
    for line in xml:
        i += 1
        if '<div class="rprtheader">' in xml[i-3]:
            iceberg[gene].append(line)
            break
    print(gene)

aclame_db = open('aclame_db.fasta', 'r')
for line in aclame_db:
    if line[0] != '>':
        continue
    line = line[1:]
    header = line.split(' # ')
    if header[0][:header[0].find(' ')] in aclame:
        for info in header:
            if 'MgeName: ' in info:
                group = info.split(': ')[1]
                aclame[header[0][:header[0].find(' ')]].append(group)
            if 'genbank:GeneID:' in info:
                ID = info.split(':')[-1]
                if ID[-1] == '\n':
                    ID = ID[:-1]
                URL = 'https://www.ncbi.nlm.nih.gov/gene/' + ID
                r = requests.get(URL)
                soup = BeautifulSoup(r.content, 'xml')
                xml = soup.prettify().split('\n')
                i = 0
                inDiv = False
                for line in xml:
                    i += 1
                    if 'id="summaryDiv"' in line:
                        inDiv = True
                    if inDiv:
                        if "</div>" in line:
                            inDiv = False
                        if 'Organism' in xml[i-5]:
                            aclame[header[0][:header[0].find(' ')]].append(line)
                            break
        print(header[0][:header[0].find(' ')])
aclame_db.close()

temp = open('temp.csv', 'w')
for gene, infoList in aclame.items():
    temp.write(gene)
    for info in infoList:
        temp.write(',' + info)
    temp.write('\n')
for gene, infoList in iceberg.items():
    temp.write(gene)
    for info in infoList:
        temp.write(',' + info)
    temp.write('\n')
for gene, infoList in plasmid.items():
    temp.write(gene)
    for info in infoList:
        temp.write(',' + info)
    temp.write('\n')
temp.close()
