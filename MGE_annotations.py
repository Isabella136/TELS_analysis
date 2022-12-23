import requests
from bs4 import BeautifulSoup

aclame = {}
iceberg = {}
plasmid = {}
all_mges = open('unique_mges.csv', 'r')
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

for gene in plasmid:
    accession = gene[gene.rfind('_')+1:]
    URL = 'https://www.ncbi.nlm.nih.gov/nuccore/' + accession
    r = requests.get(URL)
    soup = BeautifulSoup(r.content, 'xml')
    test = open("temp.txt", "w")
    test.write(soup.prettify())
    test.close()
    html = soup.prettify().split('\n')
    i = 0
    for line in html:
        i += 1
        if 'plasmid' not in line:
            continue
        if '<h1>' in html[i-2]:
            title = line.split('plasmid ')
            plasmidName= ''
            if len(title) == 1:
                title = line.split(' plasmid')
                plasmidName = title[0][title[0].rfind(' ')+1:]
            else:
                plasmidName = title[1][:title[1].find(' ')]
            if plasmidName[-1] == ',':
                plasmidName = plasmidName[:-1]
            plasmid[gene].append(plasmidName)
            break

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
                aclame[header[0]] = group
