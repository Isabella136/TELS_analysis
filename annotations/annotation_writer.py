from bs4 import BeautifulSoup
from Bio import SeqIO
import requests
import csv

with open('TELS2_MGEs_Annotations.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    reader.__next__()
    class_and_db = dict()
    for row in reader:
        class_and_db[row[0]] = [row[1], row[2]]

gene_count = 0

with open('MGEs_Classification/FASTA/iceberg_db.fasta', 'r') as fastafile:
    iceberg = dict()
    for ref in SeqIO.parse(fastafile, 'fasta'):
        if ref.id not in class_and_db:
            continue
        iceberg[str(ref.id)] = [ref.id.split('|')[2], 
                                ref.description[ref.description.find(' ')+1:(
                                                ref.description.find(',') if ',' in ref.description
                                                else ref.description.rfind('.'))]]
        
        gene_count += 1
        print(str(gene_count) + ' / ' + str(len(class_and_db)))
        
with open('MGEs_Classification/FASTA/plasmid_finder_db.fasta', 'r') as fastafile:
    plasmid = dict()
    for ref in SeqIO.parse(fastafile, 'fasta'):
        if ref.id not in class_and_db:
            continue
        plasmid[str(ref.id)] = list()

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
            plasmid[gene].append(line[11:(
                line.find(',') if ',' in ref.description
                else len(line))])
            break

    gene_count += 1
    print(str(gene_count) + ' / ' + str(len(class_and_db)))

with open('iceberg_annotations.csv', 'w') as csv_file:
    for accession in iceberg:
        csv_file.write(accession)
        csv_file.write(',' + class_and_db[accession][0])
        for annot in iceberg[accession]:
            csv_file.write(',' + annot)
        csv_file.write('\n')

with open('plasmid_annotations.csv', 'w') as csv_file:
    for accession in plasmid:
        csv_file.write(accession)
        csv_file.write(',' + class_and_db[accession][0])
        for annot in plasmid[accession]:
            csv_file.write(',' + annot)
        csv_file.write('\n')

with open('MGEs_Classification/FASTA/aclame_db.fasta', 'r') as fastafile:
    aclame = dict()
    for ref in SeqIO.parse(fastafile, 'fasta'):
        if ref.id not in class_and_db:
            continue
        aclame[str(ref.id)] = list()
        if 'MgeName: ' in ref.description:
            group = ref.description.split(': ')[5]
            aclame[ref.id].append(group[:group.find('#')-1])
        elif 'genbank:GeneID:' in ref.description:
            ID = ref.description.split(':')[-1]
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
                        aclame[ref.id].append(line)
                        break
        gene_count += 1
        print(str(gene_count) + ' / ' + str(len(class_and_db)))

with open('aclame_annotations.csv', 'w') as csv_file:
    for accession in aclame:
        csv_file.write(accession)
        csv_file.write(',' + class_and_db[accession][0])
        for annot in aclame[accession]:
            csv_file.write(',' + annot)
        csv_file.write('\n')