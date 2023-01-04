import os
import itertools
import threading
import numpy as np
from ICEtree import ICEtree

iceberg_lock = threading.Lock()

iceberg = {}
iceberg_db = open('FASTA/iceberg_db.fasta', 'r')
header = ''
for line in iceberg_db:
    if line[0] == '>':
        header = line[1:]                   #remove '>'
        header = header.split(' ')[0]       #remove newline
        continue
    sequence = line
    if sequence[-1] == '\n':
        sequence = sequence[:-1]
    iceberg[header] = sequence
iceberg_db.close()

pls = []
iceberg_pls = open('ICEaligner.pls', 'r')
for (line, lineNum) in zip(iceberg_pls, range(46059)):
    if lineNum < 5:
        continue
    pls.append(line.split(',')[:14])
iceberg_pls.close()

output = open('ICEalignMatrix.csv', 'w')
for ice in iceberg:
    output.write(',' + ice)
output.close()

empty_doube_array = [[None]*(len(iceberg)+1)]*(len(iceberg))
arr = np.array(empty_doube_array)

threads = []

def row_alignement(index1, sub_pls):
    iceberg_lock.acquire()
    iceberg_header = list(iceberg.keys())
    iceberg_lock.release()
    for line in sub_pls:
        index2 = iceberg_header.index(line[13])+1
        if arr[index1][index2] == None:
            arr[index1][index2] = line[0]
        else:
            if int(arr[index1][index2]) < int(line[0]):
                arr[index1][index2] = line[0]

pls_index = 0
for (index1,ice1) in zip(range(0,len(iceberg)), iceberg):
    arr[index1][0] = ice1
    start = pls_index
    while (pls_index < len(pls)) and (pls[pls_index][9] == ice1):
        pls_index += 1
    end = pls_index   
    threads.append(threading.Thread(target=row_alignement, args=(index1,pls[start:end])))  

for i in range(0, 552, 12):
    for t in threads[i:i+12]:
        t.start()

    for t in threads[i:i+12]:
        t.join()

    outputMatrix = open('ICEalignMatrix.csv', 'a')
    for row in arr[i:i+12]:
        outputMatrix.write('\n')
        for cell in row:
            if cell == None:
                cell = ''
            outputMatrix.write(cell + ',')
    outputMatrix.close()

similarities = {}
iceberg_header = list(iceberg.keys())
for (row, rowNum) in zip(arr, range((len(iceberg)))):
    base = int(arr[rowNum][rowNum+1])
    similarities[iceberg_header[rowNum]] = []
    for (cell, cellNum) in zip(row[1:], range((len(iceberg)))):
        if cellNum == rowNum :
            continue
        if cell == None:
            continue
        if (base < int(arr[cellNum][cellNum+1])):
            if (int(cell)/base) >= .5:
                similarities[iceberg_header[rowNum]].append(iceberg_header[cellNum])
        else:
            if (int(cell)/int(arr[cellNum][cellNum+1])) >= .5:
                similarities[iceberg_header[rowNum]].append(iceberg_header[cellNum])

output_similarities = open('ICEalignSimilarities.csv', 'w')
for ice in similarities:
    output_similarities.write(ice)
    for similar in similarities[ice]:
        output_similarities.write(',' + similar)
    output_similarities.write('\n')
output_similarities.close()

tree = ICEtree.ICEtree(similarities)