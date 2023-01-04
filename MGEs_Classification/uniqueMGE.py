mgeList = list()
allMGEs = open("all mges.csv", "r")
uniqueMGEs = open("unique mges.xlsx", "w")
i = 0
for line in allMGEs:
    i+=1
    if (i%100 == 0):
        print(i)
    for mge in (line[:-1]).split(','):
        if (mge not in mgeList) and (mge != ""):
            mgeList.append(mge)
allMGEs.close()
for mge in mgeList:
    uniqueMGEs.write(mge + "\n")
uniqueMGEs.close()