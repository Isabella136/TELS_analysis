def fileDict(fileName):
    sample = ""
    seqPlatform = ""

    #Determine organism
    if fileName[0] == 'B':
        sample = "Bovine fecal"
        fileName = fileName[2:]
    elif fileName[0] == 'H':
        sample = "Human fecal"
        fileName = fileName[2:]
    elif fileName[0] == 'M':
        sample = "Mock"
        fileName = fileName[2:]
    else: #fileName[0] = 'S'
        sample = "Soil"
        fileName = fileName[1:]

    #Determine whether V2 or XT
    if fileName[0] == 'V':
        sample = sample + " (V2)"
        fileName = fileName[2:]
    else: #fileName[0] = 'X'
        sample = sample + " (XT)"
        fileName = fileName[2:]

    #Determine probes
    if fileName[0:2] == "AM":
        sample = sample + " + ARG-MGE probe"
        seqPlatform = "TELSeq"
        fileName = fileName[2:]
    elif fileName[0] == 'A':
        sample = sample + " + ARG probe"
        seqPlatform = "TELSeq"
        fileName = fileName[1:]
    elif fileName[0] == 'M':
        sample = sample + " + MGE probe"
        seqPlatform = "TELSeq"
        fileName = fileName[1:]
    else: #fileName[0] == 'N'
        seqPlatform = "PacBio"
        return (sample, seqPlatform)

    #Determine baits
    sample = sample + fileName
    return (sample, seqPlatform)

def x_axis(fileName):
    sub_table = ""
    x_axis = ""

    #Determine organism
    if fileName[0] == 'B':
        sub_table = "Bovine"
        fileName = fileName[2:]
    elif fileName[0] == 'H':
        sub_table = "Human"
        fileName = fileName[2:]
    elif fileName[0] == 'M':
        sub_table = "Mock"
        fileName = fileName[2:]
    else: #fileName[0] = 'S'
        sub_table = "Soil"
        fileName = fileName[1:]

    #Determine whether V2 or XT
    if fileName[0] == 'V':
        x_axis = "+V2"
        fileName = fileName[2:]
    else: #fileName[0] = 'X'
        x_axis = "+XT"
        fileName = fileName[2:]

    #Determine probes
    if fileName[0:2] == "AM":
        x_axis = "TELSeq" + x_axis + "+ARG-MGE probes"
    elif fileName[0] == 'A':
        x_axis = "TELSeq" + x_axis + "+ARG probes"
    elif fileName[0] == 'M':
        x_axis = "TELSeq" + x_axis + "+MGE probes"
    else: #fileName[0] == 'N'
        x_axis = "PacBio"
    return (sub_table,x_axis)

def megares_analyzer(megaresFile):
    def sortList(mechanismsList):
        if len(mechanismsList) == 1:
            return mechanismsList
        else:
            listA = sortList(mechanismsList[0:int(len(mechanismsList)/2)])
            listB = sortList(mechanismsList[int(len(mechanismsList)/2):])
            toReturn = []
            for i in range(0,len(mechanismsList)):
                if len(listA) == 0:
                    toReturn.append(listB.pop(0))
                elif len(listB) == 0:
                    toReturn.append(listA.pop(0))
                elif listA[0][0] < listB[0][0]:
                    toReturn.append(listA.pop(0))
                else:
                    toReturn.append(listB.pop(0))
            return toReturn
    drugList = []
    otherList = []
    megares = open(megaresFile, "r")
    megares.readline()
    for line in megares:
        splitLine = line.split('|')
        tempTuple = (splitLine[2], splitLine[3])
        if splitLine[2] == "betalactams": tempTuple = ("Betalactams", splitLine[3])
        if splitLine[1] == "Drugs":
            if drugList.count(tempTuple) == 0:
                drugList.append(tempTuple)
        else:
            if otherList.count(tempTuple) == 0:
                otherList.append(tempTuple)
    megares.close()
    drugList = sortList(drugList)
    otherList = sortList(otherList)
    return (drugList, otherList)

def getGenesLength(megaresFile):
    megares = open(megaresFile, "r")
    toReturn = {}
    header = False
    name = ""
    for line in megares:
        header = not(header)
        if header: 
            name = line.split('|')[0][1:]
            continue
        toReturn.update({name:len(line)})
    megares.close()
    return toReturn

def getSampleGroup(fileName):
    sub_table = ""
    legend = ""

    #Determine organism
    if fileName[0] == 'B':
        sub_table = "Bovine"
        fileName = fileName[2:]
    elif fileName[0] == 'H':
        sub_table = "Human"
        fileName = fileName[2:]
    elif fileName[0] == 'M':
        sub_table = "Mock"
        fileName = fileName[2:]
    else: #fileName[0] = 'S'
        sub_table = "Soil"
        fileName = fileName[1:]

    #Determine whether V2 or XT
    if fileName[0] == 'V':
        sub_table = sub_table + "+V2"
        fileName = fileName[2:]
    else: #fileName[0] = 'X'
        sub_table = sub_table + "+XT"
        fileName = fileName[2:]

    #Determine probes
    if fileName[0:2] == "AM":
        legend = "TELSeq+ARG-MGE probes"
    elif fileName[0] == 'A':
        legend = "TELSeq+ARG probes"
    elif fileName[0] == 'M':
        legend = "TELSeq+MGE probes"
    else: #fileName[0] == 'N'
        legend = "PacBio"
    return (sub_table,legend)

def getSampleAndIndex(fileName):
    sub_table = ""
    index = 0

    #Determine organism
    if fileName[0] == 'B':
        sub_table = "Bovine"
        fileName = fileName[2:]
    elif fileName[0] == 'H':
        sub_table = "Human"
        fileName = fileName[2:]
    elif fileName[0] == 'M':
        sub_table = "Mock"
        fileName = fileName[2:]
    else: #fileName[0] = 'S'
        sub_table = "Soil"
        fileName = fileName[1:]

    #Determine whether V2 or XT
    if fileName[0] == 'V':
        sub_table = sub_table + "+V2"
        fileName = fileName[2:]
    else: #fileName[0] = 'X'
        sub_table = sub_table + "+XT"
        fileName = fileName[2:]

    #Determine probes
    if fileName[0:2] == "AM":
        index = 2
    elif fileName[0] == 'A':
        index = 1
    elif fileName[0] == 'M':
        index = 3
    else: #fileName[0] == 'N'
        index = 4
    return (sub_table,index)

def mgeDict(filePath):
    mgeDict = {}
    mgeFile = open(filePath, "r")
    for line in mgeFile:
        mgeDict[line.split(",")[0]] = line.split(",")[1][:-1]
    mgeFile.close()
    return mgeDict