class reads:
    def __init__(this, filename):
        this.readsDict = {}
        Reads_length_file = open("deduplicated_sequel-demultiplex." + filename + "ccs.fastq.gz_reads_length.json", "r")
        line = Reads_length_file.readline()
        readNameBool = False
        lengthBool = False
        name = ""
        length = ""
        for c in line:
            if c == '\"':
                readNameBool = not(readNameBool)
            elif c == ',':
                lengthBool = False
                readsDict.update({name:int(length)})
                name = ""
                length = ""
            elif c == ':':
                lengthBool = True
            elif c == ' ':
                continue
            elif readNameBool:
                name = name + c
            elif lengthBool:
                length = length + c
    def getLength(this, name):
        return readsDict[name]