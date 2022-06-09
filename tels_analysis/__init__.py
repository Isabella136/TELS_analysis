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
        sample = sample + " + ARG-MGE probes"
        seqPlatform = "TELSeq"
        fileName = fileName[2:]
    elif fileName[0] == 'A':
        sample = sample + " + ARG probes"
        seqPlatform = "TELSeq"
        fileName = fileName[1:]
    elif fileName[0] == 'M':
        sample = sample + " + MGE probes"
        seqPlatform = "TELSeq"
        fileName = fileName[1:]
    else: #fileName[0] == 'N'
        seqPlatform = "PacBio"
        if fileName[-1] == 'A':
            sample = sample + " + ARG probes + no baits"
        elif fileName[-2] == 'A':
            sample = sample + " + ARG-MGE probes + no baits"
        else: #fileName[-1] == 'M'
            sample = sample + " + MGE probes + no baits"
        return (sample, seqPlatform)

    #Determine baits
    sample = sample + " + bait " + fileName
    return (sample, seqPlatform)

def heatmap_x_axis(fileName):
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