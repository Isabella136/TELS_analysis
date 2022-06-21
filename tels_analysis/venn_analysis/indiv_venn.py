from matplotlib import pyplot
from matplotlib_venn import venn3, venn2

class indiv_venn:
    def __init__(this, sample):
        this.data = {
            '1':0,      '2':0,      '3':0,      '4':0,      '123':0,
            '12':0,     '13':0,     '14':0,     '23':0,     '24':0,     '34':0,
            } #where 1 is arg, 2 is arg-mge, 3 is mge, 4 is pacbio
        this.groupCount = {}
        this.sample = sample

    def addToCount(this, filepath, index):
        sampleFile = open(filepath, "r")
        for line in sampleFile:
            blank_space = 0
            gene_header_list = line.split(',')[1:]
            for gene in gene_header_list:
                if gene == "":
                    blank_space += 1
                    if blank_space == 2: break
                    continue
                gene_header = gene.split('|')
                if gene_header[4] not in this.groupCount:
                    this.groupCount[gene_header]=[False,False,False,False]
                this.groupCount[gene_header][index-1] = True

    def findFinalCount(this):
        for geneCount in this.groupCount.values():
            trueIndex = []
            for i in range(0,4):
                if geneCount[i]:
                    trueIndex.append(i)
            if len(trueIndex) == 4:
                for intersection in this.data:
                    this.data[intersection] += 1
            elif len(trueIndex) > 0:
                this.data[str(trueIndex[0])] += 1
                if len(trueIndex) > 1:
                    this.data[str(trueIndex[1])] += 1
                    this.data[str(trueIndex[0])+str(trueIndex[1])] += 1
                    if len(trueIndex) > 2:
                        this.data[str(trueIndex[2])] += 1
                        this.data[str(trueIndex[0])+str(trueIndex[2])] += 1
                        this.data[str(trueIndex[1])+str(trueIndex[2])] += 1
                        if trueIndex[2] != 4:
                            this.data[str(trueIndex[0])+str(trueIndex[1])+str(trueIndex[2])] += 1

    def makeFigure(this, outputFolder, VENN):
        def makeVenn3():
            dataForVenn = {
                '1':0,      '2':0,      '3':0,      '123':0,
                '12':0,     '13':0,     '23':0
                }
            dataForVenn['123'] = this.data['123']
            dataForVenn['12'] = this.data['12'] - dataForVenn['123']
            dataForVenn['13'] = this.data['13'] - dataForVenn['123']
            dataForVenn['23'] = this.data['23'] - dataForVenn['123']
            dataForVenn['1'] = this.data['1'] - dataForVenn['12'] - dataForVenn['13'] - dataForVenn['123']
            dataForVenn['2'] = this.data['2'] - dataForVenn['12'] - dataForVenn['23'] - dataForVenn['123']
            dataForVenn['3'] = this.data['3'] - dataForVenn['13'] - dataForVenn['23'] - dataForVenn['123']
            return dataForVenn
        def makeVenn2(index):
            dataForVenn = {
                index:0,    '4':0,       index+'4':0
                }
            dataForVenn[index+'4'] = this.data[index+'4']
            dataForVenn[index] = this.data[index] - dataForVenn[index+'4']
            dataForVenn['4'] = this.data['4'] - dataForVenn[index+'4']
            return dataForVenn
        gs_kw = dict(width_ratios=[3,1],height_ratios=[1,1,1])
        fig, axs = pyplot.subplot.mosaic([['l','ur'],['l','cr'],['l','lr']], gridspec_kw=gs_kw, constrained_layout=True)
        fig.subtitle(this.sample)

        pyplot.sca(axs['l'])
        data = makeVenn3()
        venn3(subsets=(data['1'], data['2'], data['12'], data['3'], data['13'], data['23'], data['123']), set_labels=("ARG probes", "ARG-MGE proges", "MGE probes"))

        pyplot.sca(axs['ur'])
        data = makeVenn2('1')
        venn2(subsets=(data['1'], data['4'], data['14']), set_labels=("ARG probes", "PacBio"))

        pyplot.sca(axs['cr'])
        data = makeVenn2('2')
        venn2(subsets=(data['2'], data['4'], data['24']), set_labels=("ARG-MGE probes", "PacBio"))

        pyplot.sca(axs['lr'])
        data = makeVenn2('3')
        venn2(subsets=(data['3'], data['4'], data['34']), set_labels=("MGE probes", "PacBio"))

        pyplot.gcf()
        pyplot.savefig(outputFolder + "/arg_" + this.sample + VENN)
        pyplot.close()
