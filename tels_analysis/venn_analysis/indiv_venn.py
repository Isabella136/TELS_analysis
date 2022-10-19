from matplotlib import pyplot
from matplotlib_venn import venn3, venn2, venn3_circles, venn2_circles
import numpy, os

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
        lineNum = 0
        for line in sampleFile:
            lineNum += 1
            if lineNum < 20:
                continue
            line_list = line.split(',')
            gene_header = line_list[0].split('|')
            if gene_header[4] not in this.groupCount:
                this.groupCount[gene_header[4]]=[False,False,False,False]
            this.groupCount[gene_header[4]][index-1] = True

    def findFinalCount(this):
        for geneCount in this.groupCount.values():
            trueIndex = []
            for i in range(0,4):
                if geneCount[i]:
                    trueIndex.append(i+1)
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
        fig, axs = pyplot.subplot_mosaic([['l','ur'],['l','cr'],['l','lr']], gridspec_kw=gs_kw, constrained_layout=True)
        fig.suptitle(this.sample, fontsize=20)

        pyplot.sca(axs['l'])
        data = makeVenn3()
        v = venn3(subsets=(data['1'], data['2'], data['12'], data['3'], data['13'], data['23'], data['123']), set_labels=('','','') ,set_colors=('lightcoral', 'seagreen', 'deepskyblue'), alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
        c = venn3_circles(subsets=(data['1'], data['2'], data['12'], data['3'], data['13'], data['23'], data['123']), linestyle="solid", linewidth=.5, color='black')
        if this.sample == "Bovine+XT":
            v.get_label_by_id('001').set_text('')
            pyplot.annotate(data['3'], xy=v.get_label_by_id('001').get_position(), xytext=(-0.36, -0.53), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('101').set_text('')
            pyplot.annotate(data['13'], xy=v.get_label_by_id('101').get_position(), xytext=(-0.38, -0.28), ha='center', fontsize=15)
            v.get_label_by_id('110').set_text('')
            pyplot.annotate(data['12'], xy=v.get_label_by_id('110').get_position(), xytext=(-0.17, 0.1), ha='center', fontsize=15)
            v.get_label_by_id('011').set_text('')
            pyplot.annotate(data['23'], xy=v.get_label_by_id('011').get_position() - numpy.array([.05,.02]), xytext=(-0.22, -0.55), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.22, -0.22), ha='center', fontsize=15)
        elif this.sample == "Human+V2":
            v.get_label_by_id('101').set_text('')
            pyplot.annotate(data['13'], xy=v.get_label_by_id('101').get_position() - numpy.array([0,.02]), xytext=(-0.47, -0.36), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.25, -0.12), ha='center', fontsize=15)
            v.get_label_by_id('011').set_text('')
            pyplot.annotate(data['23'], xy=v.get_label_by_id('011').get_position(), xytext=(-0.13, -0.27), ha='center', fontsize=15)
        elif this.sample == "Human+XT":
            v.get_label_by_id('101').set_text('')
            pyplot.annotate(data['13'], xy=v.get_label_by_id('101').get_position(), xytext=(.16, -.36), ha='center', fontsize=15)
            v.get_label_by_id('100').set_text('')
            pyplot.annotate(data['1'], xy=v.get_label_by_id('100').get_position(), xytext=(-.35, 0), ha='center', fontsize=15)
            v.get_label_by_id('110').set_text('')
            pyplot.annotate(data['12'], xy=v.get_label_by_id('110').get_position(), xytext=(.16, .15), ha='center', fontsize=15)
        elif this.sample == "Mock+V2":
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.475, -0.015), ha='center', fontsize=15)
            v.get_label_by_id('001').set_text('')
            pyplot.annotate(data['3'], xy=v.get_label_by_id('001').get_position(), xytext=(-0.45, -0.5), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('010').set_text('')
            pyplot.annotate(data['2'], xy=v.get_label_by_id('010').get_position() + numpy.array([.015,.01]), xytext=(0.52, 0.47), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('011').set_text('')
            pyplot.annotate(data['23'], xy=v.get_label_by_id('011').get_position(), xytext=(0, 0), ha='center', fontsize=15)
        elif this.sample == "Mock+XT":
            v.get_label_by_id('100').set_text('')
            pyplot.annotate(data['1'], xy=v.get_label_by_id('100').get_position() - numpy.array([.05,0]), xytext=(-0.65, 0), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('110').set_text('')
            pyplot.annotate(data['12'], xy=v.get_label_by_id('110').get_position(), xytext=(-.05, 0.5), ha='center', fontsize=15)
            v.get_label_by_id('010').set_text('')
            pyplot.annotate(data['2'], xy=v.get_label_by_id('010').get_position() + numpy.array([.075,0]), xytext=(0.65, 0.01), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(0, 0), ha='center', fontsize=15)
        elif this.sample == "Soil+V2":
            v.get_label_by_id('111').set_text('')
            pyplot.annotate(data['123'], xy=v.get_label_by_id('111').get_position(), xytext=(-0.16, 0.35), ha='center', fontsize=15)
            v.get_label_by_id('011').set_text('')
            pyplot.annotate(data['23'], xy=v.get_label_by_id('011').get_position(), xytext=(0.05, 0.25), ha='center', fontsize=15)
            v.get_label_by_id('010').set_text('')
            pyplot.annotate(data['2'], xy=v.get_label_by_id('010').get_position(), xytext=(0.05, 0.5), ha='center', fontsize=15)




        pyplot.sca(axs['ur'])
        data = makeVenn2('1')
        v = venn2(subsets=(data['1'], data['4'], data['14']), set_labels=('',''), set_colors=('lightcoral', 'gold'), alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
        c = venn2_circles(subsets=(data['1'], data['4'], data['14']), linestyle="solid", linewidth=.5, color='black')
        if this.sample == "Bovine+XT":
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position(), xytext=(0.77, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Human+V2":
            v.get_label_by_id('10').set_text('')
            pyplot.annotate(data['1'], xy=v.get_label_by_id('10').get_position(), xytext=(-0.78, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Human+XT":
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position() - numpy.array([.01,0]), xytext=(0.75, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Mock+V2":
            v.get_label_by_id('10').set_text('')
            v.get_label_by_id('11').set_text('')
            pyplot.annotate(data['14'], xy=v.get_label_by_id('11').get_position() - numpy.array([.03,0]), xytext=(-0.77, -.05), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Mock+XT":
            v.get_label_by_id('10').set_text('')
            pyplot.annotate(data['1'], xy=v.get_label_by_id('10').get_position() + numpy.array([.01,0]), xytext=(-0.8, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position() - numpy.array([.01,0]), xytext=(0.75, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Soil+V2":
            v.get_label_by_id('10').set_text('')
            v.get_label_by_id('11').set_text('')
            pyplot.annotate(data['14'], xy=v.get_label_by_id('11').get_position() - numpy.array([.05,0]), xytext=(-0.76, -.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Soil+XT":
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position() + numpy.array([.025,0]), xytext=(0.8, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))



        pyplot.sca(axs['cr'])
        data = makeVenn2('2')
        v = venn2(subsets=(data['2'], data['4'], data['24']), set_labels=('',''), set_colors=('seagreen', 'gold'), alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
        c = venn2_circles(subsets=(data['2'], data['4'], data['24']), linestyle="solid", linewidth=.5, color='black')
        if this.sample == "Bovine+XT":
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position(), xytext=(0.8, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Mock+V2":
            v.get_label_by_id('10').set_text('')
            pyplot.annotate(data['2'], xy=v.get_label_by_id('10').get_position() + numpy.array([.03,0]), xytext=(-0.85, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position() - numpy.array([.03,0]), xytext=(0.8, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Mock+XT":
            v.get_label_by_id('10').set_text('')
            pyplot.annotate(data['2'], xy=v.get_label_by_id('10').get_position(), xytext=(-0.81, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position(), xytext=(0.75, -0.06), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        

        pyplot.sca(axs['lr'])
        data = makeVenn2('3')
        v = venn2(subsets=(data['3'], data['4'], data['34']), set_labels=('',''), set_colors=('deepskyblue', 'gold'), alpha=.5)
        for text in v.subset_labels:
            if text != None:
                text.set_fontsize(15)
        c = venn2_circles(subsets=(data['3'], data['4'], data['34']), linestyle="solid", linewidth=.5, color='black')
        if this.sample == "Mock+V2":
            v.get_label_by_id('10').set_text('')
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position(), xytext=(0.82, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Mock+XT":
            v.get_label_by_id('10').set_text('')
            pyplot.annotate(data['3'], xy=v.get_label_by_id('10').get_position() + numpy.array([.01,0]), xytext=(-0.75, -0.055), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
            v.get_label_by_id('01').set_text('')
            pyplot.annotate(data['4'], xy=v.get_label_by_id('01').get_position() + numpy.array([.01,0]), xytext=(0.81, -0.057), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))
        elif this.sample == "Soil+V2":
            v.get_label_by_id('11').set_text('')
            pyplot.annotate(data['34'], xy=v.get_label_by_id('11').get_position() + numpy.array([.03,.23]), xytext=(.3, .5), ha='center', fontsize=15, arrowprops=dict(arrowstyle='-', color='black'))


        pyplot.gcf()
        if not(os.path.exists(outputFolder)):
            os.makedirs(outputFolder)
        pyplot.savefig(outputFolder + "/arg_" + this.sample + VENN)
        pyplot.close()
