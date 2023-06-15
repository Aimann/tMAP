#!/usr/bin/env python3
#Aidan Manning 2.21.23
#Script to determine sites of tRNA modficiation from sequencing data and generates some plots
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import re
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
import seaborn as sns

class tMAPalign(object):

    def __init__(self):

        self.eukpositions = list(['drop1',-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e5','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])
        self.eukpositionD = {q:k for q,k in enumerate(self.eukpositions[1:-1])}
        self.fixSprinzl = {'-1':'-1','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6','7':'7','8':'8','9':'9','10':'10','11':'11','12':'12','13':'13','14':'14','15':'15','16':'16','17':'17','17a':'17a','18':'18','19':'19','20':'20','20a':'20a','20b':'20b','21':'21','22':'22','23':'23','24':'24','25':'25','26':'26','27':'27','28':'28','29':'29','30':'30','31':'31','32':'32','33':'33','34':'34','35':'35','36':'36','37':'37','38':'38','39':'39','40':'40','41':'41','42':'42','43':'43','44':'44','45':'45','e1':'e11','e2':'e12','e3':'e13','e4':'e14','e5':'e15','e6':'e16','e7':'e17','e8':'e1','e9':'e2','e10':'e3','e11':'e4','e12':'e5','e13':'e27','e14':'e26','e15':'e25','e16':'e24','e17':'e23','e18':'e22','e19':'e21','46':'46','47':'47','48':'48','49':'49','50':'50','51':'51','52':'52','53':'53','54':'54','55':'55','56':'56','57':'57','58':'58','59':'59','60':'60','61':'61','62':'62','63':'63','64':'64','65':'65','66':'66','67':'67','68':'68','69':'69','70':'70','71':'71','72':'72','73':'73','74':'74','75':'75','76':'76'}

    def readSprinzl(self):
        
        if args.org == 'bact':
            positions = list(['drop1',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])
        elif args.org == 'arch':
            positions = list(['drop1',-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])
        elif args.org == 'euk':
            positions = list(['drop1',-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e5','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])
        elif args.org == 'mito':
            positions = list(['drop1',-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])

        positionD = {q:k for q,k in enumerate(positions[1:-1])}

        return positionD, positions

    def readCoverage(self, df):
        data = pd.read_csv(df, sep='\t')
        # data = data[data['Sample'].isin(list(self.conditionDict.keys()))]
        data = data[data['Feature'].str.startswith('tRNA')]
        # data = data['tRNAreadstotal'].astype()
        data = data[data['tRNAreadstotal'] >= int(args.minreads)]
        data['position'] = data['position'].map(self.fixSprinzl)
        data['Isoacceptor'] = data['Feature'].str.replace('(-[0-9]+)', '', regex=True)
        isodecoders = list(data['Feature'].drop_duplicates())
        isoacceptors = list(data['Isoacceptor'].drop_duplicates())

        return data, isodecoders, isoacceptors

    def isodecoderAlign(self, trna, alignments, positions):
        isoacceptor = re.sub('(-[0-9]+)', '', trna)
        isodecoderD = {}
        isodecoderdifferencelist = []

        for line in open(alignments):
            if line.startswith(isoacceptor):
                isodecoders = line.strip().split()
                sequence = isodecoders[1].replace('.', '')
                sequence = re.sub("[a-z]", '', sequence)
                isodecoderD[isodecoders[0]] = sequence

        isodecoderdf = pd.DataFrame.from_dict(isodecoderD, orient='index')
        isodecoderdf[positions] = isodecoderdf[0].str.split('', expand=True)
        isodecoderdf = isodecoderdf.drop([0, 'drop1', 'drop2'], axis=1)
        testvector = isodecoderdf.to_numpy()
        checkfordiffs = list((testvector[0] == testvector).all(0))

        for idx,pos in enumerate(list(isodecoderdf.columns)):
            if checkfordiffs[idx] == True:
                pass
            else:
                isodecoderdifferencelist.append(pos)
                pass
        diffvalues = ','.join(str(i) for i in isodecoderdifferencelist)

        return (isoacceptor, diffvalues) , isodecoderdf

    def compareIsodecoders(self, trna, sprinzlalign):

        comparisonList = []
        isoacceptor = re.sub('(-[0-9]+)', '', trna)
        sprinzlalign['Isoacceptor'] = sprinzlalign.index.str.replace('(-[0-9]+)', '', regex=True)
        isoacceptordf = sprinzlalign[sprinzlalign['Isoacceptor'] == isoacceptor]
        list_of_isodecoders = list(isoacceptordf.index)
        list_of_isodecoders.remove(trna)
        for isod in list_of_isodecoders:
            difflist = []
            diffnuclist = []
            queryrefdf = isoacceptordf.loc[[trna, isod]].drop(['Isoacceptor'], axis=1)
            testvector = queryrefdf.to_numpy()
            checkfordiffs = list((testvector[0] == testvector).all(0))
            for idx,pos2 in enumerate(list(queryrefdf.columns)):
                if checkfordiffs[idx] == True:
                    pass
                else:
                    ref  = str(queryrefdf.loc[trna][pos2])
                    query  = str(queryrefdf.loc[isod][pos2])
                    diffnuclist.append(ref+':'+query)
                    difflist.append(str(pos2))
            diffvals = ','.join(difflist)
            diffvalnucs = ','.join(diffnuclist)
            comparisonList.append((trna, isod, diffvals, len(difflist), diffvalnucs))
            
        return comparisonList
    
    def plotVariance(self, variance, positions, opath):
        
        varianceD = {str(q):0 for q in positions[1:-1]}

        isoacceptors = list(variance.index)

        for isoac in isoacceptors:

            values = variance.loc[isoac]['Isodecoder Variation Sprinzl'].split(',')
            if values[0] != '':
                for val in values:
                    varianceD[val] += 1
            else:
                pass

        variancedf = pd.DataFrame.from_dict(varianceD, orient='index').reset_index().reset_index()

        sprinzl_x = dict(zip(list(variancedf['index']), list(variancedf['level_0'])))

        fig, axes = plt.subplots(1,1, sharex=False, figsize=(4, 1))
        axes.bar(variancedf['level_0'], variancedf[0], width=.95, linewidth=0,color='#D3D3D3')

        # axes[idx][qdx].set_title(tRNA  + ' | ' + cond, fontsize=10, loc='left', pad=2)
        xlabels = list(variancedf['index'])
        xticks = list(variancedf['level_0'])
        axes.set_xlim(-1, len(xlabels))
        axes.set_xticks(xticks)
        axes.set_xticklabels(xlabels, fontsize=3, rotation=90, horizontalalignment='center')
        for axis in ['top', 'bottom', 'left', 'right']:
            axes.spines[axis].set_linewidth(0.2)

        trnafeatlist = [sprinzl_x.get('9') + 0.3,sprinzl_x.get('21') + 0.3,sprinzl_x.get('26') + 0.3,sprinzl_x.get('38') + 0.3,sprinzl_x.get('48') + 0.3,sprinzl_x.get('60') + 0.3]
        trnawidthlist = [(sprinzl_x.get('13') + 0.5) - (sprinzl_x.get('9') + 0.4),(sprinzl_x.get('25') + 0.5) - (sprinzl_x.get('21') + 0.5),(sprinzl_x.get('31') + 0.5) - (sprinzl_x.get('26') + 0.5),(sprinzl_x.get('43') + 0.5) - (sprinzl_x.get('38') + 0.5),(sprinzl_x.get('53') + 0.5) - (sprinzl_x.get('48') + 0.5),(sprinzl_x.get('65') + 0.5) - (sprinzl_x.get('60') + 0.5)]
        ymaxval = variancedf[0].max() + 2
        featcolorD = {'dloop':'#FF9999','ac':'#99FF99','tloop':'#99FFFF'}
        for z,ss  in enumerate(['dloop', 'dloop', 'ac', 'ac', 'tloop', 'tloop']):
            axes.axvspan(trnafeatlist[z], trnafeatlist[z] + trnawidthlist[z] + 0.3, facecolor=featcolorD.get(ss), alpha=0.2, zorder=-100)
        axes.set_ylim(0, ymaxval)
        axes.set_title(args.o, fontsize=6, loc='left', pad=2)
        axes.set_yticks([0,ymaxval])
        axes.set_yticklabels(['0',str(ymaxval)], fontsize=6)
        axes.tick_params(axis='both', which='major', pad=0.5, width=0.25,length=0.25)
        for axis in ['top', 'right']:
                axes.spines[axis].set_linewidth(0)
        axes.bar_label(axes.containers[0], fontsize=3)

        plt.tight_layout()
        plt.savefig(opath + args.o + '-isodecodervariance.pdf')
        plt.close()

        return

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--cov', required=True, help='-coverage.txt file from tRAX output')
    ap.add_argument('--alignments', required=True, help='stk file of mature tRNA alignments')
    ap.add_argument('--minreads', required=False, help='read coverage cutoff for tRNAs', default=0)
    ap.add_argument('--org', required=False, help='organism of interest; euk, bact, arch, mito; default=euk', default='euk')
    ap.add_argument('--plot', required=False, help='whether or not to calculate significant mismatch differences', action='store_true')
    ap.add_argument('--o', required=True, help='path for output files')
    args = ap.parse_args()
    
    ## Generates tMAP object
    tMAPalign = tMAPalign()

    startTime = datetime.now()

    posD, positions = tMAPalign.readSprinzl()

    ## Reads in coverage file and generates a list of isoacceptors and isodecoders
    coveragedata, isodecoderlist, isoacceptorlist = tMAPalign.readCoverage(args.cov)

    ## Checks for output directory and creates if not present
    opath = os.getcwd() + '/isodecoderanalysis/'
    isdir = os.path.isdir(opath)
    if isdir == True:
        pass
    else:
        os.mkdir(opath)
    
    ## Generates a list of isodecoders and their differences
    print('Assessing tRNA isodecoder differences...')
    list_of_isodecoder_dfs = []
    list_of_diffs = []
    for isoaccept in isoacceptorlist:
        diff_tuple, diff_df = tMAPalign.isodecoderAlign(isoaccept, args.alignments, positions)
        list_of_isodecoder_dfs.append(diff_df)
        list_of_diffs.append(diff_tuple)
    isodecoderdf_all = pd.concat(list_of_isodecoder_dfs)
    isodecoderdf_all.to_csv(opath + args.o + '-Isodecoder_alignments.csv')
    isodecoder_differences = pd.DataFrame(list_of_diffs).rename(columns={0:'Isoacceptor', 1:'Isodecoder Variation Sprinzl'}).set_index('Isoacceptor')
    isodecoder_differences.to_csv(opath + args.o + '-Isodecoder_variation.csv')

    ## Removes isodecoders not present in alignment file
    isodecoderlist = [x for x in isodecoderlist if x in list(isodecoderdf_all.index)]

    ## Extracts the differences between isodecoders
    list_of_isodecoder_diffs = []
    for isod in isodecoderlist:
        list_of_isodecoder_diffs.append(tMAPalign.compareIsodecoders(isod, isodecoderdf_all))
    flat_list_of_isodecoder_diffs = [item for sublist in list_of_isodecoder_diffs for item in sublist]
    isodecodervariationdf = pd.DataFrame(flat_list_of_isodecoder_diffs, columns=['tRNA_x', 'tRNA_y', 'Sprinzl Differences', '# of Differences', 'Nucleotide Differences'])
    isodecodervariationdf_removedup = (isodecodervariationdf[~isodecodervariationdf.filter(like='tRNA').apply(frozenset, axis=1).duplicated()].reset_index(drop=True))
    isodecodervariationdf_removedup = isodecodervariationdf_removedup.set_index('tRNA_x')
    isodecodervariationdf_removedup.to_csv(opath + args.o + '-Isodecoder_comparisons.csv')

    if args.plot == True:
        ## Generates a plot of the isodecoder differences
        print('Generating plot of isodecoder differences...')
        tMAPalign.plotVariance(isodecoder_differences, positions, opath)
    else:
        pass

    print('\nTime elasped: ', datetime.now() - startTime)