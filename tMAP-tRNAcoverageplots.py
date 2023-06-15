#!/usr/bin/env python3
#Aidan Manning 2.21.23
#Script to determine sites of tRNA modficiation from sequencing data and generates some plots
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
import seaborn as sns

class tMAPcoverageplots(object):

    def __init__(self, sampleinfo):

        self.conditionDict = dict(zip(list(sampleinfo['sample']), list(sampleinfo['condition'])))    
        self.conditions = list(sampleinfo['condition'].drop_duplicates())
        self.fixSprinzl = {'-1':'-1','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6','7':'7','8':'8','9':'9','10':'10','11':'11','12':'12','13':'13','14':'14','15':'15','16':'16','17':'17','17a':'17a','18':'18','19':'19','20':'20','20a':'20a','20b':'20b','21':'21','22':'22','23':'23','24':'24','25':'25','26':'26','27':'27','28':'28','29':'29','30':'30','31':'31','32':'32','33':'33','34':'34','35':'35','36':'36','37':'37','38':'38','39':'39','40':'40','41':'41','42':'42','43':'43','44':'44','45':'45','e1':'e11','e2':'e12','e3':'e13','e4':'e14','e5':'e15','e6':'e16','e7':'e17','e8':'e1','e9':'e2','e10':'e3','e11':'e4','e12':'e5','e13':'e27','e14':'e26','e15':'e25','e16':'e24','e17':'e23','e18':'e22','e19':'e21','46':'46','47':'47','48':'48','49':'49','50':'50','51':'51','52':'52','53':'53','54':'54','55':'55','56':'56','57':'57','58':'58','59':'59','60':'60','61':'61','62':'62','63':'63','64':'64','65':'65','66':'66','67':'67','68':'68','69':'69','70':'70','71':'71','72':'72','73':'73','74':'74','75':'75','76':'76'}

    def readCoverage(self, df):
        data = pd.read_csv(df, sep='\t')
        data = data[data['Sample'].isin(list(self.conditionDict.keys()))]

        if args.org == 'euk':
            data = data[data['Feature'].str.startswith('tRNA')]
            data['position'] = data['position'].map(self.fixSprinzl)
        elif args.org == 'mito':
            data = data[data['Feature'].str.startswith('mt-tRNA')]
            data['position'] = data['position']
        else:
            pass

        # data = data[data['tRNAreadstotal'] >= int(args.minreads)]
        data['Isoacceptor'] = data['Feature'].str.replace('(-[0-9]+)', '', regex=True)
        data['Isotype'] = data['Isoacceptor'].str.replace('-[A-Z][A-Z][A-Z]', '', regex=True).str.lstrip('tRNA').str.lstrip('-')
        isodecoders = list(data['Feature'].drop_duplicates())
        isoacceptors = list(data['Isoacceptor'].drop_duplicates())
        isotypes = list(data['Isotype'].drop_duplicates())

        return data.dropna(), isodecoders, isoacceptors, isotypes
    
    def plotCoverage(self, data, isodecoders, isotypes, opath):

        sns.set_theme(style="white",font_scale=1)

        for iso in isotypes:
            data_iso = data[data['Isotype'] == iso]
            print("Generating coverage plots for..." + iso)
            featurestoplot = [i for i in list(data_iso['Feature'].drop_duplicates()) if i in isodecoders]
            featurestoplot = sorted(featurestoplot, key=lambda x: (x.split('-')[2], int(x.split('-')[3])))

            nucD = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}
            nucDcolors = {'adenines':'#7CAE00','cytosines':'#619CFF','thymines':'#F8766D','guanines':'#CD9600','deletions':'grey', 'ref':'#D3D3D3'}

            # nucDcolors = {'adenines':'tab:green','cytosines':'tab:blue','thymines':'tab:red','guanines':'tab:orange','deletions':'grey', 'ref':'#D3D3D3'}
            # fig, axes = plt.subplots(max(2, len(self.conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*4,max(2,len(self.conditions))*1.2))
            # fig.subplots_adjust(hspace=0.6, wspace=0.2)

            fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
            fig.subplots_adjust(hspace=0.6, wspace=0.2)

            for idx, tRNA in enumerate(featurestoplot):
                trna_data = data_iso[data_iso['Feature'] == tRNA]
                trna_data = trna_data[trna_data['actualbase'] != '-']
                trna_data = trna_data[['Feature', 'Sample', 'position', 'actualbase', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']]
                sprinzl = list(trna_data['position'].drop_duplicates())
                sprinzl_x = {sprin:featurepos for featurepos,sprin in enumerate(sprinzl)}
                trna_data['xlabel'] = trna_data['position'] + ':' + trna_data['actualbase']
                xlabels = list(trna_data['xlabel'].drop_duplicates())
                xticks = [s+0.1 for s in range(len(xlabels))]
                sprinzl_xlabel = {featurepos:sprin for featurepos,sprin in enumerate(xlabels)}
                xorder = [featurepos for featurepos,sprin in enumerate(sprinzl)]
                trna_data['Sample'] = trna_data['Sample'].map(self.conditionDict)
                trna_data_mean = trna_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).sum().reset_index()

            #     findmax = trna_data_mean['adenines'] + trna_data_mean['thymines'] + trna_data_mean['cytosines'] + trna_data_mean['guanines'] + trna_data_mean['deletions']
            #     ymaxval = round(1.1 * findmax.max())
            
                trna_data_mean['indexbase'] = trna_data_mean['position'].map(sprinzl_x)

                for qdx, cond in enumerate(self.conditions):
                    trna_data_mean_condition = trna_data_mean[trna_data_mean['Sample'] == cond]

                    findmax = trna_data_mean_condition['adenines'] + trna_data_mean_condition['thymines'] + trna_data_mean_condition['cytosines'] + trna_data_mean_condition['guanines'] + trna_data_mean_condition['deletions']
                    ymaxval = max(20,round(1.1 * findmax.max()))

                    for pos in xorder:

                        trna_data_mean_pos = trna_data_mean_condition[trna_data_mean_condition['indexbase'] == pos]
                        realbase = trna_data_mean_pos.iloc[0]['actualbase']
                        realbasecounts = float(trna_data_mean_pos[nucD.get(realbase)])
                        trna_data_mean_pos = trna_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                        conditionslice = trna_data_mean_pos[trna_data_mean_pos.columns[4:-1]].T
                        conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)

                        if float(conditionslice.max()) >= 1:
                            conditionslice = conditionslice.T
                            bottom = 0
                            for base in list(conditionslice.columns):
                                axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], linewidth=1,color=nucDcolors.get(base), bottom=bottom)
                                bottom += conditionslice.loc['y'][base]
                            axes[idx][qdx].bar(pos, realbasecounts, width=.95, linewidth=1,color='#D3D3D3', bottom=bottom)
                        else:
                            axes[idx][qdx].bar(pos, realbasecounts, width=.95, linewidth=1,color='#D3D3D3')

                    axes[idx][qdx].set_title(tRNA  + ' | ' + cond, fontsize=10, loc='left', pad=2)
                    axes[idx][qdx].set_xlim(-1, len(xlabels))
                    axes[idx][qdx].set_xticks(xticks)
                    axes[idx][qdx].set_xticklabels(xlabels, fontsize=2.5, rotation=90, horizontalalignment='center')

                    axes[idx][qdx].set_ylim(0, ymaxval)
                    axes[idx][qdx].set_yticks([0,ymaxval])
                    axes[idx][qdx].set_yticklabels(['0',str(ymaxval)], fontsize=5)
                    axes[idx][qdx].tick_params(axis='x', which='major', pad=0.5, width=0.25,length=0.25)

                    for axis in ['top', 'bottom', 'left', 'right']:
                            axes[idx][qdx].spines[axis].set_linewidth(0.1)

                    # trnafeatlist = [sprinzl_x.get('9') + 0.3,sprinzl_x.get('21') + 0.3,sprinzl_x.get('26') + 0.3,sprinzl_x.get('38') + 0.3,sprinzl_x.get('48') + 0.3,sprinzl_x.get('60') + 0.3]
                    # trnawidthlist = [(sprinzl_x.get('13') + 0.5) - (sprinzl_x.get('9') + 0.4),(sprinzl_x.get('25') + 0.5) - (sprinzl_x.get('21') + 0.5),(sprinzl_x.get('31') + 0.5) - (sprinzl_x.get('26') + 0.5),(sprinzl_x.get('43') + 0.5) - (sprinzl_x.get('38') + 0.5),(sprinzl_x.get('53') + 0.5) - (sprinzl_x.get('48') + 0.5),(sprinzl_x.get('65') + 0.5) - (sprinzl_x.get('60') + 0.5)]

                    # featcolorD = {'dloop':'#FF9999','ac':'#99FF99','tloop':'#99FFFF'}
                    # for z,ss  in enumerate(['dloop', 'dloop', 'ac', 'ac', 'tloop', 'tloop']):
                    #     axes[idx][qdx].axvspan(trnafeatlist[z], trnafeatlist[z] + trnawidthlist[z] + 0.3, facecolor=featcolorD.get(ss), alpha=0.2, zorder=-100)

                    nucDrev = {'adenines':'A','cytosines':'C','thymines':'T','guanines':'G','deletions':'Del', 'ref':'Ref'}
                    legend_elements = []
                    for nuct in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions', 'ref']: 
                        nucpatch = mpatches.Patch(label=nucDrev.get(nuct), facecolor=nucDcolors.get(nuct), edgecolor='black', linewidth=0.05)
                        legend_elements.append(nucpatch)
                    axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.325), ncol=3,prop={'size': 4}, frameon=False, columnspacing=0.5, handletextpad=0.3)

            plt.savefig(opath + args.o + '-' + iso + '-mismatchcoverage.pdf', bbox_inches='tight')
            plt.close()

        return

    def plotCoverageInput(self, data, featurestoplot, opath):

        sns.set_theme(style="white",font_scale=1)

        print("Generating coverage plots for tRNAs of interest...")
        nucD = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}
        # nucDcolors = {'adenines':'tab:green','cytosines':'tab:blue','thymines':'tab:red','guanines':'tab:orange','deletions':'grey', 'ref':'#D3D3D3'}
        # nucDcolors = {'adenines':'#00BA38','cytosines':'#619CFF','thymines':'#F8766D','guanines':'#CD9600','deletions':'grey', 'ref':'#D3D3D3'}
        nucDcolors = {'adenines':'#7CAE00','cytosines':'#619CFF','thymines':'#F8766D','guanines':'#CD9600','deletions':'grey', 'ref':'#D3D3D3'}
        # fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))

        # fig, axes = plt.subplots(max(2, len(self.conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*4,max(2,len(self.conditions))*1.2))
        # fig.subplots_adjust(hspace=0.6, wspace=0.2)

        fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
        fig.subplots_adjust(hspace=0.6, wspace=0.2)

        for idx, tRNA in enumerate(featurestoplot):
            print('Generating coverage plots for...' + tRNA)
            trna_data = data[data['Feature'] == tRNA]
            trna_data = trna_data[trna_data['actualbase'] != '-']
            trna_data = trna_data[['Feature', 'Sample', 'position', 'actualbase', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']]
            sprinzl = list(trna_data['position'].drop_duplicates())
            sprinzl_x = {sprin:featurepos for featurepos,sprin in enumerate(sprinzl)}
            trna_data['xlabel'] = trna_data['position'] + ':' + trna_data['actualbase']
            xlabels = list(trna_data['xlabel'].drop_duplicates())
            xticks = [s+0.1 for s in range(len(xlabels))] 
            sprinzl_xlabel = {featurepos:sprin for featurepos,sprin in enumerate(xlabels)}
            xorder = [featurepos for featurepos,sprin in enumerate(sprinzl)]
            trna_data['Sample'] = trna_data['Sample'].map(self.conditionDict)
            trna_data_mean = trna_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).sum().reset_index()

        #     findmax = trna_data_mean['adenines'] + trna_data_mean['thymines'] + trna_data_mean['cytosines'] + trna_data_mean['guanines'] + trna_data_mean['deletions']
        #     ymaxval = round(1.1 * findmax.max())
        
            trna_data_mean['indexbase'] = trna_data_mean['position'].map(sprinzl_x)

            for qdx, cond in enumerate(self.conditions):
                trna_data_mean_condition = trna_data_mean[trna_data_mean['Sample'] == cond]

                findmax = trna_data_mean_condition['adenines'] + trna_data_mean_condition['thymines'] + trna_data_mean_condition['cytosines'] + trna_data_mean_condition['guanines'] + trna_data_mean_condition['deletions']
                ymaxval = max(20,round(1.1 * findmax.max()))

                for pos in xorder:

                    trna_data_mean_pos = trna_data_mean_condition[trna_data_mean_condition['indexbase'] == pos]
                    realbase = trna_data_mean_pos.iloc[0]['actualbase']
                    realbasecounts = float(trna_data_mean_pos[nucD.get(realbase)])
                    trna_data_mean_pos = trna_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                    conditionslice = trna_data_mean_pos[trna_data_mean_pos.columns[4:-1]].T
                    conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)

                    if float(conditionslice.max()) >= 1:
                        conditionslice = conditionslice.T
                        bottom = 0
                        for base in list(conditionslice.columns):
                            axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], linewidth=1,color=nucDcolors.get(base), edgecolor=nucDcolors.get(base), bottom=bottom)
                            bottom += conditionslice.loc['y'][base]
                        axes[idx][qdx].bar(pos, realbasecounts, color='#D3D3D3', linewidth=1, edgecolor='#D3D3D3', bottom=bottom)
                    else:
                        axes[idx][qdx].bar(pos, realbasecounts, color='#D3D3D3', linewidth=1, edgecolor='#D3D3D3')


                axes[idx][qdx].set_title(tRNA.replace('tRNA-', '')  + ' | ' + cond, fontsize=8, loc='left', pad=2)
                axes[idx][qdx].set_xlim(-1, len(xlabels))
                axes[idx][qdx].set_xticks(xticks)
                axes[idx][qdx].set_xticklabels(xlabels, fontsize=3, rotation=90, horizontalalignment='center')

                axes[idx][qdx].set_ylim(0, ymaxval)
                axes[idx][qdx].set_yticks([0,ymaxval])
                axes[idx][qdx].set_yticklabels(['0',str(ymaxval)], fontsize=5)
                axes[idx][qdx].tick_params(axis='x', which='major', pad=0.5, width=0.25,length=0.25)

                for axis in ['top', 'bottom', 'left', 'right']:
                        axes[idx][qdx].spines[axis].set_linewidth(0.2)

                trnafeatlist = [sprinzl_x.get('9') + 0.3,sprinzl_x.get('21') + 0.3,sprinzl_x.get('26') + 0.3,sprinzl_x.get('38') + 0.3,sprinzl_x.get('48') + 0.3,sprinzl_x.get('60') + 0.3]
                trnawidthlist = [(sprinzl_x.get('13') + 0.5) - (sprinzl_x.get('9') + 0.4),(sprinzl_x.get('25') + 0.5) - (sprinzl_x.get('21') + 0.5),(sprinzl_x.get('31') + 0.5) - (sprinzl_x.get('26') + 0.5),(sprinzl_x.get('43') + 0.5) - (sprinzl_x.get('38') + 0.5),(sprinzl_x.get('53') + 0.5) - (sprinzl_x.get('48') + 0.5),(sprinzl_x.get('65') + 0.5) - (sprinzl_x.get('60') + 0.5)]

                featcolorD = {'dloop':'#FF9999','ac':'#99FF99','tloop':'#99FFFF'}
                for z,ss  in enumerate(['dloop', 'dloop', 'ac', 'ac', 'tloop', 'tloop']):
                    axes[idx][qdx].axvspan(trnafeatlist[z], trnafeatlist[z] + trnawidthlist[z] + 0.3, facecolor=featcolorD.get(ss), alpha=0.2, zorder=-100)

                nucDrev = {'adenines':'A','cytosines':'C','thymines':'T','guanines':'G','deletions':'Del', 'ref':'Ref'}
                legend_elements = []
                for nuct in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions', 'ref']: 
                    nucpatch = mpatches.Patch(label=nucDrev.get(nuct), facecolor=nucDcolors.get(nuct), edgecolor='black', linewidth=0.05)
                    legend_elements.append(nucpatch)
                axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.3), ncol=3,prop={'size': 4}, frameon=False, columnspacing=0.5, handletextpad=0.3)

        plt.savefig(opath + args.o + '-mismatchcoverage-query.pdf', bbox_inches='tight')
        plt.savefig(opath + args.o + '-mismatchcoverage-query.svg', bbox_inches='tight')
        plt.close()
        return

    def plotCoverageInputMapping(self, data, featurestoplot, opath):

        sns.set_theme(style="white",font_scale=1)

        print("Generating coverage plots for tRNAs of interest...")

        # mapDcolors = {'uniquecoverage':'#D3D3D3','multitrnacoverage':'#D3D3D3','multianticodoncoverage':'#D3D3D3','multiaminocoverage':'#D3D3D3', 'coverage':'#D3D3D3'}
        mapDcolors = {'uniquecoverage':'#A781BA','multitrnacoverage':'#00BEC4','multianticodoncoverage':'#7cae00','multiaminocoverage':'#f8766d', 'coverage':'#00BFC4'}
        fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
        fig.subplots_adjust(hspace=0.6, wspace=0.2)

        # fig, axes = plt.subplots(max(2, len(self.conditions)),max(2,len(featurestoplot)),sharey='col', sharex=False, figsize=(max(2, len(featurestoplot))*4,max(2,len(self.conditions))*1.2))
        # fig.subplots_adjust(hspace=0.6, wspace=0.2)
        
        for idx, tRNA in enumerate(featurestoplot):
            print('Generating coverage plots for...' + tRNA)
            trna_data = data[data['Feature'] == tRNA]
            trna_data = trna_data[trna_data['actualbase'] != '-']
            trna_data = trna_data[['Feature', 'Sample', 'position', 'actualbase', 'coverage', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage']]
            sprinzl = list(trna_data['position'].drop_duplicates())
            sprinzl_x = {sprin:featurepos for featurepos,sprin in enumerate(sprinzl)}
            trna_data['xlabel'] = trna_data['position'] + ':' + trna_data['actualbase']
            # xlabels = list(trna_data['xlabel'].drop_duplicates())
            xlabels = list(trna_data['xlabel'].drop_duplicates())
            xticks = [s+0.1 for s in range(len(xlabels))] 
            sprinzl_xlabel = {featurepos:sprin for featurepos,sprin in enumerate(xlabels)}
            xorder = [featurepos for featurepos,sprin in enumerate(sprinzl)]
            trna_data['Sample'] = trna_data['Sample'].map(self.conditionDict)
            trna_data_mean = trna_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).mean().reset_index()

            findmax = trna_data_mean['coverage']
            ymaxval = round(1.1 * findmax.max())

            findxmax = trna_data_mean['position']
            ymaxval = round(1.1 * findmax.max())
        
            trna_data_mean['indexbase'] = trna_data_mean['position'].map(sprinzl_x)
            trna_data_mean = trna_data_mean.sort_values(by=['indexbase'])

            findxmax = trna_data_mean['indexbase']
            xmaxval = findxmax.max()

            for qdx, cond in enumerate(self.conditions):

                trna_data_mean_condition = trna_data_mean[trna_data_mean['Sample'] == cond]

                for mappability in ['uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage']:
                    if mappability == 'uniquecoverage':
                        axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=1,color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability))
                        bottom = trna_data_mean_condition[mappability]
                    else:
                        axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=1,color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability), bottom=bottom)
                        bottom += trna_data_mean_condition[mappability]

                axes[idx][qdx].set_title(tRNA.replace('tRNA-', '')  + ' | ' + cond, fontsize=8, loc='left', pad=2)
                axes[idx][qdx].set_xlim(-1, len(xlabels))
                axes[idx][qdx].set_xticks(xticks)
                axes[idx][qdx].set_xticklabels(xlabels, fontsize=3, rotation=90, horizontalalignment='center')

                axes[idx][qdx].set_ylim(0, ymaxval)
                axes[idx][qdx].set_yticks([0,ymaxval])
                axes[idx][qdx].set_yticklabels(['0',str(ymaxval)], fontsize=5)
                axes[idx][qdx].tick_params(axis='x', which='major', pad=0.5, width=0.25,length=0.25)

                for axis in ['top', 'bottom', 'left', 'right']:
                        axes[idx][qdx].spines[axis].set_linewidth(0.2)

                # trnafeatlist = [sprinzl_x.get('9') + 0.3,sprinzl_x.get('21') + 0.3,sprinzl_x.get('26') + 0.3,sprinzl_x.get('38') + 0.3,sprinzl_x.get('48') + 0.3,sprinzl_x.get('60') + 0.3]
                # trnawidthlist = [(sprinzl_x.get('13') + 0.5) - (sprinzl_x.get('9') + 0.4),(sprinzl_x.get('25') + 0.5) - (sprinzl_x.get('21') + 0.5),(sprinzl_x.get('31') + 0.5) - (sprinzl_x.get('26') + 0.5),(sprinzl_x.get('43') + 0.5) - (sprinzl_x.get('38') + 0.5),(sprinzl_x.get('53') + 0.5) - (sprinzl_x.get('48') + 0.5),(sprinzl_x.get('65') + 0.5) - (sprinzl_x.get('60') + 0.5)]

                # featcolorD = {'dloop':'#FF9999','ac':'#99FF99','tloop':'#99FFFF'}
                # for z,ss  in enumerate(['dloop', 'dloop', 'ac', 'ac', 'tloop', 'tloop']):
                #     axes[idx][qdx].axvspan(trnafeatlist[z], trnafeatlist[z] + trnawidthlist[z] + 0.3, facecolor=featcolorD.get(ss), alpha=0.2, zorder=-100)

                # nucDrev = {'adenines':'A','cytosines':'C','thymines':'T','guanines':'G','deletions':'Del', 'ref':'Ref'}
                legend_elements = []
                for nuct in ['uniquecoverage', 'multitrnacoverage']:#, 'multianticodoncoverage', 'multiaminocoverage']: 
                    nucpatch = mpatches.Patch(label=nuct.replace('coverage', ''), facecolor=mapDcolors.get(nuct), edgecolor='black', linewidth=0.05)
                    legend_elements.append(nucpatch)
                axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.325), ncol=1,prop={'size': 4}, frameon=False, columnspacing=0.5, handletextpad=0.3)

        plt.savefig(opath + args.o + '-mappabilitycoverage-query.pdf', bbox_inches='tight')
        # plt.savefig(opath + args.o + '-mappabilitycoverage-query.svg', bbox_inches='tight')
        plt.close()
        return

    def plotCoverageMapping(self, data, isodecoders, isotypes, opath):

        sns.set_theme(style="white",font_scale=1)
        mapDcolors = {'uniquecoverage':'#A781BA','multitrnacoverage':'#00BEC4','multianticodoncoverage':'#7cae00','multiaminocoverage':'#f8766d', 'coverage':'#00BFC4'}
        # mapDcolors = {'uniquecoverage':'#D3D3D3','multitrnacoverage':'#D3D3D3','multianticodoncoverage':'#7CAE00','multiaminocoverage':'#F8766D', 'coverage':'#D3D3D3'}

        for iso in isotypes:
            data_iso = data[data['Isotype'] == iso]

            print("Generating coverage plots for..." + iso)
            featurestoplot = [i for i in list(data_iso['Feature'].drop_duplicates()) if i in isodecoders]
            featurestoplot = sorted(featurestoplot, key=lambda x: (x.split('-')[2], int(x.split('-')[3])))

            fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
            fig.subplots_adjust(hspace=0.6, wspace=0.2)

            for idx, tRNA in enumerate(featurestoplot):
                trna_data = data[data['Feature'] == tRNA]
                trna_data = trna_data[trna_data['actualbase'] != '-']
                trna_data = trna_data[['Feature', 'Sample', 'position', 'actualbase', 'coverage', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage']]
                sprinzl = list(trna_data['position'].drop_duplicates())
                sprinzl_x = {sprin:featurepos for featurepos,sprin in enumerate(sprinzl)}
                trna_data['xlabel'] = trna_data['position'] + ':' + trna_data['actualbase']
                # xlabels = list(trna_data['xlabel'].drop_duplicates())
                xlabels = list(trna_data['xlabel'].drop_duplicates())
                xticks = [s+0.1 for s in range(len(xlabels))] 
                sprinzl_xlabel = {featurepos:sprin for featurepos,sprin in enumerate(xlabels)}
                xorder = [featurepos for featurepos,sprin in enumerate(sprinzl)]
                trna_data['Sample'] = trna_data['Sample'].map(self.conditionDict)
                trna_data_mean = trna_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).mean().reset_index()

                findmax = trna_data_mean['coverage']
                ymaxval = round(1.1 * findmax.max())

                findxmax = trna_data_mean['position']
                ymaxval = round(1.1 * findmax.max())
            
                trna_data_mean['indexbase'] = trna_data_mean['position'].map(sprinzl_x)
                trna_data_mean = trna_data_mean.sort_values(by=['indexbase'])

                findxmax = trna_data_mean['indexbase']
                xmaxval = findxmax.max()

                for qdx, cond in enumerate(self.conditions):

                    trna_data_mean_condition = trna_data_mean[trna_data_mean['Sample'] == cond]

                    for mappability in ['uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage']:
                        if mappability == 'uniquecoverage':
                            axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=1,color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability))
                            bottom = trna_data_mean_condition[mappability]
                        else:
                            axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=1,color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability), bottom=bottom)
                            bottom += trna_data_mean_condition[mappability]


                    axes[idx][qdx].set_title(tRNA.replace('tRNA-', '')  + ' | ' + cond, fontsize=8, loc='left', pad=2)
                    axes[idx][qdx].set_xlim(-1, len(xlabels))
                    axes[idx][qdx].set_xticks(xticks)
                    axes[idx][qdx].set_xticklabels(xlabels, fontsize=2.5, rotation=90, horizontalalignment='center')

                    axes[idx][qdx].set_ylim(0, ymaxval)
                    axes[idx][qdx].set_yticks([0,ymaxval])
                    axes[idx][qdx].set_yticklabels(['0',str(ymaxval)], fontsize=5)
                    axes[idx][qdx].tick_params(axis='x', which='major', pad=0.5, width=0.25,length=0.25)

                    for axis in ['top', 'bottom', 'left', 'right']:
                            axes[idx][qdx].spines[axis].set_linewidth(0.2)

                    # trnafeatlist = [sprinzl_x.get('9') + 0.3,sprinzl_x.get('21') + 0.3,sprinzl_x.get('26') + 0.3,sprinzl_x.get('38') + 0.3,sprinzl_x.get('48') + 0.3,sprinzl_x.get('60') + 0.3]
                    # trnawidthlist = [(sprinzl_x.get('13') + 0.5) - (sprinzl_x.get('9') + 0.4),(sprinzl_x.get('25') + 0.5) - (sprinzl_x.get('21') + 0.5),(sprinzl_x.get('31') + 0.5) - (sprinzl_x.get('26') + 0.5),(sprinzl_x.get('43') + 0.5) - (sprinzl_x.get('38') + 0.5),(sprinzl_x.get('53') + 0.5) - (sprinzl_x.get('48') + 0.5),(sprinzl_x.get('65') + 0.5) - (sprinzl_x.get('60') + 0.5)]

                    # featcolorD = {'dloop':'#FF9999','ac':'#99FF99','tloop':'#99FFFF'}
                    # for z,ss  in enumerate(['dloop', 'dloop', 'ac', 'ac', 'tloop', 'tloop']):
                    #     axes[idx][qdx].axvspan(trnafeatlist[z], trnafeatlist[z] + trnawidthlist[z] + 0.3, facecolor=featcolorD.get(ss), alpha=0.2, zorder=-100)

                    # nucDrev = {'adenines':'A','cytosines':'C','thymines':'T','guanines':'G','deletions':'Del', 'ref':'Ref'}
                    legend_elements = []
                    for nuct in ['uniquecoverage', 'multitrnacoverage']:#, 'multianticodoncoverage', 'multiaminocoverage']: 
                        nucpatch = mpatches.Patch(label=nuct.replace('coverage', ''), facecolor=mapDcolors.get(nuct), edgecolor='black', linewidth=0.05)
                        legend_elements.append(nucpatch)
                    axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.325), ncol=1,prop={'size': 4}, frameon=False, columnspacing=0.5, handletextpad=0.3)

            # plt.tight_layout()
            plt.savefig(opath + args.o + '-' + iso + '-mappabilitycoverage.pdf', bbox_inches='tight')
            # plt.savefig(opath + args.o + '-' + iso + '-mappabilitycoverage.svg', bbox_inches='tight')
            plt.close()
        return

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--cov', required=True, help='-coverage.txt file from tRAX output')
    ap.add_argument('--samples', required=True, help='samples.txt file from tRAX to group samples')
    ap.add_argument('--mode', required=True, help='mismatch or mappability; if you want to color the base coverage based on nucleotide misincorporations or read mappability',default='mismatch')
    ap.add_argument('--trnas', required=False, help='tRNA isodecoders of interest', nargs='+', default=[])
    ap.add_argument('--org', required=False, help='euk, bact, arch, or mito', default='euk')
    ap.add_argument('--o', required=True, help='path for output files')
    args = ap.parse_args()
    
    ## Generates tMAP object
    tMAPcoverageplots = tMAPcoverageplots(pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq']))

    startTime = datetime.now()

    ## Reads in coverage file and generates a list of isoacceptors and isodecoders
    coveragedata, isodecoderlist, isoacceptorlist, isotypelist = tMAPcoverageplots.readCoverage(args.cov)

    ## Checks for output directory and creates if not present
    opath = os.getcwd() + '/coverageplots/'
    isdir = os.path.isdir(opath)
    if isdir == True:
        pass
    else:
        os.mkdir(opath)

    if len(args.trnas) == 0 and args.mode == 'mismatch':
        ## Generates coverage plots for each isotype
        tMAPcoverageplots.plotCoverage(coveragedata, isodecoderlist, isotypelist, opath)
    elif len(args.trnas) > 0 and args.mode == 'mismatch':
        ## Generates coverage plots for each input tRNA
        tMAPcoverageplots.plotCoverageInput(coveragedata, args.trnas, opath)
    elif len(args.trnas) == 0 and args.mode == 'coverage':
        ## Generates coverage plots for each isotype
        tMAPcoverageplots.plotCoverageMapping(coveragedata, isodecoderlist, isotypelist, opath)
    elif len(args.trnas) > 0 and args.mode == 'coverage':
        ## Generates coverage plots for each input tRNA
        tMAPcoverageplots.plotCoverageInputMapping(coveragedata, args.trnas, opath)

    print('\nTime elasped: ', datetime.now() - startTime)