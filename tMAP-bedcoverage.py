#!/usr/bin/env python3
#Aidan Manning 2.21.23
#Script to determine sites of feat modficiation from sequencing data and generates some plots
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
import pybedtools as bedtools
import pysam
from re import search
from multiprocessing import Pool as ProcessPool
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class tMAPbedcoverage(object):

    def __init__(self, opath, sampleinfo):

        self.reverseComplement = {'A':'T','C':'G','G':'C','T':'A','-':'-', 'N':'N'}
        lociList = []
        for line in open(opath + args.o + '-regionsofinterest.bed'):
            line = line.strip().split('\t')

            dsplit = line[0].split('::')
            ddsplit = dsplit[1].split('(')
            dddsplit = ddsplit[1].split(')')
            possplit = ddsplit[0].split(':')
            possplit2 = possplit[1].split('-')
            featurename = dsplit[0].strip().replace(' ', '')
            chromosome = dsplit[1]
            strand = dddsplit[0]
            start = possplit2[0]
            end = possplit2[1]

            featureseq = line[1].strip()
            if search('feat', featurename):
                featureseq = featureseq + 'CCA'
            else:
                featureseq = featureseq

            for j,nuc in enumerate(featureseq):
                lociList.append((featurename.strip(), j + 1, nuc))

        self.locidf = pd.DataFrame.from_records(lociList, columns=['Feature', 'position', 'actualbase'])
        self.conditionDict = dict(zip(list(sampleinfo['sample']), list(sampleinfo['condition'])))    
        self.conditions = list(sampleinfo['condition'].drop_duplicates())

        return
    
    def reformatSequences(self, loci):

        lociList = []
        for line in open(loci):
            line = line.strip().split('\t')

            dsplit = line[0].split('::')
            ddsplit = dsplit[1].split('(')
            dddsplit = ddsplit[1].split(')')
            possplit = ddsplit[0].split(':')
            possplit2 = possplit[1].split('-')
            featurename = dsplit[0].strip().replace(' ', '')
            chromosome = dsplit[1]
            strand = dddsplit[0]
            start = possplit2[0]
            end = possplit2[1]

            featureseq = line[1].strip()
            if search('feat', featurename):
                featureseq = featureseq + 'CCA'
            else:
                featureseq = featureseq

            for j,nuc in enumerate(featureseq):
                lociList.append((featurename.strip(), j + 1, nuc))

        regiondf = pd.DataFrame.from_records(lociList, columns=['Feature', 'position', 'actualbase'])

        return regiondf
    
    def getReadCovClip(self, alignments):

        nucD = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines', '-':'deletions'}
        samplename = Path(alignments).stem

        samfile = pysam.AlignmentFile(alignments, "rb")
        print('Performing read pileup for... ' + samplename)

        dataList = []

        for feats in open(args.regions):

            feats = feats.strip().split('\t')

            strt = int(feats[1])
            ed = int(feats[2])
            if search('feat', feats[3]) and feats[5] == '+':
                ed = int(feats[2]) + 3
            elif search('miR', feats[3]) and feats[5] == '+':
                ed = int(feats[2]) + 1
            elif search('feat', feats[3]) and feats[5] == '-':
                strt = int(feats[1]) - 3
            elif search('miR', feats[3]) and feats[5] == '-':
                strt = int(feats[1]) - 1

            featlen = (ed - strt) + 1
            featdata = {q: {z: 0 for z in ['coverage','adenines', 'cytosines', 'thymines', 'guanines', 'deletions', 'readstarts', 'readends']} for q in range(1, featlen)}

            if feats[5] == '+':
                for read in samfile.fetch(feats[0], strt, ed):
                    # if not read.is_secondary or read.is_supplementary or read.is_unmapped and read.mapping_qality > 30:
                    if not read.is_supplementary or read.is_unmapped and read.mapping_qality > 30:
                        
                        ref_start = read.reference_start
                        ref_end = read.reference_end
                        read_cigart = read.cigartuples
                        
                        clip_start = ref_start - strt
                        clip_end = ed - ref_end

                        if clip_start + args.clipping >= 0 and clip_end + args.clipping >= 0 and len(read_cigart) == 1:

                            queryseq = read.query_sequence

                            if clip_start < 0 and clip_end >= 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):]
                                featend = featstart + readlength
                            elif clip_start < 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):clip_end]
                                featend = featstart + readlength
                            elif clip_start >= 0 and clip_end < 0:
                                queryseq = queryseq[:clip_end]

                                featstart = (ref_start - strt) + 1
                            else:
                                featstart = (ref_start - strt) + 1


                            readlength = len(queryseq)

                            featend = featstart + readlength
                            currcount = 0

                            for idx in range(featstart, featend):
                                querybase = queryseq[currcount].upper()
                                featdata[idx]['coverage'] += 1
                                if idx == featstart and clip_start >= 0:
                                    featdata[idx]['readstarts'] += 1
                                elif idx == featend - 1 and clip_end >= 0:
                                    featdata[idx]['readends'] += 1
                                else:
                                    pass

                                if querybase == 'N':
                                        currcount += 1
                                else:
                                    try:
                                        featdata[idx][nucD.get(querybase)] += 1
                                        currcount += 1
                                    except:
                                        print(queryseq[currcount])

                        elif clip_start + args.clipping >= 0 and clip_end + args.clipping >= 0 and len(read_cigart) > 1:
                        # elif ref_start >= strt and ref_end <= ed and len(read_cigart) > 1:

                            queryseq = read.query_sequence
                            quals = read.query_qualities
                            fullstring = []
                            querypos = 0

                            for tup in read_cigart:
                                operation = int(tup[0])
                                oplen = int(tup[1])
                                if operation == 0:
                                    fullstring.append(queryseq[querypos:(oplen + querypos)])
                                    querypos += oplen
                                elif operation == 1:
                                    quals.pop(querypos)
                                    querypos += oplen
                                elif operation == 2:
                                    fullstring.append('-'*oplen)

                            queryseq = ''.join(fullstring)

                            if clip_start < 0 and clip_end >= 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):]
                                featend = featstart + readlength
                            elif clip_start < 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):clip_end]
                                featend = featstart + readlength
                            elif clip_start >= 0 and clip_end < 0:
                                queryseq = queryseq[:clip_end]

                                featstart = (ref_start - strt) + 1
                            else:
                                featstart = (ref_start - strt) + 1


                            readlength = len(queryseq)

                            featend = featstart + readlength
                            currcount = 0

                            for idx in range(featstart, featend):
                                querybase = queryseq[currcount].upper()
                                featdata[idx]['coverage'] += 1
                                if idx == featstart and clip_start >= 0:
                                    featdata[idx]['readstarts'] += 1
                                elif idx == featend - 1 and clip_end >= 0:
                                    featdata[idx]['readends'] += 1
                                else:
                                    pass

                                if querybase == 'N':
                                        currcount += 1
                                else:
                                    try:
                                        featdata[idx][nucD.get(querybase)] += 1
                                        currcount += 1
                                    except:
                                        print(queryseq[currcount])
                        else:
                            pass
                    else:
                        pass

            elif feats[5] == '-':
                for read in samfile.fetch(feats[0], strt, ed):
                    # if not read.is_secondary or read.is_supplementary or read.is_unmapped and read.mapping_qality > 30:
                    if not read.is_supplementary or read.is_unmapped and read.mapping_qality > 30:
                        ref_start = read.reference_start
                        ref_end = read.reference_end
                        read_cigart = read.cigartuples

                        clip_start = ref_start - strt
                        clip_end = ed - ref_end

                        if clip_start + args.clipping >= 0 and clip_end + args.clipping >= 0 and len(read_cigart) == 1:


                            queryseq = read.query_sequence
                            # print(queryseq, clip_start, clip_end, 'old')

                            if clip_start < 0 and clip_end >= 0:
                                featstart = ed - ref_end + 1
                                queryseq = queryseq[abs(clip_start):]

                            elif clip_start < 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):clip_end]

                            elif clip_start >= 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[:clip_end]
                            else:
                                featstart = ed - ref_end + 1


                            readlength = len(queryseq)

                            featend = featstart + readlength

                            queryseq = queryseq[::-1]

                            currcount = 0

                            for idx in range(featstart, featend):
                                # idx = featend - idx

                                querybase = self.reverseComplement.get(queryseq[currcount].upper())
                                featdata[idx]['coverage'] += 1
                                if idx == featstart and clip_end >= 0:
                                    featdata[idx]['readstarts'] += 1
                                elif idx == featend - 1 and clip_start >= 0:
                                    featdata[idx]['readends'] += 1
                                else:
                                    pass

                                if querybase == 'N':
                                    currcount += 1
                                else:
                                    try:
                                        featdata[idx][nucD.get(querybase)] += 1
                                        currcount += 1
                                    except:
                                        print(queryseq[currcount])

                        elif clip_start + args.clipping >= 0 and clip_end + args.clipping >= 0 and len(read_cigart) > 1:


                            queryseq = read.query_sequence
                            quals = read.query_qualities
                            fullstring = []
                            querypos = 0

                            for tup in read_cigart:
                                operation = int(tup[0])
                                oplen = int(tup[1])
                                if operation == 0:
                                    fullstring.append(queryseq[querypos:(oplen + querypos)])
                                    querypos += oplen
                                elif operation == 1:
                                    quals.pop(querypos)
                                    querypos += oplen
                                elif operation == 2:
                                    fullstring.append('-'*oplen)

                            queryseq= ''.join(fullstring)

                            if clip_start < 0 and clip_end >= 0:
                                featstart = ed - ref_end + 1
                                queryseq = queryseq[abs(clip_start):]

                            elif clip_start < 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[abs(clip_start):clip_end]

                            elif clip_start >= 0 and clip_end < 0:
                                featstart = 1
                                queryseq = queryseq[:clip_end]
                            else:
                                featstart = ed - ref_end + 1


                            readlength = len(queryseq)

                            featend = featstart + readlength

                            queryseq = queryseq[::-1]

                            currcount = 0

                            for idx in range(featstart, featend):
                                # idx = featend - idx
                                # print(queryseq[currcount])
                                querybase = self.reverseComplement.get(queryseq[currcount].upper())

                                featdata[idx]['coverage'] += 1
                                if idx == featstart and clip_end >= 0:
                                    featdata[idx]['readstarts'] += 1
                                elif idx == featend - 1 and clip_start >= 0:
                                    featdata[idx]['readends'] += 1
                                else:
                                    pass

                                if querybase == 'N':
                                    currcount += 1
                                else:
                                    try:
                                        featdata[idx][nucD.get(querybase)] += 1
                                        currcount += 1
                                    except:
                                        print(queryseq[currcount])
                        else:
                            pass
                    else:
                        pass

            featdf = pd.DataFrame.from_dict(featdata, orient='index')
            featdf['Feature'] = feats[3].strip()
            featdf['Sample'] = samplename

            featdf = featdf.reset_index().rename(columns={'index':'position'})
            dataList.append(featdf)

        df = pd.concat(dataList)

        merger = self.locidf.merge(df, on=['Feature', 'position'], how='left').rename(columns={'-':'deletions'})

        nucList = ['A', 'C', 'T', 'G']
        dfL = []
        # print('Performing Mismatch % Calculation for... ' + samplename)
        for base in nucList:
            nuc_df = merger[merger['actualbase'] == base]
            # nuc_df['coverage'] = nuc_df['C'] + nuc_df['G'] + nuc_df['T'] + nuc_df['A'] + nuc_df['D']
            if base == 'A':
                nuc_df['mismatchedbases'] = nuc_df['cytosines'] + nuc_df['guanines'] + nuc_df['thymines'] + nuc_df['deletions']
            elif base == 'C':
                nuc_df['mismatchedbases'] = nuc_df['adenines'] + nuc_df['guanines'] + nuc_df['thymines'] + nuc_df['deletions']
            elif base == 'G':
                nuc_df['mismatchedbases'] = nuc_df['adenines'] + nuc_df['cytosines'] + nuc_df['thymines'] + nuc_df['deletions']
            elif base == 'T':
                nuc_df['mismatchedbases'] = nuc_df['adenines'] + nuc_df['guanines'] + nuc_df['cytosines'] + nuc_df['deletions']
            dfL.append(nuc_df)
        mmdf = pd.concat(dfL)

        full_merge = merger.merge(mmdf, on=['Feature', 'Sample', 'position', 'actualbase','coverage', 'adenines', 'cytosines', 'guanines', 'thymines', 'deletions','readstarts', 'readends'], how='left').fillna(0)
        # full_merge = full_merge[full_merge['coverage'] >= 20]
        # full_merge['% Mismatch'] = (full_merge['mismatchedbases'] / full_merge['coverage'])*100
        # full_merge = full_merge[full_merge['% Mismatch'] >= 20]
        return full_merge

    def getFeatLengths(self, d, feats):
        
        featL = []

        for feat in feats:
            d_f = d[d['Feature'] == feat]
            length = d_f['position'].max()
            featL.append(length*0.8)

        return featL

    def plotCoverageInput(self, data, opath):

        sns.set_theme(style="white",font_scale=1)

        print("Generating coverage plots for regions of interest...")
        nucD = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}
        # nucDcolors = {'adenines':'tab:green','cytosines':'tab:blue','thymines':'tab:red','guanines':'tab:orange','deletions':'grey', 'ref':'#D3D3D3'}
        nucDcolors = {'adenines':'#7CAE00','cytosines':'#619CFF','thymines':'#F8766D','guanines':'#CD9600','deletions':'grey', 'ref':'#D3D3D3'}

        featurestoplot = list(data['Feature'].drop_duplicates())

        getfeatlengths = tMAPbedcov.getFeatLengths(data, featurestoplot)


        # fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
        if max(getfeatlengths) > 100:
            fig, axes = plt.subplots(max(2, len(self.conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*2,max(2,len(self.conditions))*0.8), gridspec_kw={'width_ratios': getfeatlengths})
        else:
            fig, axes = plt.subplots(max(2, len(self.conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*2,max(2,len(self.conditions))*0.8), gridspec_kw={'width_ratios': getfeatlengths})

        fig.subplots_adjust(hspace=0.8, wspace=0.3)

        for qdx, feat in enumerate(featurestoplot):

            feat_data = data[data['Feature'] == feat]
            feat_data = feat_data[feat_data['actualbase'] != '-']
            feat_data = feat_data[['Feature', 'Sample', 'position', 'actualbase', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']]

            sprinzl = list(feat_data['position'].drop_duplicates())
            sprinzl_x = {sprin:featurepos for featurepos,sprin in enumerate(sprinzl)}
            feat_data['xlabel'] = feat_data['position'].astype(str) + ':' + feat_data['actualbase'].replace('T','U')
            xlabels = list(feat_data['xlabel'].drop_duplicates())
            xticks = [s+0.1 for s in range(len(xlabels))] 
            sprinzl_xlabel = {featurepos:sprin for featurepos,sprin in enumerate(xlabels)}
            xorder = [featurepos for featurepos,sprin in enumerate(sprinzl)]
            feat_data['Sample'] = feat_data['Sample'].map(self.conditionDict)
            feat_data_mean = feat_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).sum().reset_index()

        #     findmax = feat_data_mean['adenines'] + feat_data_mean['thymines'] + feat_data_mean['cytosines'] + feat_data_mean['guanines'] + feat_data_mean['deletions']
        #     ymaxval = round(1.1 * findmax.max())
        
            feat_data_mean['indexbase'] = feat_data_mean['position'].map(sprinzl_x)

            for idx, cond in enumerate(self.conditions):
                feat_data_mean_condition = feat_data_mean[feat_data_mean['Sample'] == cond]

                findmax = feat_data_mean_condition['adenines'] + feat_data_mean_condition['thymines'] + feat_data_mean_condition['cytosines'] + feat_data_mean_condition['guanines'] + feat_data_mean_condition['deletions']
                ymaxval = max(20,round(1.1 * findmax.max()))

                for pos in xorder:

                    feat_data_mean_pos = feat_data_mean_condition[feat_data_mean_condition['indexbase'] == pos]
                    realbase = feat_data_mean_pos.iloc[0]['actualbase']
                    realbasecounts = float(feat_data_mean_pos[nucD.get(realbase)])
                    feat_data_mean_pos = feat_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                    conditionslice = feat_data_mean_pos[feat_data_mean_pos.columns[4:-1]].T
                    conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)

                    if float(conditionslice.max()) >= 1:
                        conditionslice = conditionslice.T
                        bottom = 0
                        for base in list(conditionslice.columns):
                            axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], linewidth=0,color=nucDcolors.get(base), bottom=bottom)
                            bottom += conditionslice.loc['y'][base]
                        axes[idx][qdx].bar(pos, realbasecounts, width=1, linewidth=0,color='#D3D3D3', bottom=bottom)
                    else:
                        axes[idx][qdx].bar(pos, realbasecounts, width=1, linewidth=0,color='#D3D3D3')


                axes[idx][qdx].set_title(feat.replace('mmu-', '')  + ' | ' + cond, fontsize=8, loc='left', pad=2)
                axes[idx][qdx].set_xlim(-1, len(xlabels))
                axes[idx][qdx].set_xticks(xticks)
                axes[idx][qdx].set_xticklabels(xlabels, fontsize=3, rotation=90, horizontalalignment='center')

                axes[idx][qdx].set_ylim(0, ymaxval)
                axes[idx][qdx].set_yticks([0,ymaxval])
                axes[idx][qdx].set_yticklabels(['0',str(ymaxval)], fontsize=5)
                axes[idx][qdx].tick_params(axis='x', which='major', pad=0.5, width=0.25,length=0.25)

                for axis in ['top', 'bottom', 'left', 'right']:
                        axes[idx][qdx].spines[axis].set_linewidth(0.2)


                nucDrev = {'adenines':'A','cytosines':'C','thymines':'T','guanines':'G','deletions':'Del', 'ref':'Ref'}
                # legend_elements = []
                # for nuct in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions', 'ref']: 
                #     nucpatch = mpatches.Patch(label=nucDrev.get(nuct), color=nucDcolors.get(nuct))
                #     legend_elements.append(nucpatch)
                # # axes[idx][qdx].legend(handles = legend_elements,  loc='center right', bbox_to_anchor=(1.15, 0.5), ncol=1,prop={'size': 4}, frameon=False, columnspacing=0.8)
                # # if getfeatlengths[qdx] < 10:
                # #     axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.175), ncol=6,prop={'size': 4}, frameon=False, columnspacing=0.8)
                # #     axes[idx][qdx].set_title(feat  + ' | ' + cond, fontsize=10, loc='left', pad=6)
                # # else:
                # axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.325), ncol=3,prop={'size': 4}, frameon=False, columnspacing=0.8)


                legend_elements = []
                for nuct in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions', 'ref']: 
                    nucpatch = mpatches.Patch(label=nucDrev.get(nuct).replace('T', 'U'), facecolor=nucDcolors.get(nuct), edgecolor='black', linewidth=0.05)
                    legend_elements.append(nucpatch)
                if getfeatlengths[qdx] < 25:
                    axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.485), ncol=6,prop={'size': 3.8}, frameon=False, columnspacing=0.2, handletextpad=0.3)
                    axes[idx][qdx].set_title(feat.replace('mmu-','')  + ' | ' + cond, fontsize=6, loc='left', pad=2)
                else:
                    axes[idx][qdx].legend(handles = legend_elements,  loc='upper right', bbox_to_anchor=(1, 1.415), ncol=3,prop={'size': 4}, frameon=False, columnspacing=0.5, handletextpad=0.3)

        # plt.tight_layout()
        plt.savefig(opath + args.o + '-mismatchcoverage.pdf', bbox_inches='tight')
        plt.close()
        return

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--fasta', required=True, help='Genome fasta file; needs an index file (.fai), can be generated with samtools faidx fasta.fa')
    ap.add_argument('--regions', required=True, help='bed file containing genomic region of interest')
    ap.add_argument('--threads', required=True, help='number of threads to process the data on; Default=8', default=1)
    ap.add_argument('--alignments', required=True, help='bam file containing read alignments; needs index files (.fai), can be generated with samtools index alignments.bam', nargs='+')
    ap.add_argument('--o', required=True, help='path for output file')
    # ap.add_argument('--annotation', required=False,help='gtf file of genome feature annotations. Format should conatin gene_name "GENE_NAME"; gene_biotype "BIOTYPE";')
    # ap.add_argument('--scale', required=False,help='-sizefactors.txt file from tRAX. Can also generate a custom list')
    ap.add_argument('--clipping', required=False,help='number of bases a read can extend beyond the feature start/end to be included in the analysis; default=0', default=0, type=int)
    ap.add_argument('--samples', required=False, help='samples.txt file from tRAX to group samples')
    ap.add_argument('--plot', required=False, help='whether or not to generate coverage plots for each of the features', action='store_true')
    args = ap.parse_args()

    startTime = datetime.now()

    lociofinterest = bedtools.BedTool(args.regions)
    fasta = bedtools.BedTool(args.fasta)

    opath = os.getcwd() + '/bedcoverage/'
    isdir = os.path.isdir(opath)
    if isdir == True:
        pass
    else:
        os.mkdir(opath)

    print('Generating fasta sequences for the regions of interest...')
    lociofinterest.sequence(fi=fasta, name=True, s=True,tab=True, fo=opath + args.o + '-regionsofinterest.bed')

    ## Generates tMAP object
    tMAPbedcov = tMAPbedcoverage(opath, pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq']))

    with ProcessPool(processes=int(args.threads)) as pool:

        # samplebedcoverage = pool.map(tMAPbedcov.getReadCov2, args.alignments)
        samplebedcoverage = pool.map(tMAPbedcov.getReadCovClip, args.alignments)

    samplebedcoverage_all = pd.concat(samplebedcoverage)
    samplebedcoverage_all = samplebedcoverage_all[['Sample','Feature','position', 'actualbase','coverage','mismatchedbases','readstarts', 'readends','adenines', 'cytosines', 'thymines', 'guanines','deletions']]

    # samplebedcoverage_all = samplebedcoverage_all.sort_values(by=['Sample','Feature','chr','position'])
    samplebedcoverage_all.to_csv(opath + args.o + '-bedcoverage.txt', sep='\t', index=False)

    if args.plot == True:
        tMAPbedcov.plotCoverageInput(samplebedcoverage_all, opath)

    print('\nTime elasped: ', datetime.now() - startTime)