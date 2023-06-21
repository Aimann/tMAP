#!/usr/bin/env python3
#Aidan C. Manning 6.21.23
#Script to determine sites of tRNA modficiation from sequencing data and generates some plots
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import re
from re import search
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools as bedtools
import pysam
import scipy.special as sc
import joblib
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as smm
import warnings
from multiprocessing import Pool as ProcessPool
from pathlib import Path
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
warnings.filterwarnings("ignore")

class tMAPgetmismatch(object):

    def __init__(self):
        return

    def readCoverage(self, df):
        fixSprinzl = {'-1':'-1','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6','7':'7','8':'8','9':'9','10':'10','11':'11','12':'12','13':'13','14':'14','15':'15','16':'16','17':'17','17a':'17a','18':'18','19':'19','20':'20','20a':'20a','20b':'20b','21':'21','22':'22','23':'23','24':'24','25':'25','26':'26','27':'27','28':'28','29':'29','30':'30','31':'31','32':'32','33':'33','34':'34','35':'35','36':'36','37':'37','38':'38','39':'39','40':'40','41':'41','42':'42','43':'43','44':'44','45':'45','e1':'e11','e2':'e12','e3':'e13','e4':'e14','e5':'e15','e6':'e16','e7':'e17','e8':'e1','e9':'e2','e10':'e3','e11':'e4','e12':'e5','e13':'e27','e14':'e26','e15':'e25','e16':'e24','e17':'e23','e18':'e22','e19':'e21','46':'46','47':'47','48':'48','49':'49','50':'50','51':'51','52':'52','53':'53','54':'54','55':'55','56':'56','57':'57','58':'58','59':'59','60':'60','61':'61','62':'62','63':'63','64':'64','65':'65','66':'66','67':'67','68':'68','69':'69','70':'70','71':'71','72':'72','73':'73','74':'74','75':'75','76':'76'}

        ## Reads in the tRAX coverage file and returns a dataframe with the coverage data with some extra annotations for downstream analysis
        data = pd.read_csv(df, sep='\t')
        if args.org == 'euk':
            data = data[data['Feature'].str.startswith('tRNA')]
            data['position'] = data['position'].map(fixSprinzl)
        elif args.org == 'mito':
            data = data[data['Feature'].str.startswith('mt-tRNA')]
            data['position'] = data['position']
        else:
            pass
        data['Isoacceptor'] = data['Feature'].str.replace('(-[0-9]+)', '', regex=True)
        isodecoders = list(data['Feature'].drop_duplicates())
        isoacceptors = list(data['Isoacceptor'].drop_duplicates())

        return data, isodecoders, isoacceptors

    def isodecoderAlign(self, trna, isocov):
        ## Aligns the isodecoders of each isoacceptor group and returns a dataframe with the aligned sequences
        isoacceptor = re.sub('(-[0-9]+)', '', trna)
        isodecoderdifferencelist = []

        sprinzl = list(isocov['position'].drop_duplicates())
        isoacceptseq = isocov[['Feature', 'position', 'actualbase']].drop_duplicates()
        isopiv = isoacceptseq.pivot(index='Feature', columns='position', values='actualbase')
        isopiv = isopiv[sprinzl]

        testvector = isopiv.to_numpy()
        checkfordiffs = list((testvector[0] == testvector).all(0))

        for idx,pos in enumerate(list(isopiv.columns)):
            if checkfordiffs[idx] == True:
                pass
            else:
                isodecoderdifferencelist.append(pos)
                pass
        diffvalues = ','.join(str(i) for i in isodecoderdifferencelist)

        return (isoacceptor, diffvalues) , isopiv

    def compareIsodecoders(self, trna, sprinzlalign):
         ## Compares the isodecoders of each isoacceptor group
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

    def variationFalsePositives(self, cnf, isovar):
        ## Checks for potential false positives in the mismatch data by checking if the mismatch position is also a position of variation in the isoacceptor group
        cnf['varbarcode'] = cnf['Feature'].astype(str) + '|' + cnf['position']
        test_isodecoders = list(cnf['varbarcode'].drop_duplicates())
        vartrnaL = {trna:[] for trna in test_isodecoders}
        varnucL = {trna:[] for trna in test_isodecoders}
        for isod in test_isodecoders:
            sub_cnf = cnf[cnf['varbarcode'] == isod]
            sub_isovarx = isovar[isovar['tRNA_x'] == isod.split('|')[0]]
            sub_isovary = isovar[isovar['tRNA_y'] == isod.split('|')[0]]
            mergevar = pd.concat([sub_isovarx, sub_isovary])
            postest = isod.split('|')[1]
            for row2 in mergevar.iterrows():
                var_positions = row2[1]['Sprinzl Differences'].split(',')
                var_nucleotides = row2[1]['Nucleotide Differences'].split(',')
                var_trnax = row2[1]['tRNA_x']
                var_trnay = row2[1]['tRNA_y']
                if postest in var_positions:
                    mmindex = var_positions.index(str(postest))
                    if var_trnax == isod.split('|')[0]:
                        vartrnaL[isod].append(var_trnay)
                        varnucL[isod].append(var_nucleotides[mmindex])
                    else:
                        vartrnaL[isod].append(var_trnax)
                        varnucL[isod].append(var_nucleotides[mmindex][::-1])
                else:
                    pass
        vartrnaL = {k: list(set(v)) for k, v in vartrnaL.items()}
        varnucL = {k: list(set(v)) for k, v in varnucL.items()}
        cnf['Potential Multimapping tRNAs'] = cnf['varbarcode'].map(vartrnaL)
        pmt = [','.join(x) for x in list(cnf['Potential Multimapping tRNAs'])]
        cnf['Potential Multimapping tRNAs'] = pmt
        cnf['Base Alterations'] = cnf['varbarcode'].map(varnucL)
        ba = [','.join(x) for x in list(cnf['Base Alterations'])]
        cnf['Base Alterations'] = ba
        cnf['% Multimap'] = 100*(cnf['multimapcoverage'] / cnf['readcoverage'])
        cnf['MultiDiff'] = abs(cnf['% Mismatch'] - cnf['% Multimap'])
        cnf['Filter'] = np.where(((cnf['MultiDiff'] < args.multimapcutoff) & (cnf['Potential Multimapping tRNAs'].str.contains('tRNA-'))) | ((cnf['% Multimap'] >= args.multimapcutoff) & (cnf['Potential Multimapping tRNAs'].str.contains('tRNA-'))), 'Potential Artifact', 'Pass')
        cnf = cnf.drop(['varbarcode', 'MultiDiff'], axis=1)
        cnf = cnf[['Feature', 'Sample', 'position', 'actualbase','adenines', 'thymines', 'cytosines','guanines', 'deletions','readcoverage', 'uniquecoverage', 'multimapcoverage', '% Multimap', 'mismatchcounts', '% Mismatch','confidence', 'Potential Multimapping tRNAs', 'Base Alterations', 'Filter']]

        return cnf

    def betacalc(val1, val2):
        return 1 - sc.betainc(val1 + args.alpha, val2 + args.beta, args.editfrac)

    def filterSignificant(row):
        ## Filters for significant differences in the mismatch data
        if row['pval'] <= 0.05 and abs(row['percentdiff']) >= 12:
            return 1
        else:
            return 0

    def betafunction(betadf):
        ## Calculates the beta distribution for the confidence values
        betadf['% Mismatch'] = 100*(betadf['mismatchcounts'] / (betadf['mismatchcounts'] + betadf['refcounts']))
        betadf['confidence'] = betadf.apply(lambda row: tMAPgetmismatch.betacalc(row['mismatchcounts'], row['refcounts']), axis=1).fillna(0)
        return betadf
    
    def generateConfidentSites(self, coverage, nucs):

        nucleotideDict = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}

        ## Generates confident sites for each nucleotide of interest
        confidentmismatches = []
        nonconfidentmismatches = []

        for nucleotide in nucs:
            print("Generating confidence values for - " + nucleotide + "'s...")

            ## Filters for nucleotides of interest and filters out positions with only multimapping mismatches (helps remove false positives)
            coverage_nuc = coverage[coverage['actualbase'] == nucleotide]
            coverage_nuc['multimapcoverage'] = coverage_nuc['multitrnacoverage'] + coverage_nuc['multianticodoncoverage'] + coverage_nuc['multiaminocoverage']
            # coverage_nuc = coverage_nuc[coverage_nuc['mismatchedbases'] != coverage_nuc['multimapcoverage']]

            ## Filters out positions with low coverage
            coverage_nuc = coverage_nuc[['Feature', 'Sample', 'position', 'actualbase', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions', 'coverage', 'uniquecoverage', 'multimapcoverage']]
            coverage_nuc['readcoverage'] = coverage_nuc['adenines'] + coverage_nuc['thymines'] + coverage_nuc['cytosines'] + coverage_nuc['guanines'] + coverage_nuc['deletions']
            coverage_nuc = coverage_nuc[coverage_nuc['readcoverage'] >= int(args.minreads)]
            coverage_nuc['sf'] = coverage_nuc['readcoverage'] / coverage_nuc['coverage']
            coverage_nuc['uniquecoverage'] = coverage_nuc['sf'] * coverage_nuc['uniquecoverage']
            coverage_nuc['multimapcoverage'] = coverage_nuc['sf'] * coverage_nuc['multimapcoverage']
            coverage_nuc = coverage_nuc.drop(['coverage', 'sf'], axis=1)

            ## Adds a column with the barcode for each position
            coverage_nuc['barcode'] = coverage_nuc['Feature'] + '|' + coverage_nuc['position'] + '|' + coverage_nuc['actualbase'] + '|' + coverage_nuc['Sample']
            coverage_nuc = coverage_nuc.set_index('barcode')

            ## Gets the counts for the reference base and the mismatched bases
            realbasecounts = coverage_nuc[nucleotideDict.get(nucleotide)]
            coverage_nuc = coverage_nuc.drop([nucleotideDict.get(nucleotide)], axis=1)
            mismatchbasecounts = coverage_nuc[coverage_nuc.columns[4:8]].T.sum()

            ## Calculates the confidence values for each position
            betainput = pd.concat([realbasecounts, mismatchbasecounts],axis=1).rename(columns={nucleotideDict.get(nucleotide):'refcounts', 0:'mismatchcounts'})

            betaout = tMAPgetmismatch.betafunction(betainput)

            ## Fixes the confidence values for positions with 100% mismatch
            betaout['confidence'] = np.where((betaout['% Mismatch'] == 100) & (betaout['confidence'] == 0), 1, betaout['confidence'])

            ## Merge the confidence values with the coverage data
            coverage_nuc_beta = coverage_nuc.merge(betaout, left_index=True, right_index=True).dropna()

            ## Filters out positions with low confidence
            mismatch_pos = coverage_nuc_beta[coverage_nuc_beta['confidence'] >= 0.95].rename(columns={'refcounts':nucleotideDict.get(nucleotide)})
            confidentmismatches.append(mismatch_pos)

            ## Grabs the soft mismatches, e.g. those which are just below the confidence threshold
            mismatch_neg_soft = coverage_nuc_beta[coverage_nuc_beta['confidence'] < 0.95].rename(columns={'refcounts':nucleotideDict.get(nucleotide)})
            mismatch_neg_soft = mismatch_neg_soft[mismatch_neg_soft['actualbase'] == nucleotide]
            nonconfidentmismatches.append(mismatch_neg_soft)

        return confidentmismatches, nonconfidentmismatches

    def generateConfidentSitesPredict(self, confsites, nucs):

        prediction_list = []
        for nuc in nucs:
            confsites_input = confsites[confsites['actualbase'] == nuc]

            if nuc == 'A':
                confsites_input_subset = confsites_input[['barcode', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']].set_index('barcode')
                confsites_input_subset = confsites_input_subset.div(confsites_input_subset.sum(axis=1), axis=0) * 100
                class_model = joblib.load('/mnt/d/Users/Aidan/bioinformatics/Collabs/Tools/tMAP/EukModel/out/tMAP-Euk-' + nuc + 'model.pkl')
                y_pred = class_model.predict(confsites_input_subset)
                class_probabilities = class_model.predict_proba(confsites_input_subset).tolist()
                outdf = pd.DataFrame(list(zip(list(confsites_input_subset.index), y_pred)),columns=['barcode','Predicted Modification']).set_index('barcode')
                classD = {cl:z for z,cl in enumerate(list(class_model.classes_))}
                probs = []
                for q,prediction in enumerate(list(outdf.index)):
                    predictionval = class_probabilities[q][classD.get(outdf.loc[prediction]['Predicted Modification'])]
                    probs.append(predictionval)
                outdf['Prediction Probability'] = probs
                prediction_list.append(outdf.reset_index())
                
            elif nuc == 'C':
                confsites_input['Predicted Modification'] = 'm3C'
                outdf = confsites_input[['barcode','Predicted Modification', 'confidence']].rename(columns={'confidence':'Prediction Probability'})
                prediction_list.append(outdf)
                
            elif nuc == 'T':
                confsites_input['Predicted Modification'] = 'acp3U'
                outdf = confsites_input[['barcode','Predicted Modification', 'confidence']].rename(columns={'confidence':'Prediction Probability'})
                prediction_list.append(outdf)
                
            elif nuc == 'G':
                confsites_input_subset = confsites_input[['barcode', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']].set_index('barcode')
                confsites_input_subset = confsites_input_subset.div(confsites_input_subset.sum(axis=1), axis=0) * 100
                class_model = joblib.load('/mnt/d/Users/Aidan/bioinformatics/Collabs/Tools/tMAP/EukModel/out/tMAP-Euk-' + nuc + 'model.pkl')
                y_pred = class_model.predict(confsites_input_subset)
                class_probabilities = class_model.predict_proba(confsites_input_subset).tolist()
                outdf = pd.DataFrame(list(zip(list(confsites_input_subset.index), y_pred)),columns=['barcode','Predicted Modification']).set_index('barcode')
                classD = {cl:z for z,cl in enumerate(list(class_model.classes_))}
                probs = []
                for q,prediction in enumerate(list(outdf.index)):
                    predictionval = class_probabilities[q][classD.get(outdf.loc[prediction]['Predicted Modification'])]
                    probs.append(predictionval)
                outdf['Prediction Probability'] = probs
                prediction_list.append(outdf.reset_index())
                
            else:
                pass
            
        final_predictions = pd.concat(prediction_list)
        final_output = confsites.merge(final_predictions, on=['barcode']).set_index('barcode')

        return final_output

    def getMismatchDifferences(confsites, nonconfsites, pwc):
        
        ## Creates a directory to store the output files
        opath = os.getcwd() + '/mismatchdifferences/'
        isdir = os.path.isdir(opath)
        if isdir == True:
            pass
        else:
            os.mkdir(opath)

        multimappinginfo = confsites[list(confsites.columns)[-3:]].reset_index().rename(columns={'index':'barcode'})
        confsites = confsites[list(confsites.columns)[:-3]]
        
        sampleinfo = pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq'])
        conditionDict = dict(zip(list(sampleinfo['sample']), list(sampleinfo['condition'])))   
        ## Merges the confident and non-confident mismatch dataframes
        confsites = pd.concat([nonconfsites, confsites])
        confsites['Sample'] = confsites['Sample'].map(conditionDict)

        ## Gets the pairwise combinations of conditions
        pairwisecombinationlist = list(zip(pwc['condition1'], pwc['condition2']))
        # widths = [3]*len(pairwisecombinationlist)
        # widths.append(0.5)


        # fig, axes = plt.subplots(1,max(2,len(pairwisecombinationlist) + 1), figsize=(1.5*max(2,len(pairwisecombinationlist) * 2),1.75), gridspec_kw={'width_ratios': widths})
        # fig.subplots_adjust(wspace=0.6)
        for idx, comparison in enumerate(pairwisecombinationlist):
            print('Getting significant misincorporation differences for... ' + comparison[0] + '_v_' + comparison[1])
            df1 = confsites[confsites['Sample'] == comparison[0]]
            df2 = confsites[confsites['Sample'] == comparison[1]]

            ## Gets the positions that overlap between the two conditions and filters out positions with below the confidence threshold, default=0.95
            overlapsites = df1.merge(df2, on=['Feature', 'position', 'actualbase']).dropna()
            overlapsites['barcode'] = overlapsites['Feature'] + '|' + overlapsites['position'] + '|' + overlapsites['actualbase']
            overlapsites['maxconfidence'] = np.select([overlapsites['confidence_x'] >= 0.95 , overlapsites['confidence_y'] >= 0.95], [overlapsites['confidence_x'], overlapsites['confidence_y']], default=0)
            overlapsites = overlapsites[overlapsites['maxconfidence'] >= 0.95]

            df1['barcode'] = df1['Feature'] + '|' + df1['position'] + '|' + df1['actualbase']
            df2['barcode'] = df2['Feature'] + '|' + df2['position'] + '|' + df2['actualbase']

            ## Calculates the log2 fold change and percent difference for each position
            pvalDict = {}
            for base in overlapsites['barcode'].unique():
                df1_base = df1[df1['barcode'] == base]
                df2_base = df2[df2['barcode'] == base]

                l2fc = np.log2(df2_base['% Mismatch'].mean() / df1_base['% Mismatch'].mean())
                percentdiff = df2_base['% Mismatch'].mean() - df1_base['% Mismatch'].mean()
                cond1_cov = df1_base['readcoverage'].sum()
                cond2_cov  = df2_base['readcoverage'].sum()
                cond1_mismatches = df1_base['mismatchcounts'].sum()
                cond2_mismatches = df2_base['mismatchcounts'].sum()
                cond1_mismatchpercent = df1_base['% Mismatch'].mean()
                cond2_mismatchpercent = df2_base['% Mismatch'].mean()
                cond1_confidence = df1_base['confidence'].max()
                cond2_confidence = df2_base['confidence'].max()

                ## Calculates the p-value for each position
                tstat, pval = stats.ttest_ind(df2_base['% Mismatch'], df1_base['% Mismatch'])
                pvalDict[base] = [pval, tstat, l2fc, percentdiff, cond1_cov, cond1_mismatches, cond1_mismatchpercent, cond1_confidence, cond2_cov, cond2_mismatches, cond2_mismatchpercent, cond2_confidence]

            dfpval = pd.DataFrame.from_dict(pvalDict, orient='index').reset_index().rename(columns={'index':'barcode', 0:'pval', 1:'tstat', 2:'l2fc', 3:'percentdiff', 4:comparison[0] + '_cov', 5:comparison[0] + '_mismatches', 6:comparison[0] + '_mismatchpercent', 7:comparison[0] + '_confidence', 8:comparison[1] + '_cov', 9:comparison[1] + '_mismatches', 10:comparison[1] + '_mismatchpercent', 11:comparison[1] + '_confidence'}).fillna(1)
            dfpval[['Feature', 'position', 'referencebase']] = dfpval['barcode'].str.split('|', expand=True)
            dfpval['padj'] = smm.multipletests(dfpval['pval'], method='fdr_bh')[1]
            dfpval['padj'] = dfpval['padj'].fillna(1)
            dfpval = dfpval[list(dfpval.columns)[-4:] + list(dfpval.columns)[1:-4]].sort_values(by=['percentdiff','padj']).set_index('Feature')

            ## Writes the output to a csv file
            dfpval.to_csv(opath + comparison[1] + '_v_' + comparison[0] + '-mismatchdifferences.csv')

        #     ## Plots the scatter plot

        #     tMAPgetmismatches.plotDifferenceScatter(dfpval, comparison[1], comparison[0], axes[idx])

        # tMAPgetmismatches.plotLegend(axes[idx + 1])

        # # plt.tight_layout(w_pad=0.38) 
        # plt.savefig(opath + args.o + '-significantmismatchscatter.pdf', bbox_inches='tight')
        # plt.savefig(opath + args.o + '-significantmismatchscatter.svg', bbox_inches='tight')
        # plt.close()

        return

class tMAPalign(object):

    def __init__(self):

        self.eukpositions = list(['drop1',-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e5','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,'drop2'])
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

        return positions

    def gettRNAs(self, alignments):

        isoacceptorlist = []
        isodecoderlist = []

        for line in open(alignments):

            if line.startswith('tRNA-'):
                isodecoders = line.strip().split()
                isoacceptor = re.sub('(-[0-9]+)', '', isodecoders[0])
                isodecoder = isodecoders[0]
                isoacceptorlist.append(isoacceptor)
                isodecoderlist.append(isodecoder)

        return list(set(isoacceptorlist)), list(set(isodecoderlist))

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
                        nucleotide_possible = ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions']
                        trna_data_mean_pos = trna_data_mean_condition[trna_data_mean_condition['indexbase'] == pos]
                        realbase = trna_data_mean_pos.iloc[0]['actualbase']
                        realbasecounts = float(trna_data_mean_pos[nucD.get(realbase)])
                        nucleotide_possible.remove(nucD.get(realbase))
                        trna_data_mean_pos = trna_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                        conditionslice = trna_data_mean_pos[nucleotide_possible].T
                        conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)

                        if float(conditionslice.max()) >= 1:
                            conditionslice = conditionslice.T
                            bottom = 0
                            for base in list(conditionslice.columns):
                                axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], linewidth=0, width=1, color=nucDcolors.get(base), bottom=bottom)
                                bottom += conditionslice.loc['y'][base]
                            axes[idx][qdx].bar(pos, realbasecounts, linewidth=0, width=1, color='#D3D3D3', bottom=bottom)
                        else:
                            axes[idx][qdx].bar(pos, realbasecounts, linewidth=0, width=1, color='#D3D3D3')

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
                    nucleotide_possible = ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions']
                    trna_data_mean_pos = trna_data_mean_condition[trna_data_mean_condition['indexbase'] == pos]
                    realbase = trna_data_mean_pos.iloc[0]['actualbase']
                    realbasecounts = float(trna_data_mean_pos[nucD.get(realbase)])
                    nucleotide_possible.remove(nucD.get(realbase))
                    trna_data_mean_pos = trna_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                    conditionslice = trna_data_mean_pos[nucleotide_possible].T
                    conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)


                    if float(conditionslice.max()) >= 1:
                        conditionslice = conditionslice.T
                        bottom = 0
                        for base in list(conditionslice.columns):
                            axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], linewidth=0, width=1, color=nucDcolors.get(base), edgecolor=nucDcolors.get(base), bottom=bottom)
                            bottom += conditionslice.loc['y'][base]
                        axes[idx][qdx].bar(pos, realbasecounts, color='#D3D3D3', linewidth=0, width=1, edgecolor='#D3D3D3', bottom=bottom)
                    else:
                        axes[idx][qdx].bar(pos, realbasecounts, color='#D3D3D3', linewidth=0, width=1, edgecolor='#D3D3D3')


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
        plt.savefig(opath + args.o + '-mismatchcoverage-query.jpg', bbox_inches='tight', dpi=600)

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
                        axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=0, width=1, color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability))
                        bottom = trna_data_mean_condition[mappability]
                    else:
                        axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=0, width=1, color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability), bottom=bottom)
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
                            axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=0, width=1, color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability))
                            bottom = trna_data_mean_condition[mappability]
                        else:
                            axes[idx][qdx].bar(trna_data_mean_condition['position'], trna_data_mean_condition[mappability], linewidth=0, width=1, color=mapDcolors.get(mappability), edgecolor=mapDcolors.get(mappability), bottom=bottom)
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

class tMAPbedcoverage(object):

    def __init__(self, opath):

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
            if search('tRNA', featurename):
                featureseq = featureseq + 'CCA'
            else:
                featureseq = featureseq

            for j,nuc in enumerate(featureseq):
                lociList.append((featurename.strip(), j + 1, nuc))

        self.locidf = pd.DataFrame.from_records(lociList, columns=['Feature', 'position', 'actualbase'])

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

    def getFeatLengths(d, feats):
        
        featL = []

        for feat in feats:
            d_f = d[d['Feature'] == feat]
            length = d_f['position'].max()
            featL.append(length*0.8)

        return featL

    def plotCoverageInputbed(self, data, opath, sampleinfo):
        conditionDict = dict(zip(list(sampleinfo['sample']), list(sampleinfo['condition']))) 
        conditions = list(sampleinfo['condition'].drop_duplicates())

        sns.set_theme(style="white",font_scale=1)

        print("Generating coverage plots for regions of interest...")
        nucD = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}
        # nucDcolors = {'adenines':'tab:green','cytosines':'tab:blue','thymines':'tab:red','guanines':'tab:orange','deletions':'grey', 'ref':'#D3D3D3'}
        nucDcolors = {'adenines':'#7CAE00','cytosines':'#619CFF','thymines':'#F8766D','guanines':'#CD9600','deletions':'grey', 'ref':'#D3D3D3'}

        featurestoplot = list(data['Feature'].drop_duplicates())

        getfeatlengths = tMAPbedcoverage.getFeatLengths(data, featurestoplot)


        # fig, axes = plt.subplots(max(2, len(featurestoplot)),max(2,len(self.conditions)), sharex=False, figsize=(1.875*max(2,len(self.conditions))*2, max(2, len(featurestoplot))*1.2))
        if max(getfeatlengths) > 100:
            fig, axes = plt.subplots(max(2, len(conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*2,max(2,len(conditions))*0.8), gridspec_kw={'width_ratios': getfeatlengths})
        else:
            fig, axes = plt.subplots(max(2, len(conditions)),max(2,len(featurestoplot)),sharey=False, sharex=False, figsize=(max(2, len(featurestoplot))*2,max(2,len(conditions))*0.8), gridspec_kw={'width_ratios': getfeatlengths})

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
            feat_data['Sample'] = feat_data['Sample'].map(conditionDict)
            feat_data_mean = feat_data.groupby(['Feature', 'Sample', 'position', 'actualbase']).sum().reset_index()

        #     findmax = feat_data_mean['adenines'] + feat_data_mean['thymines'] + feat_data_mean['cytosines'] + feat_data_mean['guanines'] + feat_data_mean['deletions']
        #     ymaxval = round(1.1 * findmax.max())
        
            feat_data_mean['indexbase'] = feat_data_mean['position'].map(sprinzl_x)

            for idx, cond in enumerate(conditions):
                feat_data_mean_condition = feat_data_mean[feat_data_mean['Sample'] == cond]

                findmax = feat_data_mean_condition['adenines'] + feat_data_mean_condition['thymines'] + feat_data_mean_condition['cytosines'] + feat_data_mean_condition['guanines'] + feat_data_mean_condition['deletions']
                ymaxval = max(20,round(1.1 * findmax.max()))

                for pos in xorder:
                    nucleotide_possible = ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions']
                    feat_data_mean_pos = feat_data_mean_condition[feat_data_mean_condition['indexbase'] == pos]
                    realbase = feat_data_mean_pos.iloc[0]['actualbase']
                    realbasecounts = float(feat_data_mean_pos[nucD.get(realbase)])
                    nucleotide_possible.remove(nucD.get(realbase))
                    feat_data_mean_pos = feat_data_mean_pos.drop([nucD.get(realbase)], axis=1)
                    conditionslice = feat_data_mean_pos[nucleotide_possible].T
                    conditionslice = conditionslice.rename(columns={conditionslice.columns[0]:'y'}).sort_values(by=['y'], ascending=False)

                    if float(conditionslice.max()) >= 1:
                        conditionslice = conditionslice.T
                        bottom = 0
                        for base in list(conditionslice.columns):
                            axes[idx][qdx].bar(pos, conditionslice.loc['y'][base], width=1, linewidth=0,color=nucDcolors.get(base), bottom=bottom)
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

    def define_parser():
        
        parser = argparse.ArgumentParser(os.path.basename(sys.argv[0]),formatter_class=argparse.RawDescriptionHelpFormatter,)
        subparsers = parser.add_subparsers(dest='module', required=True, description='tMAP modules')
        align_parser = subparsers.add_parser("align", help='Gets the base-level variation exising amongst isoacceptor groups (ex. all Arg-TCT sequence variation), and the pairwise sequence differences between tRNA isodecoders')
        mismatch_parser = subparsers.add_parser("getmismatch", help='Used identify potential sites of base modification present in tRNA sequencing data')
        coverage_parser = subparsers.add_parser("plotcoverage", help='Plots the coverage of tRNA sequencing data across the tRNA transcript to visualize misincorporations')
        bedregion_parser = subparsers.add_parser("bedregion", help='Generates a tRAX-like coverage output file of user defined regions provided in bed format')
        
        align_parser.add_argument('--stk', required=True, help='stk file of mature tRNA alignments')
        align_parser.add_argument('--o', required=True, help='path for output files')
        align_parser.add_argument('--minreads', required=False, help='read coverage cutoff for tRNAs', default=0)
        align_parser.add_argument('--org', required=False, help='organism of interest; euk, bact, arch, mito; default=euk', default='euk')

        mismatch_parser.add_argument('--cov', required=True, help='-coverage.txt file from tRAX output')
        mismatch_parser.add_argument('--o', required=True, help='path for output files')
        mismatch_parser.add_argument('--positions', required=False, help='tRNA positions of interest', nargs='+', default=[])
        mismatch_parser.add_argument('--exppairs', required=False, help='pairs.txt file from tRAX to perform pairwise comaprisons')
        mismatch_parser.add_argument('--samples', required=False, help='samples.txt file from tRAX to group samples for pairwise comaprisons')
        mismatch_parser.add_argument('--org', required=False, help='euk, bact, arch, or mito', default='euk')
        mismatch_parser.add_argument('--alpha', required=False, help='pseudocounts to add to mismatch counts for beta function; default=0', default=0)
        mismatch_parser.add_argument('--beta', required=False, help='pseudocounts to add to reference counts for beta function; default=0', default=0)
        mismatch_parser.add_argument('--editfrac', required=False, help='minimum fraction of mismatched reads for beta function; default=0.05', default=0.05)
        mismatch_parser.add_argument('--minreads', required=False, help='read coverage cutoff for tRNAs; default=10', default=5)
        mismatch_parser.add_argument('--multimapcutoff', required=False, help='Percent of multimapping reads covering a site to mark as artefact', default=20)
        mismatch_parser.add_argument('--predict', required=False, help='whether or not to calculate significant mismatch differences', action='store_true')

        coverage_parser.add_argument('--cov', required=True, help='-coverage.txt file from tRAX output')
        coverage_parser.add_argument('--samples', required=True, help='samples.txt file from tRAX to group samples')
        coverage_parser.add_argument('--mode', required=True, help='mismatch or mapping; if you want to color the base coverage based on nucleotide misincorporations or read mappability',default='mismatch')
        coverage_parser.add_argument('--o', required=True, help='path for output files')
        coverage_parser.add_argument('--trnas', required=False, help='tRNA isodecoders of interest', nargs='+', default=[])
        coverage_parser.add_argument('--org', required=False, help='euk, bact, arch, or mito', default='euk')

        bedregion_parser.add_argument('--fasta', required=True, help='Genome fasta file; needs an index file (.fai), can be generated with samtools faidx fasta.fa')
        bedregion_parser.add_argument('--regions', required=True, help='bed file containing genomic region of interest')
        bedregion_parser.add_argument('--alignments', required=True, help='bam file containing read alignments; needs index files (.fai), can be generated with samtools index alignments.bam', nargs='+')
        bedregion_parser.add_argument('--o', required=True, help='path for output file')
        bedregion_parser.add_argument('--threads', required=False, help='number of threads to process the data on; Default=1', default=1)
        bedregion_parser.add_argument('--clipping', required=False,help='number of bases a read can extend beyond the feature start/end to be included in the analysis; default=2', default=2, type=int)
        bedregion_parser.add_argument('--samples', required=False, help='samples.txt file from tRAX to group samples; required if using --plot', default='')
        bedregion_parser.add_argument('--plot', required=False, help='whether or not to generate coverage plots for each of the features', action='store_true')
        bedregion_parser.add_argument('--split', required=False, help='whether or not using BED12 format; for transcripts with multiple exons will concatenate them together', action='store_true')

        import textwrap
        parser.epilog = textwrap.dedent(
            f"""\
            tMAP modules:\n
            {align_parser.format_usage()}
            {mismatch_parser.format_usage()}
            {coverage_parser.format_usage()}
            {bedregion_parser.format_usage()}
            """
        )
        return parser
    
    parser = define_parser()
    try:
        args = parser.parse_args()
    except:
        # pass
        parser.print_help()
        sys.exit(0)

    def tMAPalignstk():
        alignstk = tMAPalign()

        startTime = datetime.now()

        positions = alignstk.readSprinzl()


        ## Reads in coverage file and generates a list of isoacceptors and isodecoders
        isoacceptorlist, isodecoderlist = alignstk.gettRNAs(args.stk)

        ## Checks for output directory and creates if not present
        opath = os.getcwd() + '/align/'
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
            diff_tuple, diff_df = alignstk.isodecoderAlign(isoaccept, args.stk, positions)
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
            list_of_isodecoder_diffs.append(alignstk.compareIsodecoders(isod, isodecoderdf_all))
        flat_list_of_isodecoder_diffs = [item for sublist in list_of_isodecoder_diffs for item in sublist]
        isodecodervariationdf = pd.DataFrame(flat_list_of_isodecoder_diffs, columns=['tRNA_x', 'tRNA_y', 'Sprinzl Differences', '# of Differences', 'Nucleotide Differences'])
        isodecodervariationdf_removedup = (isodecodervariationdf[~isodecodervariationdf.filter(like='tRNA').apply(frozenset, axis=1).duplicated()].reset_index(drop=True))
        isodecodervariationdf_removedup = isodecodervariationdf_removedup.set_index('tRNA_x')
        isodecodervariationdf_removedup.to_csv(opath + args.o + '-Isodecoder_comparisons.csv')

        print('\nTime elasped: ', datetime.now() - startTime)
        return
    
    def tMAPgetmismatches():
        getmismatches = tMAPgetmismatch()
        startTime = datetime.now()

        ## Reads in coverage file and generates a list of isoacceptors and isodecoders
        coveragedata, isodecoderlist, isoacceptorlist = getmismatches.readCoverage(args.cov)

        ## Checks for output directory and creates if not present
        opath = os.getcwd() + '/mismatchconfidence/'
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
            isoacceptcoverage = coveragedata[coveragedata['Isoacceptor'] == isoaccept].dropna()
            diff_tuple, diff_df = getmismatches.isodecoderAlign(isoaccept, isoacceptcoverage)
            list_of_isodecoder_dfs.append(diff_df)
            list_of_diffs.append(diff_tuple)
        isodecoderdf_all = pd.concat(list_of_isodecoder_dfs)
        isodecoderdf_all.to_csv(opath + args.o + '-Isodecoder_alignments.csv')
        isodecoder_differences = pd.DataFrame(list_of_diffs).rename(columns={0:'Isoacceptor', 1:'Isodecoder Variation Sprinzl'}).set_index('Isoacceptor')
        isodecoder_differences.to_csv(opath + args.o + '-Isodecoder_variation.csv')

        ## Extracts the differences between isodecoders
        list_of_isodecoder_diffs = []
        for isod in isodecoderlist:
            list_of_isodecoder_diffs.append(getmismatches.compareIsodecoders(isod, isodecoderdf_all))
        flat_list_of_isodecoder_diffs = [item for sublist in list_of_isodecoder_diffs for item in sublist]
        isodecodervariationdf = pd.DataFrame(flat_list_of_isodecoder_diffs, columns=['tRNA_x', 'tRNA_y', 'Sprinzl Differences', '# of Differences', 'Nucleotide Differences'])
        isodecodervariationdf_removedup = (isodecodervariationdf[~isodecodervariationdf.filter(like='tRNA').apply(frozenset, axis=1).duplicated()].reset_index(drop=True))
        isodecodervariationdf_removedup = isodecodervariationdf_removedup.set_index('tRNA_x')
        isodecodervariationdf_removedup.to_csv(opath + args.o + '-Isodecoder_comparisons.csv')

        confidentsites, nonconfidentsites = getmismatches.generateConfidentSites(coveragedata, ['A', 'C', 'G', 'T'])

        ## Check isodecoder variation for potential false positives
        # conf_varfilt, nonconf_varfilt = tMAPgetmismatches.variationFalsePositives(pd.concat(confidentsites), pd.concat(nonconfidentsites), isodecodervariationdf_removedup.reset_index())

        if len(args.positions) > 0:
            confidentsitesall = pd.concat(confidentsites)
            confidentsitesall  = confidentsitesall[confidentsitesall['position'].isin(args.positions)]
            nonconfidentsitesall = pd.concat(nonconfidentsites)
            nonconfidentsitesall = nonconfidentsitesall[nonconfidentsitesall['position'].isin(args.positions)]
            conf_varfilt = getmismatches.variationFalsePositives(confidentsitesall, isodecodervariationdf_removedup.reset_index())
            conf_varfilt.to_csv(opath + args.o + '-confidentmismatches-query.csv')
            # confidentsitesall.to_csv(opath + args.o + '-confidentmismatches-query.csv')
            nonconfidentsitesall['% Multimap'] = 100*(nonconfidentsitesall['multimapcoverage'] / nonconfidentsitesall['readcoverage'])
            nonconfidentsitesall = nonconfidentsitesall[['Feature', 'Sample', 'position', 'actualbase','adenines', 'thymines', 'cytosines','guanines', 'deletions','readcoverage', 'uniquecoverage', 'multimapcoverage', '% Multimap', 'mismatchcounts', '% Mismatch','confidence']]
            nonconfidentsitesall.to_csv(opath + args.o + '-nonconfidentmismatches-query.csv')
        else:
            confidentsitesall = pd.concat(confidentsites)
            nonconfidentsitesall = pd.concat(nonconfidentsites)
            conf_varfilt = getmismatches.variationFalsePositives(confidentsitesall, isodecodervariationdf_removedup.reset_index())
            conf_varfilt.to_csv(opath + args.o + '-confidentmismatches.csv')
            # confidentsitesall.to_csv(opath + args.o + '-confidentmismatches.csv')
            nonconfidentsitesall['% Multimap'] = 100*(nonconfidentsitesall['multimapcoverage'] / nonconfidentsitesall['readcoverage'])
            nonconfidentsitesall = nonconfidentsitesall[['Feature', 'Sample', 'position', 'actualbase','adenines', 'thymines', 'cytosines','guanines', 'deletions','readcoverage', 'uniquecoverage', 'multimapcoverage', '% Multimap', 'mismatchcounts', '% Mismatch','confidence']]
            nonconfidentsitesall.to_csv(opath + args.o + '-nonconfidentmismatches.csv')

        if args.predict == True and len(args.positions) == 0:
            confidencepredictions = getmismatches.generateConfidentSitesPredict(confidentsitesall.reset_index(), ['A', 'C', 'G', 'T'])
            confidencepredictions.to_csv(opath + args.o + '-confidentmismatches.csv')
        elif args.predict == True and len(args.positions) > 0:
            confidencepredictions = getmismatches.generateConfidentSitesPredict(confidentsitesall.reset_index(), ['A', 'C', 'G', 'T'])
            confidencepredictions.to_csv(opath + args.o + '-confidentmismatches-query.csv')
        else:
            pass

        if args.exppairs != None:
            print('Performing pairwise comparisons based on input samples/pairs file...')
            tMAPgetmismatch.getMismatchDifferences(conf_varfilt, nonconfidentsitesall, pd.read_csv(args.exppairs, delim_whitespace=True, names=['condition1','condition2']))
        else:
            pass

        print('\nTime elasped: ', datetime.now() - startTime)
        return
    
    def tMAPplotcoverage():
        coverageplots = tMAPcoverageplots(pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq']))

        startTime = datetime.now()

        ## Reads in coverage file and generates a list of isoacceptors and isodecoders
        coveragedata, isodecoderlist, isoacceptorlist, isotypelist = coverageplots.readCoverage(args.cov)

        ## Checks for output directory and creates if not present
        opath = os.getcwd() + '/coverageplots/'
        isdir = os.path.isdir(opath)
        if isdir == True:
            pass
        else:
            os.mkdir(opath)

        if len(args.trnas) == 0 and args.mode == 'mismatch':
            ## Generates coverage plots for each isotype
            coverageplots.plotCoverage(coveragedata, isodecoderlist, isotypelist, opath)
        elif len(args.trnas) > 0 and args.mode == 'mismatch':
            ## Generates coverage plots for each input tRNA
            coverageplots.plotCoverageInput(coveragedata, args.trnas, opath)
        elif len(args.trnas) == 0 and args.mode == 'mapping':
            ## Generates coverage plots for each isotype
            coverageplots.plotCoverageMapping(coveragedata, isodecoderlist, isotypelist, opath)
        elif len(args.trnas) > 0 and args.mode == 'mapping':
            ## Generates coverage plots for each input tRNA
            coverageplots.plotCoverageInputMapping(coveragedata, args.trnas, opath)

        print('\nTime elasped: ', datetime.now() - startTime)
        return
    
    def tMAPbedregion():

        startTime = datetime.now()

        opath = os.getcwd() + '/bedcoverage/'
        isdir = os.path.isdir(opath)
        if isdir == True:
            pass
        else:
            os.mkdir(opath)

        lociofinterest = bedtools.BedTool(args.regions)
        fasta = bedtools.BedTool(args.fasta)

        print('Generating fasta sequences for the regions of interest...')
        if args.split == True:
            lociofinterest.sequence(fi=fasta, name=True, s=True, tab=True, split=True, fo=opath + args.o + '-regionsofinterest.bed')
        else:
            lociofinterest.sequence(fi=fasta, name=True, s=True,tab=True, fo=opath + args.o + '-regionsofinterest.bed')

        tMAPbedcov = tMAPbedcoverage(opath)

        # ## Generates tMAP object
        # tMAPbedcov = tMAPbedcoverage(opath, pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq']))

        with ProcessPool(processes=int(args.threads)) as pool:

            # samplebedcoverage = pool.map(tMAPbedcov.getReadCov2, args.alignments)
            samplebedcoverage = pool.map(tMAPbedcov.getReadCovClip, args.alignments)

        samplebedcoverage_all = pd.concat(samplebedcoverage)
        samplebedcoverage_all = samplebedcoverage_all[['Sample','Feature','position', 'actualbase','coverage','mismatchedbases','readstarts', 'readends','adenines', 'cytosines', 'thymines', 'guanines','deletions']]

        # samplebedcoverage_all = samplebedcoverage_all.sort_values(by=['Sample','Feature','chr','position'])
        samplebedcoverage_all.to_csv(opath + args.o + '-bedcoverage.txt', sep='\t', index=False)

        if args.plot == True:
            tMAPbedcov.plotCoverageInputbed(samplebedcoverage_all, opath, pd.read_csv(args.samples, delim_whitespace=True, names=['sample','condition','fastq']))

        print('\nTime elasped: ', datetime.now() - startTime)
        return
    
    if args.module == 'align':
        tMAPalignstk()
    elif args.module == 'getmismatch':
        tMAPgetmismatches()
    elif args.module == 'plotcoverage':
        tMAPplotcoverage()
    elif args.module == 'bedregion':
        tMAPbedregion()