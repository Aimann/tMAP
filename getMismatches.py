#!/usr/bin/env python3
#Aidan Manning 6.15.23
#Script to determine sites of tRNA modficiation from sequencing data 
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import re
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import scipy.special as sc

class getMismatches(object):

    def __init__(self):
        ## Define some useful global variables
        self.nucleotideDict = {'A':'adenines','C':'cytosines','T':'thymines','G':'guanines'}
        self.isotypes = ['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','iMet','Phe','Pro','Ser','Thr','Trp','Tyr','Val']

        ## Fix Sprinzl numbering for eukaryotes
        self.fixSprinzl = {'-1':'-1','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6','7':'7','8':'8','9':'9','10':'10','11':'11','12':'12','13':'13','14':'14','15':'15','16':'16','17':'17','17a':'17a','18':'18','19':'19','20':'20','20a':'20a','20b':'20b','21':'21','22':'22','23':'23','24':'24','25':'25','26':'26','27':'27','28':'28','29':'29','30':'30','31':'31','32':'32','33':'33','34':'34','35':'35','36':'36','37':'37','38':'38','39':'39','40':'40','41':'41','42':'42','43':'43','44':'44','45':'45','e1':'e11','e2':'e12','e3':'e13','e4':'e14','e5':'e15','e6':'e16','e7':'e17','e8':'e1','e9':'e2','e10':'e3','e11':'e4','e12':'e5','e13':'e27','e14':'e26','e15':'e25','e16':'e24','e17':'e23','e18':'e22','e19':'e21','46':'46','47':'47','48':'48','49':'49','50':'50','51':'51','52':'52','53':'53','54':'54','55':'55','56':'56','57':'57','58':'58','59':'59','60':'60','61':'61','62':'62','63':'63','64':'64','65':'65','66':'66','67':'67','68':'68','69':'69','70':'70','71':'71','72':'72','73':'73','74':'74','75':'75','76':'76'}

    def readCoverage(self, df):
        ## Reads in the tRAX coverage file and returns a dataframe with the coverage data with some extra annotations for downstream analysis
        data = pd.read_csv(df, sep='\t')
        if args.org == 'euk':
            data = data[data['Feature'].str.startswith('tRNA')]
            data['position'] = data['position'].map(self.fixSprinzl)
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

    def betacalc(self, val1, val2, edit_frac):
        return 1 - sc.betainc(val1 + args.alpha, val2 + args.beta, edit_frac)

    def betafunction(self, betadf):
        ## Calculates the beta distribution for the confidence values
        betadf['% Mismatch'] = 100*(betadf['mismatchcounts'] / (betadf['mismatchcounts'] + betadf['refcounts']))
        betadf['confidence'] = betadf.apply(lambda row: getmismatches.betacalc(row['mismatchcounts'], row['refcounts'], args.editfrac), axis=1).fillna(0)
        return betadf
    
    def generateConfidentSites(self, coverage, nucs):
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
            coverage_nuc['sizefactor'] = coverage_nuc['readcoverage'] / coverage_nuc['coverage']
            coverage_nuc['uniquecoverage'] = coverage_nuc['sizefactor'] * coverage_nuc['uniquecoverage']
            coverage_nuc['multimapcoverage'] = coverage_nuc['sizefactor'] * coverage_nuc['multimapcoverage']
            coverage_nuc = coverage_nuc.drop(['coverage', 'sizefactor'], axis=1)

            ## Adds a column with the barcode for each position
            coverage_nuc['barcode'] = coverage_nuc['Feature'] + '|' + coverage_nuc['position'] + '|' + coverage_nuc['actualbase'] + '|' + coverage_nuc['Sample']
            coverage_nuc = coverage_nuc.set_index('barcode')

            ## Gets the counts for the reference base and the mismatched bases
            realbasecounts = coverage_nuc[self.nucleotideDict.get(nucleotide)]
            coverage_nuc = coverage_nuc.drop([self.nucleotideDict.get(nucleotide)], axis=1)
            mismatchbasecounts = coverage_nuc[coverage_nuc.columns[4:8]].T.sum()

            ## Calculates the confidence values for each position
            betainput = pd.concat([realbasecounts, mismatchbasecounts],axis=1).rename(columns={self.nucleotideDict.get(nucleotide):'refcounts', 0:'mismatchcounts'})
            betaout = getmismatches.betafunction(betainput)

            ## Fixes the confidence values for positions with 100% mismatch
            betaout['confidence'] = np.where((betaout['% Mismatch'] == 100) & (betaout['confidence'] == 0), 1, betaout['confidence'])

            ## Merge the confidence values with the coverage data
            coverage_nuc_beta = coverage_nuc.merge(betaout, left_index=True, right_index=True).dropna()

            ## Filters out positions with low confidence
            mismatch_pos = coverage_nuc_beta[coverage_nuc_beta['confidence'] >= 0.95].rename(columns={'refcounts':self.nucleotideDict.get(nucleotide)})
            confidentmismatches.append(mismatch_pos)

            ## Grabs the soft mismatches, e.g. those which are just below the confidence threshold
            mismatch_neg_soft = coverage_nuc_beta[coverage_nuc_beta['confidence'] < 0.95].rename(columns={'refcounts':self.nucleotideDict.get(nucleotide)})
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

    
if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--cov', required=True, help='-coverage.txt file from tRAX output')
    ap.add_argument('--alpha', required=False, help='pseudocounts to add to mismatch counts for beta function; default=0', default=0)
    ap.add_argument('--beta', required=False, help='pseudocounts to add to reference counts for beta function; default=0', default=0)
    ap.add_argument('--editfrac', required=False, help='minimum fraction of mismatched reads for beta function; default=0.1', default=0.05)
    ap.add_argument('--minreads', required=False, help='read coverage cutoff for tRNAs', default=0)
    ap.add_argument('--positions', required=False, help='tRNA positions of interest', nargs='+', default=[])
    ap.add_argument('--org', required=True, help='euk, bact, arch, or mito', default='euk')
    ap.add_argument('--samples', required=False, help='samples.txt file from tRAX to group samples')
    ap.add_argument('--predict', required=False, help='whether or not to calculate significant mismatch differences', action='store_true')
    ap.add_argument('--multimapcutoff', required=False, help='Percent of multimapping reads covering a site to mark as artefact', default=20)
    ap.add_argument('--o', required=True, help='path for output files')
    args = ap.parse_args()
    
    ## Generates mismatch detection object
    getmismatches = getMismatches()

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
    # conf_varfilt, nonconf_varfilt = tMAPconfidence.variationFalsePositives(pd.concat(confidentsites), pd.concat(nonconfidentsites), isodecodervariationdf_removedup.reset_index())

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


    print('\nTime elasped: ', datetime.now() - startTime)