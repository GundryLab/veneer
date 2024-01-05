#! /usr/bin/python

import pandas as pd
import re

def printLog(log):
    return(log)

def findSCM(seq):
    seq = re.sub('\[|\]', '', seq)
    seq = re.split('\.', seq )[1] +  re.split('\.', seq )[2]
    l = re.findall('n[A|C|D|E|F|G|H|I|K|L|M|N|O|Q|R|S|T|U|V|W|Y|a|c|d|e|f|g|h|i|k|l|m|n|o|q|r|s|t|u|v|w|y][C|S|T|V|c|s|t|v]', seq)
    if l:
        return(l)
    else:
        return(0)

def findDeamination(seq):
    seq = re.sub('\[|\]', '', seq)
    seq = re.split('\.', seq )[1] +  re.split('\.', seq )[2]
    l = re.findall('n[g|G][C|S|T|V|c|s|t|v]', seq)
    if l:
        return(l)
    else:
        return (0)


def convertSeq(seq):
    seq = re.split('\.', seq )[1]
    seq = seq.upper()
    return(seq)

def makeMotifRpt(motifs):
    rows = []
    tot = motifs['S'] + motifs['C'] + motifs['T'] + motifs['V']
    rows.append({'col1':tot, 'col2': 'Total PSMs w/ Consensus Motif'})
    rows.append({'col1': motifs['T'], 'col2': 'PSMs w/ nXT Consensus Motif'})
    rows.append({'col1': "{:.2f}".format((motifs['T']/tot)*100) , 'col2':'% nXT'})
    rows.append({'col1': motifs['S'], 'col2': 'PSMs w/ nXS Consensus Motif'})
    rows.append({'col1': "{:.2f}".format((motifs['S']/tot)*100) , 'col2':'% nXS'})
    rows.append({'col1': motifs['C'], 'col2': 'PSMs w/ nXC Consensus Motif'})
    rows.append({'col1': "{:.2f}".format((motifs['C']/tot)*100) , 'col2':'% nXC'})
    rows.append({'col1': motifs['V'], 'col2': 'PSMs w/ nXV Consensus Motif'})
    rows.append({'col1': "{:.2f}".format((motifs['V']/tot)*100) , 'col2':'% nXV'})
    df = pd.DataFrame(rows)
    return(df)

def makeSpecRpt(numSCMprots, numNSBprots, numSCMpsms, numNSBpsms, numOnePSM, totSCMnG ):
    rows = []
    numProts = numSCMprots + numNSBprots
    numPSMs = numSCMpsms + numNSBpsms
#    tot = motifs['NXS'] + motifs['NXC'] + motifs['NXT'] + motifs['NXV']
    rows.append({'col1':'All Proteins', 'col2': numProts, 'col3':'', 'col4':''})
    rows.append({'col1':'Proteins with SCM', 'col2': numSCMprots, 'col3': "{:.2f}".format((numSCMprots/numProts)*100), 'col4': '% of all proteins'})
    rows.append({'col1':'Non-Specific Binding Proteins', 'col2': numNSBprots, 'col3': "{:.2f}".format((numNSBprots/numProts)*100), 'col4': '% of all proteins'})
    rows.append({'col1':'All PSMs', 'col2': numPSMs, 'col3':'', 'col4':''})
    rows.append({'col1':'PSMs with SCM', 'col2': numSCMpsms, 'col3': "{:.2f}".format((numSCMpsms/numPSMs)*100), 'col4': '% of all PSMs'})
    rows.append({'col1':'Non-Specific Binding PSMs', 'col2': numNSBpsms, 'col3': "{:.2f}".format((numNSBpsms/numPSMs)*100), 'col4': '% of all PSMs'})
    rows.append({'col1':'Sequon Proteins with just one PSM', 'col2': '', 'col3':'', 'col4':''})
    rows.append({'col1':'NSB Proteins with just one PSM', 'col2': '', 'col3':'', 'col4':''})
    rows.append({'col1':'Proteins with just one SCM PSM', 'col2': numOnePSM, 'col3':'', 'col4':''})
    rows.append({'col1':'# PSM with n-G-S/T/C/V in SCM', 'col2': totSCMnG, 'col3':'', 'col4':''})
    df = pd.DataFrame(rows)
    return(df)

def getFilterInfo():
    dict = {}
    with open('ref/filter.csv', 'r') as fo:
        fo.readline() # ignore header
        for line in fo:
            line = line.rstrip()
            (accession, PredSi, SignalP, Phobius, SPC, cirfess, topology) = line.split(',')
            dict[accession] = {}
            dict[accession]['PredSi'] = PredSi
            dict[accession]['SignalP'] = SignalP
            dict[accession]['Phobius'] = Phobius
            dict[accession]['SPC'] = int(SPC)
            dict[accession]['cirfess'] = int(cirfess)
            dict[accession]['topology'] = topology
    fo.close()
    return dict

def filterable(d, accession):
    if accession not in d.keys():
        return 0

    score = 0
    if d[accession]['cirfess'] > 0 or d[accession]['PredSi'] == 'Yes' or d[accession]['SignalP'] == 'Yes' or d[accession]['Phobius'] == 'Yes' or d[accession]['SPC'] >= 3:
        score += 2
    if d[accession]['topology'] == 'Yes':
        score += 3
    return score

    # if d[accession]['cirfess'] > 0:
    #     return 1
    # elif d[accession]['PredSi'] == 'Yes' or d[accession]['PredSi'] == 'Yes' or d[accession]['PredSi'] == 'Yes':
    #     return 1
    # elif d[accession]['SPC'] >= 3:
    #     return 1
    # else:
    #     return 0


def cScIFTING(df):
    filterInfo = getFilterInfo()
    prtc = []
    reagent = {'P21163' : 0, 'P22629' : 0, 'P00761' : 0}
    proteins = {}
    miape = []
    miape_fields = {'Annotated Sequence', 'Modifications', 'Master Protein Accessions', '# Missed Cleavages', 'Charge', 'DeltaScore', 'm/z [Da]', 'XCorr'}

    totPSMs = 0
    totMotif = {'S':0, 'C':0, 'T':0, 'V':0}
    totOneSCMPSM = 0

    for psm in df.to_dict('records'):
#        print(psm.keys())
        MPA = psm['Master Protein Accessions']
        #MIAPE stuff
        m = { key:value for key,value in psm.items() if key in miape_fields}
        miape.append(m)
        #main loop
        if MPA == 'Master Protein Accessions' or MPA == '' or pd.isna(MPA):
#        if MPA == 'Master Protein Accessions' or MPA == '' or MPA == 'NA':
            next
        elif re.search('PRTC', MPA):
            prtc.append(psm)
            next
        elif re.search('P21163', MPA):
            reagent['P21163'] += 1
            next
        elif re.search('P22629', MPA):
            reagent['P22629'] += 1
            next
        elif re.search('P00761', MPA):
            reagent['P00761'] += 1
            next
        else:
            totPSMs += 1
            for accession in MPA.split('; '):  #I should probably split on ; and then strip the spaces
                if accession not in proteins:
                    proteins[accession] = {}
                    proteins[accession]['countPSMs'] = 0
                    proteins[accession]['countSCM'] = 0
                    proteins[accession]['hasSCM'] = 0
                    proteins[accession]['countNG'] = 0
                    proteins[accession]['exclusive'] = 0
                    proteins[accession]['countPeps'] = 0
                    proteins[accession]['PSMs'] = []
                    proteins[accession]['T'] = 0
                    proteins[accession]['S'] = 0
                    proteins[accession]['C'] = 0
                    proteins[accession]['V'] = 0
                    proteins[accession]['Level of Evidence'] = ''
                    proteins[accession]['Confidence Tier'] = ''                    
                    proteins[accession]['Peptides'] = {}
                pepSeq = convertSeq(psm['Annotated Sequence'])
                proteins[accession]['countPSMs'] += 1
                motifs = findSCM(psm['Annotated Sequence'])
                if motifs:
                    psm['hasSCM'] = 1
                    for motif in motifs:
                        motif.upper()
                        proteins[accession][motif[-1].upper()] += 1
                        totMotif[motif[-1].upper()] += 1
                else:
                    psm['hasSCM'] = 0
                proteins[accession]['countSCM'] += psm['hasSCM']
                nG = findDeamination(psm['Annotated Sequence'])
                if nG:
                    psm['hasNG'] = 1
                else:
                    psm['hasNG'] = 0
                proteins[accession]['countNG'] += psm['hasNG']
                if not re.search(';', MPA):
                    proteins[accession]['exclusive'] += 1
                if re.search('-', accession):
                    (noIso, iso) = accession.split('-')
                else:
                    (noIso, iso) = (accession, '')
                proteins[accession]['MPAnoIso'] = noIso
#                proteins[accession]['MPAiso'] = iso
                psm['pepSeq'] = pepSeq
                psm['annSeq'] = psm['Annotated Sequence']
                psm['MPAnonSplit'] = MPA
                psm['numMPA'] = len(MPA.split(';'))
                psm['MPA'] = accession
                psm['MPAnoIso'] = noIso
                # create a new PSM to insert rather than a refeerence
                # the psm will change if there is a compound MPA and if we don't
                # add a new PSM, the PSM will change in all proteins
                freshPSM = psm.copy()
                proteins[accession]['PSMs'].append(freshPSM)
#                proteins[accession]['PSMs'].append(psm)
                if pepSeq not in proteins[accession]['Peptides']:
                    proteins[accession]['Peptides'][pepSeq] = {}
                    proteins[accession]['Peptides'][pepSeq]['origMPA'] = MPA
                    proteins[accession]['Peptides'][pepSeq]['countPSMs'] = 0
                    proteins[accession]['Peptides'][pepSeq]['countSCM'] = 0
                    proteins[accession]['Peptides'][pepSeq]['countNG'] = 0
#                    proteins[accession]['Peptides'][pepSeq]['AnnSeq'] = psm['Annotated Sequence']
                    proteins[accession]['countPeps'] += 1
                proteins[accession]['Peptides'][pepSeq]['countPSMs'] += 1
                proteins[accession]['Peptides'][pepSeq]['countSCM'] += psm['hasSCM']
                proteins[accession]['Peptides'][pepSeq]['countNG'] += psm['hasNG']

    nsbprots = []
    scmprots = []
    nsbpeps = []
    scmpeps = []
    nsbpsms = []
    scmpsms = []

    totSCMnG = 0

    for accession in proteins:
        prot={}
        pepCount = 0
        for pepSeq in proteins[accession]['Peptides']:
            pep={}
            pepCount += 1
            pep['pepSeq'] = pepSeq
            pep['MPA'] = accession
#            pep['AnnSeq'] = proteins[accession]['Peptides'][pepSeq]['AnnSeq']
            pep['MPAnoIso'] = proteins[accession]['MPAnoIso']
#            pep['MPAiso'] = proteins[accession]['MPAiso']
            pep['MPAnonSplit'] = proteins[accession]['Peptides'][pepSeq]['origMPA']
            pep['numMPA'] = len(pep['MPAnonSplit'].split(';'))
            pep['protPSMs'] = proteins[accession]['countPSMs']
#            pep['protPctPSMs'] = proteins[accession]['countPSMs'] / totPSMs #proteins[accession]['totPSM']
            pep['protPctPSMs'] = round( proteins[accession]['Peptides'][pepSeq]['countPSMs'] / proteins[accession]['countPSMs'], 2) #/ totPSMs #proteins[accession]['totPSM']
            pep['protExclusive'] = proteins[accession]['exclusive']
            pep['protPctExclusive'] = round( proteins[accession]['exclusive'] / proteins[accession]['countPSMs'], 2 )
            pep['protPSMsSCM'] = proteins[accession]['countSCM']
            pep['protPctPSMsSCM'] = round( proteins[accession]['countSCM'] / proteins[accession]['countPSMs'], 2)
            pep['protSCMonePSM'] = 1 if proteins[accession]['countSCM'] == 1 else 0
            pep['pepPSM'] = proteins[accession]['Peptides'][pepSeq]['countPSMs']
            pep['pepPSMwSCM'] = proteins[accession]['Peptides'][pepSeq]['countSCM']
            pep['pctPepPSMwSCM'] = round( pep['pepPSMwSCM'] / pep['pepPSM'], 2 )
            if pep['pepPSMwSCM']:
                pep['hasSCM'] = 1
            else:
                pep['hasSCM'] = 0
            pep['countNG'] = proteins[accession]['Peptides'][pepSeq]['countNG']
            freshPep = pep.copy()
            if proteins[accession]['countSCM'] and filterable(filterInfo, proteins[accession]['MPAnoIso']):
                scmpeps.append(freshPep)
            else:
                nsbpeps.append(freshPep)
        prot['MPA'] = accession
        prot['MPAnoIso'] = proteins[accession]['MPAnoIso']
#        prot['MPAiso'] = proteins[accession]['MPAiso']
        prot['numPep'] = pepCount
#        prot['totPSM'] = totPSMs #proteins[accession]['totPSM']
        prot['numPSM'] = proteins[accession]['countPSMs']  #proteins[accession]['totPSM']
#        prot['pctPSM'] = proteins[accession]['countPSMs'] / totPSMs #proteins[accession]['totPSM']
        prot['psmExclusive'] = proteins[accession]['exclusive']
        prot['pctExclusive'] = round( proteins[accession]['exclusive'] / proteins[accession]['countPSMs'], 2 )
        prot['PSMwSCM'] = proteins[accession]['countSCM']
        prot['pctPSMwSCM'] = round( proteins[accession]['countSCM'] / proteins[accession]['countPSMs'], 2 )
        if proteins[accession]['countSCM'] == 1:
            prot['SCMonePSM'] = 1
            totOneSCMPSM += 1
        else:
            prot['SCMonePSM'] = 0
#        prot['SCMonePSM'] = 1 if proteins[accession]['countSCM'] == 1 else 0
        # if proteins[accession]['countNG']:
        #     prot['hasNG'] = 1
        # else:
        #     prot['hasNG'] = 0
        prot['countNG'] = proteins[accession]['countNG']
        numMotifs = proteins[accession]['T'] + proteins[accession]['S'] + proteins[accession]['C'] + proteins[accession]['V']
        if numMotifs == 0 :
            prot['countNXT'] = 0
            prot['pctNXT'] = 0
            prot['countNXS'] = 0
            prot['pctNXS'] = 0
            prot['countNXC'] = 0
            prot['pctNXC'] = 0
            prot['countNXV'] = 0
            prot['pctNXV'] = 0
        else :
            prot['countNXT'] = proteins[accession]['T']
            prot['pctNXT'] = proteins[accession]['T'] / numMotifs
            prot['countNXS'] = proteins[accession]['S']
            prot['pctNXS'] = proteins[accession]['S'] / numMotifs
            prot['countNXC'] = proteins[accession]['C']
            prot['pctNXC'] = proteins[accession]['C'] / numMotifs
            prot['countNXV'] = proteins[accession]['V']
            prot['pctNXV'] = proteins[accession]['V'] / numMotifs

        if proteins[accession]['countSCM'] == 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 2 :
            prot['Level of Evidence'] = '2'
            prot['Confidence Tier'] = 'Low'
        elif proteins[accession]['countSCM'] == 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 3 :
            prot['Level of Evidence'] = '3'
            prot['Confidence Tier'] = 'Low'
        elif proteins[accession]['countSCM'] == 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 5 :
            prot['Level of Evidence'] = '2+3'
            prot['Confidence Tier'] = 'Low'
        elif proteins[accession]['countSCM'] == 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 0 :
            prot['Level of Evidence'] = '0'
            prot['Confidence Tier'] = 'Low'
        elif proteins[accession]['countSCM'] > 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 0 :
            prot['Level of Evidence'] = '1'
            prot['Confidence Tier'] = 'Low'
        elif proteins[accession]['countSCM'] > 0 and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 2 :
            prot['Level of Evidence'] = '1+2'
            prot['Confidence Tier'] = 'Medium'
        elif proteins[accession]['countSCM'] and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 3 :
            prot['Level of Evidence'] = '1+3'
            prot['Confidence Tier'] = 'Medium'
        elif proteins[accession]['countSCM'] and filterable(filterInfo, proteins[accession]['MPAnoIso']) == 5 :
            prot['Level of Evidence'] = '1+2+3'
            prot['Confidence Tier'] = 'High'
        if proteins[accession]['countSCM'] and filterable(filterInfo, proteins[accession]['MPAnoIso']):
            scmprots.append(prot)
            # this is to count the number of nG for all SMC proteins that we report in the
            # output spreadsheet, on the veneer results screen, and in CSC_Log
            totSCMnG += prot['countNG']
        else:
            nsbprots.append(prot)

        for psm in proteins[accession]['PSMs']:
            if proteins[accession]['countSCM'] and filterable(filterInfo, proteins[accession]['MPAnoIso']):
                scmpsms.append(psm)
            else:
                nsbpsms.append(psm)

    dfscmprot = pd.DataFrame(scmprots)
    dfnsbprot = pd.DataFrame(nsbprots)
    dfscmpep = pd.DataFrame(scmpeps)
    dfnsbpep = pd.DataFrame(nsbpeps)
    dfscmpsm = pd.DataFrame(scmpsms)
    dfnsbpsm = pd.DataFrame(nsbpsms)
    dfmiape = pd.DataFrame(miape)
    dfmotif = makeMotifRpt(totMotif)
    dfspecificity = makeSpecRpt(len(scmprots), len(nsbprots), len(scmpsms), len(nsbpsms), totOneSCMPSM, totSCMnG )
    r = []
#    if reagent['P21163'] > 0:
    r.append( {'Total PSMs': reagent['P21163'], 'Accession': 'P21163', 'Reagent': 'PNGase F'})
#    if reagent['P22629'] > 0:
    r.append( {'Total PSMs': reagent['P22629'], 'Accession': 'P22629', 'Reagent': 'Streptavidin'})
#    if reagent['P00761'] > 0:
    r.append( {'Total PSMs': reagent['P00761'], 'Accession': 'P00761', 'Reagent': 'Trypsin'})
    dfreagent = pd.DataFrame(r)

    return( dfscmprot, dfnsbprot, dfscmpep, dfnsbpep, dfscmpsm, dfnsbpsm, dfmiape, dfreagent, dfmotif, dfspecificity )

#xl = pd.read_excel('/home/jack/work/Projects/src/VeneerNG/ref/test.xlsx', engine='openpyxl')
#xl = pd.read_excel('/home/jack/work/visun/Done/vSMC (5)/KB32.xlsx', engine='openpyxl')
#xl = pd.read_excel('/home/jack/work/Projects/veneer/ref/0945_Simone_20200815_100_120_KB32.xlsx', engine='openpyxl')
#output = cScIFTING(xl)
