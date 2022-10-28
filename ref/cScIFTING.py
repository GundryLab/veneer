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
    l = re.findall('nG', seq)
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

def makeSpecRpt(numSCMprots, numNSBprots, numSCMpsms, numNSBpsms, numOnePSM ):
    rows = []
    numProts = numSCMprots + numNSBprots
    numPSMs = numSCMpsms + numNSBpsms
#    tot = motifs['NXS'] + motifs['NXC'] + motifs['NXT'] + motifs['NXV']
    rows.append({'col1':'All Proteins', 'col2': numProts, 'col3':'', 'col4':''})
    rows.append({'col1':'Proteins with Consensus Motif', 'col2': numSCMprots, 'col3': "{:.2f}".format((numSCMprots/numProts)*100), 'col4': '% of all proteins'})
    rows.append({'col1':'Non-Specific Binding Proteins', 'col2': numNSBprots, 'col3': "{:.2f}".format((numNSBprots/numProts)*100), 'col4': '% of all proteins'})
    rows.append({'col1':'All PSMs', 'col2': numPSMs, 'col3':'', 'col4':''})
    rows.append({'col1':'PSMs with Consensus Motif', 'col2': numSCMpsms, 'col3': "{:.2f}".format((numSCMpsms/numPSMs)*100), 'col4': '% of all PSMs'})
    rows.append({'col1':'Non-Specific Binding PSMs', 'col2': numNSBpsms, 'col3': "{:.2f}".format((numNSBpsms/numPSMs)*100), 'col4': '% of all PSMs'})
    rows.append({'col1':'Proteins with just one SCM PSM', 'col2': numOnePSM, 'col3':'', 'col4':''})
    df = pd.DataFrame(rows)
    return(df)

def cScIFTING(df):
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
        if MPA == 'Master Protein Accessions' or MPA == '':
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
                    proteins[accession]['Peptides'] = {}
                pepSeq = convertSeq(psm['Annotated Sequence'])
                proteins[accession]['countPSMs'] += 1
                motifs = findSCM(psm['Annotated Sequence'])
                if motifs:
                    psm['hasSCM'] = 1
                    for motif in motifs:
                        motif.upper()
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
            pep['protPctPSMs'] = proteins[accession]['countPSMs'] / totPSMs #proteins[accession]['totPSM']
            pep['protExclusive'] = proteins[accession]['exclusive']
            pep['protPctExclusive'] = proteins[accession]['exclusive'] / proteins[accession]['countPSMs']
            pep['protPSMsSCM'] = proteins[accession]['countSCM']
            pep['protPctPSMsSCM'] = proteins[accession]['countSCM'] / proteins[accession]['countPSMs']
            pep['protSCMonePSM'] = 1 if proteins[accession]['countSCM'] == 1 else 0
            pep['pepPSM'] = proteins[accession]['Peptides'][pepSeq]['countPSMs']
            pep['pepPSMwSCM'] = proteins[accession]['Peptides'][pepSeq]['countSCM']
            pep['pctPepPSMwSCM'] = pep['pepPSMwSCM'] / pep['pepPSM']
            if pep['pepPSMwSCM']:
                pep['hasSCM'] = 1
            else:
                pep['hasSCM'] = 0
            pep['hasNG'] = proteins[accession]['Peptides'][pepSeq]['countNG']
            freshPep = pep.copy()
            if proteins[accession]['countSCM']:
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
        prot['pctExclusive'] = proteins[accession]['exclusive'] / proteins[accession]['countPSMs']
        prot['PSMwSCM'] = proteins[accession]['countSCM']
        prot['pctPSMwSCM'] = proteins[accession]['countSCM'] / proteins[accession]['countPSMs']
        if proteins[accession]['countSCM'] == 1:
            prot['SCMonePSM'] = 1
            totOneSCMPSM += 1
        else:
            prot['SCMonePSM'] = 0
#        prot['SCMonePSM'] = 1 if proteins[accession]['countSCM'] == 1 else 0
        if proteins[accession]['countNG']:
            prot['hasNG'] = 1
        else:
            prot['hasNG'] = 0
        if proteins[accession]['countSCM']:
            scmprots.append(prot)
        else:
            nsbprots.append(prot)

        for psm in proteins[accession]['PSMs']:
            if proteins[accession]['countSCM']:
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
    dfspecificity = makeSpecRpt(len(scmprots), len(nsbprots), len(scmpsms), len(nsbpsms), totOneSCMPSM )
    r = []
#    if reagent['P21163'] > 0:
    r.append( {'Total PSMs': reagent['P21163'], 'Accession': 'P21163', 'Reagent': 'PNGase F'})
#    if reagent['P22629'] > 0:
    r.append( {'Total PSMs': reagent['P22629'], 'Accession': 'P22629', 'Reagent': 'Streptavidin'})
#    if reagent['P00761'] > 0:
    r.append( {'Total PSMs': reagent['P00761'], 'Accession': 'P00761', 'Reagent': 'Trypsin'})
    dfreagent = pd.DataFrame(r)

    return( dfscmprot, dfnsbprot, dfscmpep, dfnsbpep, dfscmpsm, dfnsbpsm, dfmiape, dfreagent, dfmotif, dfspecificity)
