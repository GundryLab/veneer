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

def badPSMs(df):
    uniqIDs = {}
    badIDs = []
    for psm in df.to_dict('records'):
        id = str(psm['First Scan']) + '|' + psm['Spectrum File']
        if id in uniqIDs.keys():
            if convertSeq(psm['Annotated Sequence']) != uniqIDs[id]['AnnSeq'] or psm['Master Protein Accessions'] != uniqIDs[id]['MPA']:
                badIDs.append(id)
            else:
                continue
        else:
            uniqIDs[id]= {}
            uniqIDs[id]['AnnSeq'] = convertSeq( psm['Annotated Sequence'] )
            uniqIDs[id]['MPA'] = psm['Master Protein Accessions']
    return(badIDs)

################################################################################
#
# In many of our experiments, we use two search engines in Proteome Discoverer
# to find PSMs.  This means that one detection event sometimes turns into two
# PSMs in the results.  "Duplicate" PSMs can be identified by having the same
# Scan Number and Spectrum File.  Specifically, we do not wish to double cocunt
# them for the purposes of PSM counts, spontaneous deamidation, etc.  But we do
# wish to report them in the results.
#
# To accomplish this, as we loop through the PSMs in the input file, I keep
# which track of Scan Number/Spectrum File we've seen.  If it is the first time,
# I assign it a value of 1 (and set it in the psm['value'] field) for the
# purpose of counting PSMs/Peptides/Proteins. If I have seen the PSM before,
# then the value is set to zero for the purposes of counting, etc.
#
# Also, some PSMs found by both search engines may disagree on the sequence or
# on the MPA assignment.  If the two search engines cannot agree on both the
# sequence and MPA, thne I throw the diasagreeing PSMs  out.
################################################################################

def cScIFTING(df):
    # If you run this script from python, Pandas reads the data and the "First Scan
    # field has a type of int.  If you read it in R using the readxl library the
    # type is float.  I don't think it matters, but to be safe since I have
    # fewer troubles running this from python/pandas, I convert to int.
    df[["First Scan"]] = df[["First Scan"]].astype(int)
    # now that that is out of the way, let's begin
    prtc = []
    reagent = {'P21163' : 0, 'P22629' : 0, 'P00761' : 0}
    proteins = {}
    miape = []
    miape_fields = {'Annotated Sequence', 'Modifications', 'Master Protein Accessions', '# Missed Cleavages', 'Charge', 'DeltaScore', 'm/z [Da]', 'XCorr'}

    totPSMs = 0
    totMotif = {'S':0, 'C':0, 'T':0, 'V':0}
    totOneSCMPSM = 0
    uniqPSMs = {}

    baddies = badPSMs(df)

    for psm in df.to_dict('records'):
        id = str(psm['First Scan']) + '|' + psm['Spectrum File']
        # remove the bad PSMs from the process
        if id in baddies:
            continue
        # This is where I keep track of the PSM ids to see if it is a Duplicate
        # or not.  Everything runs as normal.  The only difference is the value
        # field.  If it is the first time we have seen the PSM, then the value
        # is 1 and the summing fields (numPSM, numPSMwSCM, etc) tabulate as
        # normal.  If the PSM has been seen before, do not count it in the
        # summing fields by setting the value to 0
        if id in uniqPSMs.keys():
            psm['value'] = 0
        else:
            psm['value'] = 1
            uniqPSMs[id] = 1

        MPA = psm['Master Protein Accessions']
        #MIAPE stuff
        m = { key:value for key,value in psm.items() if key in miape_fields}
        miape.append(m)
        #main loop
        # once again, there is a difference between python/pandas and R/readxl.
        # if reading from python/pandas, the blank MPAs are nan.  From
        # R/readxl, it is empty string - ''
        if MPA == 'Master Protein Accessions' or MPA == '' or pd.isna(MPA):
            continue
        elif re.search('PRTC', MPA):
            prtc.append(psm)
            continue
        elif re.search('P21163', MPA):
            reagent['P21163'] += psm['value']
            continue
        elif re.search('P22629', MPA):
            reagent['P22629'] += psm['value']
            continue
        elif re.search('P00761', MPA):
            reagent['P00761'] += psm['value']
            continue
        else:
            totPSMs += psm['value']
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
                proteins[accession]['countPSMs'] += psm['value']
                motifs = findSCM(psm['Annotated Sequence'])
                if motifs:
                    psm['hasSCM'] = 1
                    for motif in motifs:
                        motif.upper()
                        totMotif[motif[-1].upper()] += psm['value']
                else:
                    psm['hasSCM'] = 0
                proteins[accession]['countSCM'] += psm['hasSCM'] * psm['value']
                nG = findDeamination(psm['Annotated Sequence'])
                if nG:
                    psm['hasNG'] = 1
                else:
                    psm['hasNG'] = 0
                proteins[accession]['countNG'] += psm['hasNG'] * psm['value']
                if not re.search(';', MPA):
                    proteins[accession]['exclusive'] += psm['value']
                if re.search('-', accession):
                    (noIso, iso) = accession.split('-')
                else:
                    (noIso, iso) = (accession, '')
                proteins[accession]['MPAnoIso'] = noIso
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
                if pepSeq not in proteins[accession]['Peptides']:
                    proteins[accession]['Peptides'][pepSeq] = {}
                    proteins[accession]['Peptides'][pepSeq]['origMPA'] = MPA
                    proteins[accession]['Peptides'][pepSeq]['countPSMs'] = 0
                    proteins[accession]['Peptides'][pepSeq]['countSCM'] = 0
                    proteins[accession]['Peptides'][pepSeq]['countNG'] = 0
#                    proteins[accession]['Peptides'][pepSeq]['AnnSeq'] = psm['Annotated Sequence']
                    proteins[accession]['countPeps'] += 1
                # I multiply the hasSCM and hasNG fields by the value field
                # (which is 0 or 1) in order not to overcount the hasSCM and
                # hasNG fields
                proteins[accession]['Peptides'][pepSeq]['countPSMs'] += psm['value']
                proteins[accession]['Peptides'][pepSeq]['countSCM'] += psm['hasSCM'] * psm['value']
                proteins[accession]['Peptides'][pepSeq]['countNG'] += psm['hasNG'] * psm['value']

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
            pep['countNG'] = proteins[accession]['Peptides'][pepSeq]['countNG']
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
        # if proteins[accession]['countNG']:
        #     prot['hasNG'] = 1
        # else:
        #     prot['hasNG'] = 0
        prot['countNG'] = proteins[accession]['countNG']
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
#    return( dfscmprot, dfnsbprot, dfscmpep, dfnsbpep, dfscmpsm, dfnsbpsm, dfmiape, dfreagent, dfmotif, dfspecificity)
#xl = pd.read_excel('/home/jack/work/visun/Done/vSMC (5)/KB32.xlsx', engine='openpyxl')
#xl = pd.read_excel('/home/jack/work/Projects/veneer/ref/sample2.xlsx', engine='openpyxl')
#output = cScIFTING(xl)
