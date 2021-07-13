#!/usr/bin/env python
##assign_positions.py
##written 1/15/19
##Groves Dixon

#import modules
import time
import argparse
from sys import exit
import pandas as pd


##############################
###### DEFINE FUNCTIONS ######
##############################


def read_input(infileName):
    print('\nReading in {}...'.format(infileName))
    x = pd.read_csv(infileName, delimiter="\t", low_memory=False, comment="#", names = ['chr', 'pos', 'nM', 'nU', 'sample'])
    # x=x.rename(index=int, columns={0: 'chr', 1: "pos"})
    print('done.\n')
    print(x.head())
    x.chr = x.chr.astype(str)
    return(x)


def assign_promoter(strand, start, stop, promoterSize):
    if strand == '+':
        pStart = start - 1 - promoterSize
        if pStart < 0:
            pStart = 0
        promoterCoords = [pStart, start-1]
    elif strand == '-':
        pStart = stop + 1 + promoterSize
        promoterCoords = [stop+1, pStart]
    return(promoterCoords)



def read_gff(inputGff, parentFeature, blockFeature, addBlocks, parentParseString, blockParseString, promoterSize, viewParentless):
    """
    """
    print("\nReading in GFF {}...".format(inputGff))
    gffDict = {}
    parentList = []
    blocksAssigned = 0
    totalBlocks = 0
    skippedList = []
    noParentList = []
    if addBlocks:
        addBlockDict = {}
        for ab in addBlocks:
            addBlockDict[ab] = 0
    with open(inputGff, 'r') as infile:
        for line in infile:

            #skip comment lines and blank lines
            if line[0]=="#" or not line.strip("\n"):
                continue

            #capture lines indicated by lineIndicator (eg "gene" in 3rd column)
            lineString = line
            line=line.strip("\n").split("\t")
            
            #grab gene components
            scaff = line[0]
            leftBound = int(line[3])
            rightBound = int(line[4])
            strand = line[6]
            description = line[8]

            #record information if this is a parent line
            if line[2]==parentFeature:
                parentId = description.split(parentParseString)[1].split(descriptionDelimiter)[0].strip()
                parentList.append(parentId)
                promoterCoords = assign_promoter(strand, leftBound, rightBound, promoterSize)
                gffDict[parentId] = {
                'scaff' : scaff,
                'leftBound' : int(leftBound), 
                'rightBound' : int(rightBound), 
                'strand' : strand, 
                'description' : description,
                'exonCoords' : [],
                'intronCoords' : [],
                'promoterCoords' : promoterCoords
                }
                if addBlocks:
                    gffDict[parentId]['addBlocks'] = {}
                continue
            
            elif line[2]==blockFeature:
                totalBlocks += 1
                #do some checks on this block line
                if blockParseString not in description:
                    print("WARNING. The given parsing string, {}, was not found for this line:")
                    print(lineString)
                    print("Suggest checking your GFF and make sure if fits with given arguments")
                    print("Skipping this line.")
                    skippedList.append(line)
                    continue
                #set parse the parent ID and check it has been encountered
                parentId = description.split(blockParseString)[1].split(descriptionDelimiter)[0].strip().strip('"')
                try:
                    parentDat = gffDict[parentId]
                except KeyError:
                    noParentList.append(line)
                    continue
                if scaff != parentDat['scaff']:
                    print("WARNING. The scaffold for this block line does not match its parent line:")
                    print(lineString)
                    print('skipping this line')
                    skippedList.append(line)
                    continue
                if strand != parentDat['strand']:
                    print("WARNING. The strand for this block line does not match its parent line:")
                    print(lineString)
                    print('skipping this line')
                    skippedList.append(line)
                    continue
                gffDict[parentId]['exonCoords'].append([leftBound, rightBound])
                blocksAssigned += 1
            elif addBlocks:
                if line[2] in addBlocks:
                    ab = line[2]
                    addBlockDict[ab] += 1
                    try:
                        gffDict[parentId]['addBlocks'][ab].append([leftBound, rightBound])
                    except KeyError:
                        gffDict[parentId]['addBlocks'][ab] = [[leftBound, rightBound]]
            else:
                pass

    
    #print results summary
    print("...\nDone reading GFF.")
    print("\tFound {} total lines indicated with '{}' in the third column".format(totalBlocks, blockFeature))
    print("\t{} of these were skipped because they didn't have a recorded parent with a '{}' tag".format(len(noParentList), parentFeature))
    if len(noParentList) > 0:
        print("\tThese could be features from tRNAs or pseudogenes")
    else:
        print("\tGood. Description lines for all blocks had the expected parentID parsing string")
    print("\tTo see these skipped lines, execute the same command with '-viewParentless True'")
    if viewParentless:
        print("\tHere are the parentless lines:")
        for pl in noParentList:
            print("\t".join(pl))
    if addBlocks:
        print("\n\tAlso found and recorded this many features tagged with additional block strings:")
        print("\tstring:\tfound:")
        for ab in addBlocks:
            print("\t{}\t{}".format(ab, addBlockDict[ab]))
    if len(skippedList) > 0:
        print("\t{} more were skipped because their data didn't match their parents".format(len(skippedList)))
    print("\t---\n\t{} total 'parent' features indicated by '{}' were found".format(len(parentList), parentFeature))
    print("\t{} total block lines indicated by '{}' matched these parents and were recorded".format(blocksAssigned, blockFeature))
    if blocksAssigned==totalBlocks:
        print("\tGood. All blocks were assigned to a parent")
    return(parentList, gffDict)



def assign_intron_coords(parentList, gffDict):
    print("\nDesignating exon and intron boundaries...")
    count = 0
    for pid in parentList:
        count += 1
        exonList = gffDict[pid]['exonCoords']
        exonDf = pd.DataFrame(exonList, columns=['start', 'stop'])
        strand = gffDict[pid]['strand']
        if strand=='+':
            exonDf = exonDf.sort_values('start', ascending=True)
            dodge = 1
        elif strand=='-':
            exonDf = exonDf.sort_values('start', ascending=False)
            exonDf.columns = ['stop', 'start']
            dodge = -1
        else:
            exit('error. Strand notation incorrect')
        exonDf['exon'] = range(1, len(exonDf)+1)
        exonDf = exonDf.set_index('exon')
        gffDict[pid]['exonCoords'] = exonDf

        #now set up intron coords
        starts = list(exonDf['start'])
        stops = list(exonDf['stop'])
        icoords = []


        #assign the start and stop positions of each intron
        #and arrange them so that the left bound is in first column and right in second
        if strand == '+':
            for i in range(len(stops)-1):
                istart = stops[i]+dodge
                istop = starts[i+1]-dodge
                icoords.append([istart, istop])
            intronDf = pd.DataFrame(icoords, columns=['start', 'stop'])
        elif strand == '-':
            for i in range(len(stops)-1):
                istart = stops[i]+dodge
                istop = starts[i+1]-dodge
                icoords.append([istop, istart])
            intronDf = pd.DataFrame(icoords, columns=['stop', 'start'])


        intronDf['intron'] = range(1, len(intronDf)+1)
        intronDf = intronDf.set_index('intron')
        gffDict[pid]['intronCoords'] = intronDf

        
        # if pid=='SPU_000889-tr':
        #     print('----------')
        #     print(gffDict[pid]['description'])
        #     print(gffDict[pid]['strand'])
        #     print()
        #     print(gffDict[pid]['exonCoords'])
        #     print()
        #     print(gffDict[pid]['intronCoords'])


    print('done.')
    return(gffDict)


def format_promoter(pid, promoterCoords, subdat, geneLeftBound, geneRightBound, strand):
    """Function to format lines for sites within promoters"""
    psub = pd.DataFrame(
                subdat.ix[
                    (subdat['pos'] >= promoterCoords[0]) &
                    (subdat['pos'] <= promoterCoords[1])
                ]
            )
    psub['parentGene'] = pid.strip('"')
    psub['feature'] = 'promoter'
    psub['strand'] = strand
    if strand == '+':
        psub['tssDist'] = geneLeftBound - psub['pos']
    elif strand == '-':
        psub['tssDist'] = psub['pos'] - geneRightBound
    else:
        exit('strand error')
    # print('=========')
    # print(strand)
    # print(promoterCoords)
    # print(psub)
    return(psub)


def format_intron_exon(pid, subdat, exonCoords, geneLeftBound, geneRightBound, strand, typeString):
    """Function to format lines for sites within exons or introns"""
    modEdat = pd.DataFrame()
    exonNumbers = exonCoords.index.values
    #get sites within bounds of each exon
    for i in exonNumbers:
        esub = pd.DataFrame(
                    subdat.ix[
                        (subdat['pos']>=exonCoords.ix[i][0]) &
                        (subdat['pos']<=exonCoords.ix[i][1])
                    ]
                )
        esub['parentGene'] = pid.strip('"')
        esub['feature'] = '{}_{}'.format(typeString, i)
        esub['strand'] = strand
        if strand == '+':
            esub['tssDist'] = esub['pos'] - geneLeftBound
        elif strand == '-':
            esub['tssDist'] = geneRightBound - esub['pos']
        modEdat = pd.concat([modEdat, esub])
    # print('=========')
    # print(strand)
    # print(exonCoords)
    # print(exonNumbers)
    # print(esub)
    return(modEdat)


def format_addBlocks(pid, addBlocks, addBlockDict, subdat, geneLeftBound, geneRightBound, strand):
    """Function to deal with the additional blocks.
    These were recorded to allow multiple blocks for each type
    """
    modAdat = pd.DataFrame()
    for ab in addBlocks:
        try:
            addCoordList = addBlockDict[ab]
        except KeyError:
            continue
        for i in range(len(addCoordList)):
            addCoords = addCoordList[i]
            asub = pd.DataFrame(
                    subdat.ix[
                        (subdat['pos']>=addCoords[0]) &
                        (subdat['pos']<=addCoords[1])
                    ]
                )
            asub['parentGene'] = pid.strip('"')
            asub['feature'] = '{}_{}'.format(ab, i+1)
            asub['strand'] = strand
            if strand == '+':
                asub['tssDist'] = asub['pos'] - geneLeftBound
            elif strand == '-':
                asub['tssDist'] = geneRightBound - asub['pos']
            modAdat = pd.concat([modAdat, asub])
    return(modAdat)






def assign_positions(dat, parentList, gffDict, addBlocks, outputName):
    """Loop through parent genes recordered from GFF and build up
    table of sites that fall into each gene's features.

    Sites that fall within multiple features will get a line for each
    feature they fall within.

    """
    print('\nAssigning positions to gene features...')
    count = 0
    catCount = 0
    modDat = pd.DataFrame()
    for pid in parentList:
        count += 1
        if count % 100 == 0:
            print('\t{} out of {} parent genes complete'.format(count, len(parentList)))
        pdat = gffDict[pid]
        chrom = pdat['scaff']
        strand = pdat['strand']
        geneLeftBound = pdat['leftBound']
        geneRightBound = pdat['rightBound']
        promoterCoords = pdat['promoterCoords']
        exonCoords = pdat['exonCoords']
        intronCoords = pdat['intronCoords']

        #designate leftmost and rightmost coordinates for this gene
        if strand == '+':
            leftMost = promoterCoords[0] #the promoter left bound
            rightMost = geneRightBound
        elif strand == '-':
            rightMost = promoterCoords[1] #the promoter right bound
            leftMost = geneLeftBound


        #subset the input table for this gene
        subdat = dat.ix[
            (dat['chr']==chrom) &
            (dat['pos'] >= leftMost) &
            (dat['pos'] <= rightMost)]

        #iterate though features and build up modified table

        #promoter
        pdat = format_promoter(pid, promoterCoords, subdat, geneLeftBound, geneRightBound, strand)
        modDat = pd.concat([modDat, pdat])

        #exons
        edat = format_intron_exon(pid, subdat, exonCoords, geneLeftBound, geneRightBound, strand, 'exon')
        modDat = pd.concat([modDat, edat])

        #introns
        idat = format_intron_exon(pid, subdat, intronCoords, geneLeftBound, geneRightBound, strand, 'intron')
        modDat = pd.concat([modDat, idat])

        #additional blocks
        if addBlocks:
            addBlockDict = gffDict[pid]['addBlocks']
            if len(addBlockDict.keys())>0:
                adat = format_addBlocks(pid, addBlocks, addBlockDict, subdat, geneLeftBound, geneRightBound, strand)
                modDat = pd.concat([modDat, adat])



    #add intergenic as any sites that were not assigned so far
    print('\nCalling intergenic sites as any that are still unassigned...')
    allSites = dat['chr'].map(str) + dat['pos'].map(str)
    if modDat.shape[0]>0:
        assignedSites = modDat['chr'].map(str) + modDat['pos'].map(str)
    else:
        assignedSites = []
    unassigned = ~allSites.isin(assignedSites)
    unSub = pd.DataFrame(dat.ix[unassigned])
    unSub['parentGene'] = 'NA'
    unSub['feature'] = 'intergenic'
    unSub['strand'] = 'NA'
    unSub['tssDist'] = 'NA'
    modDat = pd.concat([modDat, unSub])

    #re-order the columns
    posList = ['chr','pos', 'parentGene', 'feature',  'strand', 'tssDist']
    datList = []
    for cname in list(modDat.columns.values):
        if cname not in posList:
            posList.append(cname)
    outDat = modDat[posList]

    #write out
    print('\nWriting Results to file {}...'.format(outputName))
    outDat.to_csv(outputName, sep='\t', index=False, header = False)



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    This script assigns feature positions for genomic sites based on a gff.

    Input is a table with chromosomes (or scaffolds) in column 1, and positions in column 2 plus a gff.
    Additional columns in input table do not matter.

    Output is a modified table with genomic feature positions.

    Positions that are overlapped by multiple features will appear in multiple lines.

    White space and quotation marks will be stripped from gene IDs if present.

    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input_table', help = 'The the input table to add coordinates to (any tab-delimited table with scaffold and position as first two columns)')
    parser.add_argument('-gff', required = True, dest = 'input_gff', help = 'The the input gff')
    parser.add_argument('-o', required = True, dest = 'output_name', help = 'Name for output bed file')
    parser.add_argument('-parentFeature', required = False, default='mRNA', dest = 'parent', help = 'The type of feature you want to use as blocks for each line in the bedfile (default = exon)')
    parser.add_argument('-blockFeatures', required = False, default='exon', dest = 'block', help = 'The type of feature you want to use as blocks for each line in the bedfile (default = exon)')
    parser.add_argument('-addBlocks', required = False, default=False, dest = 'additional_blocks', help = 'comma-separated list of any additional strings in column three you want to record (eg: three_prime_utr); default = False')
    parser.add_argument('-parentId', required = False, default="ID=transcript:", dest = 'parent_id_string', help = 'The string to use to parse out the parent ID from the parent lines.')
    parser.add_argument('-blockId', required = False, default="Parent=transcript:", dest = 'block_id_string', help = 'The string to use to parse out the parent ID from the block lines. The geneID string will appear between this string and a semicolon. Default is ("Parent=transcript:") as expected for an ensemble gff3')
    parser.add_argument('-promoterSize', required = False, default=1000, dest = 'promoter_size', help = 'Size of promoters in bp (default = 1000)')
    parser.add_argument('-viewParentless', required = False, default=False, dest = 'view_parentless', help = 'Change this to true if you only want to view exons that lack a parent line (default=False)')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    inputGff = args.input_gff
    infileName = args.input_table
    outputName = args.output_name
    parentFeature = args.parent
    blockFeature = args.block
    parentParseString = args.parent_id_string
    blockParseString = args.block_id_string
    descriptionDelimiter = ";"
    promoterSize = int(args.promoter_size)
    addBlocks = args.additional_blocks
    if addBlocks:
        addBlocks = addBlocks.split(',')
    viewParentless = args.view_parentless



    #---- RUN FUNCTIONS ----#
    dat = read_input(infileName)
    parentList, gffDict = read_gff(inputGff, parentFeature, blockFeature, addBlocks, parentParseString, blockParseString, promoterSize, viewParentless)
    assign_intron_coords(parentList, gffDict)
    assign_positions(dat, parentList, gffDict, addBlocks, outputName)



    #return time to run
    print("Done.")
    Time = time.time() - Start_time
    print('\nTime took to run: {}'.format(Time))

