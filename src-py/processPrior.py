import os
import sys
import vbem
import click
import pandas as pd
import numpy as np
import eqtools

sep = os.path.sep

@cli.command()
@click.option('--samples', required=True, help='Specify comma seperated sample Names')
@click.option('--outdir', required=True, help='location for dumping alphas')
def createPrior(sampleDirs, outdir):
    '''
    quantification output file and equivalence class file parsing code, and calling vbem module
    some part from RapClust https://github.com/COMBINE-lab/RapClust/blob/master/bin/RapClust
    input:
        sampleDirs: path to the directory of the samples
    output: 
        NONE
    side effect: 
        writes prior file to outdir
    '''
    # Create a logger object.
    import logging
    logger = logging.getLogger('shoal')

    # Initialize coloredlogs.
    import coloredlogs
    coloredlogs.install(level='DEBUG')

    sep = os.path.sep

    #eqfiles = [sep.join([sd, 'aux_info', 'eq_classes.txt']) for sd in sampleDirs]
    # Read in all of the equivalence classes, summing the counts
    #eqCollection = eqtools.EquivCollection(True)
    #for eqfile in eqfiles: 
    #    eqCollection.fromFile(eqfile)
    #    logger.info("Imported file: {}; # eq = {}".format(eqfile, len(eqCollection.eqClasses)))
    #
    #eqCollection.normalizeAuxWeights()

    logger.info("Reading in the quant files")
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampleDirs]
    readVals = [] 
    quant = None
    for sffile in sffiles: #Adds the values from different samples into one data frame
        print("Importing file: {};".format(sffile))
        if quant is None:
            quant = pd.read_table(sffile)
            quant.set_index('Name', inplace=True)
            readVals.append(quant['NumReads'].values)
        else:
            quant2 = pd.read_table(sffile)
            quant2.set_index('Name', inplace=True)
            quant += quant2
            readVals.append(quant2['NumReads'].values)
    readVals = np.stack(readVals, axis=-1)
    logger.info(readVals)
    logger.info(readVals.shape)
    means = readVals.mean(axis=1)
    stddev = np.sqrt(readVals.var(axis=1))
    factors = np.array([min(1.0, stddev[i]/(means[i] + 1e-5)) for i in xrange(len(means))])
    names = quant.index
    d = pd.DataFrame({'Name' : names, 'mean' : means, 'stddev' : stddev, 'factor' : factors}).set_index('Name')
    logger.info(d)
    d.to_csv(outdir + '/prior.tsv', sep='\t')

 
def VBEMcaller(sampleDirs, outdir, priorFile, weight, out, test):
    '''
    quantification output file and equivalence class file parsing code, and calling vbem module
    some part from RapClust https://github.com/COMBINE-lab/RapClust/blob/master/bin/RapClust
    input:
        sampleDirs: path to the directory of the samples
    output:
        alphas: estimated counts after VBEM optimization
    '''
    # Create a logger object.
    import logging
    logger = logging.getLogger('shoal')

    # Initialize coloredlogs.
    import coloredlogs
    coloredlogs.install(level='DEBUG')

    sep = os.path.sep

    eqfiles = [sep.join([sd, 'aux_info', 'eq_classes.txt']) for sd in sampleDirs]

    eqCollection = eqtools.EquivCollection(True)
    for eqfile in eqfiles: #parsing eqClass file
        eqCollection.fromFile(eqfile)
        logger.info("Imported file: {}; # eq = {}".format(eqfile, len(eqCollection.eqClasses)))

    eqCollection.normalizeAuxWeights()
    
    alphasIn = []
    prior = []
    tnames = []
    firstSamp = True
    numSamp = 0
    eqClasses = {}
    auxsDict = {}
    for eqfile in eqfiles: #parsing eqClass file
        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            print("Importing file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
            if firstSamp:
                for i in range(numTran):
                    tnames.append(ifile.readline().rstrip())
            else:
                for i in range(numTran):
                    ifile.readline()

            for i in range(numEq):
                toks = list(map(float, ifile.readline().rstrip().split('\t')))
                nt = int(toks[0])
                trAux = sorted(zip(toks[1:1+nt], toks[1+nt:-1]))
                tids = tuple([int(x) for x,_ in trAux])
                count = toks[-1]
                auxs = [x*count for _,x in trAux]
                if tids in auxsDict:
                    auxsDict[tids] = ([sum(x) for x in zip(auxs, auxsDict[tids][0])], auxsDict[tids][1]+count)
                else:
                    auxsDict[tids] = (auxs, count)
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count
            firstSamp = False
    for key, val in auxsDict.items():
        auxsDict[key] = [ x/(numSamp*val[1]) for x in val[0] ] #caluclates probabilty of fragment given transcript on all samples

    print("Length of tnames = {}".format(len(tnames)))
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampleDirs]

    readVals = [] 
    quant = None
    for sffile in sffiles: #Adds the values from different samples into one data frame
        print("Importing file: {};".format(sffile))
        if quant is None:
            quant = pd.read_table(sffile)
            quant.set_index('Name', inplace=True)
            length = quant['Length'].to_dict()
            readVals.append(quant['NumReads'].values)
        else:
            quant2 = pd.read_table(sffile)
            quant2.set_index('Name', inplace=True)
            quant += quant2
            readVals.append(quant2['NumReads'].values)
    readVals = np.stack(readVals, axis=-1)
    print(readVals)
    print(readVals.shape)
    means = readVals.mean(axis=1)
    var = readVals.var(axis=1)
    factors = np.array([min(1.0, var[i]/(means[i] + 1e-5)) for i in xrange(len(means))])
    names = quant.index
    means = {names[i] : means[i] for i in xrange(len(means))}
    var = {names[i] : var[i] for i in xrange(len(var))}
    factors = {names[i] : factors[i] for i in xrange(len(var))}
        
    alphas  = {}
    txpLength = {}
    for key, val in quant.iterrows():
        alphas[key] = val['NumReads']
        txpLength[key] = val['EffectiveLength']

#    print ("Setting prior")
#    if priorFile == None:
#        for txp in tnames:
#            prior.append(weight *  txpLength[txp])
#    else:
#        with open(priorFile) as pFile:
#            pData = pd.read_table(pFile).set_index('Name')
#            for txp in tnames:
#                prior.append(weight * pData['NumReads'][txp])

    for txp in tnames:
        alphasIn.append(alphas[txp])
        txpLength[txp] = txpLength[txp] / numSamp             #txp length can be different in different samples? Taking average acroos samples

    if test:
        prior = []
        for txp in tnames:
            flatPrior = 1e-3 * txpLength[txp]
            prior.append((factors[txp] * flatPrior) + (1.0 - factors[txp]) * (weight * means[txp]))


        firstSamp = True
        numSamp = 0
        eqClasses = {}
        auxsDict = {}
        for eqfile in eqfiles[-1:]: #parsing eqClass file
            with open(eqfile) as ifile:
                numSamp += 1
                numTran = int(ifile.readline().rstrip())
                numEq = int(ifile.readline().rstrip())
                print("Importing file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
                if firstSamp:
                    for i in range(numTran):
                        tnames.append(ifile.readline().rstrip())
                else:
                    for i in range(numTran):
                        ifile.readline()

                for i in range(numEq):
                    toks = list(map(float, ifile.readline().rstrip().split('\t')))
                    nt = int(toks[0])
                    trAux = sorted(zip(toks[1:1+nt], toks[1+nt:-1]))
                    tids = tuple([int(x) for x,_ in trAux])
                    count = toks[-1]
                    auxs = [x*count for _,x in trAux]
                    if tids in auxsDict:
                        auxsDict[tids] = ([sum(x) for x in zip(auxs, auxsDict[tids][0])], auxsDict[tids][1]+count)
                    else:
                        auxsDict[tids] = (auxs, count)
                    if tids in eqClasses:
                        eqClasses[tids] += count
                    else:
                        eqClasses[tids] = count
                firstSamp = False
        for key, val in auxsDict.items():
            auxsDict[key] = [ x/(numSamp*val[1]) for x in val[0] ] #caluclates probabilty of fragment given transcript on all samples

        quant = pd.read_table(sffiles[-1]).set_index('Name')
        print ("All files imported; Starting VBEM")
        alphasOut = vbem.runVBEM(eqClasses.keys(), np.array(alphasIn), auxsDict, eqClasses, np.array(prior))
        print ("Optimization Complete; Dumping output")
        vbem.dumpSFfile(sep.join([outdir, sampleDirs[-1].split(sep)[-1]+'_{}.sf'.format(weight)]),
                        tnames,
                        quant['Length'].to_dict(),
                        quant['EffectiveLength'].to_dict(),#txpLength,
                        alphasOut)

    elif priorFile == None:
        print("Dumping Prior")
        vbem.dumpSFfile(sep.join([outdir, out]),
                        tnames,
                        length,
                        txpLength,
                        alphasIn)
    else:
        print ("Setting prior")
        with open(priorFile) as pFile:
            pData = pd.read_table(pFile).set_index('Name')
            for txp in tnames:
                prior.append(weight * pData['NumReads'][txp])


        print ("All files imported; Starting VBEM")
        alphasOut = vbem.runVBEM(eqClasses.keys(), np.array(alphasIn), auxsDict, eqClasses, np.array(prior))
        print ("Optimization Complete; Dumping output")
        vbem.dumpSFfile(sep.join([outdir, sampleDirs[-1].split(sep)[-1]+'_{}.sf'.format(weight)]),
                        tnames,
                        length,
                        txpLength,
                        alphasOut)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--samples', required=True, help='Specify comma seperated sample Names')
@click.option('--basepath', required=True, help='basePath of the analysis i.e. directory where yaml resides')
@click.option('--outdir', required=True, help='location for dumping alphas')
@click.option('--prior', required=False, help='Optional: path to .sf file whose alpha to be used as prior else default 1/EffectiveLength')
@click.option('--weight', required=False, type=float, help='necessary when you use prior')
@click.option('--test', is_flag=True)
@click.option('--out', required=False, help='name given to prior file')
def processQuant(samples, basepath, outdir, prior, weight, out, test):
    '''
    Process multisample quantification files
    '''
    # A list of all sample directories
    snames = samples.strip().split(',')
    sampleDirs = [sep.join([basepath, 'salmonData', 'quant', sample]) for sample in snames]

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            print("The output directory already exists, and is not a directory!")
            sys.exit(1)
    else:
        # create it
        os.makedirs(outdir)
    createPrior(sampleDirs, outdir)
    VBEMcaller(sampleDirs, outdir, prior, weight, out, test)

if __name__ == "__main__":
    processQuant()
