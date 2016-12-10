import gzip
import struct
import os
import logging
import logging.handlers
import sys
import json
import pandas as pd

def getBootstraps(quantDir):
    logging.basicConfig(level=logging.INFO)

    bootstrapFile = os.path.sep.join([quantDir, "aux_info", "bootstrap", "bootstraps.gz"])
    nameFile = os.path.sep.join([quantDir, "aux_info", "bootstrap", "names.tsv.gz"])
    if not os.path.isfile(bootstrapFile):
       logging.error("The required bootstrap file {} doesn't appear to exist".format(bootstrapFile))
       sys.exit(1)
    if not os.path.isfile(nameFile):
       logging.error("The required transcript name file {} doesn't appear to exist".format(nameFile))
       sys.exit(1)

    with gzip.open(nameFile) as nf:
        txpNames = nf.read().strip().split('\t')

    ntxp = len(txpNames)
    logging.info("Expecting bootstrap info for {} transcripts".format(ntxp))

    with open(os.path.sep.join([quantDir, "aux_info", "meta_info.json"])) as fh:
        meta_info = json.load(fh)

    if meta_info['samp_type'] == 'gibbs':
        s = struct.Struct('<' + 'd' * ntxp)
    elif meta_info['samp_type'] == 'bootstrap':
        s = struct.Struct('@' + 'd' * ntxp)
    else:
        logging.error("Unknown sampling method: {}".format(meta_info['samp_type']))
        sys.exit(1)

    numBoot = 0
##############################
    sffile = os.path.sep.join([quantDir, "quant.sf"])
    print("Importing file: {};".format(sffile))
    quant = pd.read_table(sffile)
    print("importing bootstraps")
    eLen = quant['EffectiveLength'].values
    tempDf = quant.copy()
##############################
    listBootstrap = []
    # Now, iterate over the bootstrap samples and write each
    with gzip.open(bootstrapFile) as bf:
        denom = sum(tempDf['NumReads'].values / eLen)
        while True:
            try:
                x = s.unpack_from(bf.read(s.size))
                listBootstrap.append(x)
                numBoot += 1
##############################
                if numBoot<1001:
                    logging.info("dumping {} bootstrap".format(numBoot))
                    for index, numRead in enumerate(list(x)):
                        tempDf.set_value(index, 'NumReads', numRead)
                        tempDf.set_value(index, 'TPM', 1000000 * (numRead/eLen[index]) / denom )
                    bootOutFile = os.path.sep.join([quantDir, "aux_info", str(numBoot)+".sf"])
                    tempDf.to_csv(bootOutFile, index=False, sep='\t')
##############################
            except:
                logging.info("read {} all bootstrap values".format(numBoot))
                break

    logging.info("read bootstraps successfully.")
    return listBootstrap
