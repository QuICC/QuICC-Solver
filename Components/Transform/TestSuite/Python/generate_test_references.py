#!/bin/env python

import argparse
import numpy as np
import subprocess as sproc
import os
from multiprocessing import Pool, Lock

def countModes(db):
    return 2**(db % 100)

def splitId(db, i):
    return db + (i+1)/10000

def createSplitMeta(fpath, db):
    nModes = countModes(db)
    meta = np.genfromtxt(f'{fpath}/P_id{db}_meta.dat')
    inData = np.genfromtxt(f'{fpath}/P_id{db}_in.dat')

    # Write split files
    for i in range(0, nModes):
        # Write meta data
        with open(f'{fpath}/P_id{splitId(db, i)}_meta.dat', 'w') as f:
            f.write(str(4) + "\n") 
            f.write(str(int(meta[1])) + "\n") 
            f.write(str(int(meta[2])) + "\n") 
            f.write(str(1) + "\n") 
            f.write(str(int(i))) 

        # Write input files
        np.savetxt(f'{fpath}/P_id{splitId(db, i)}_in.dat', inData[:,i])

def combineSplitReference(fpath, db):
    nModes = countModes(db)

    data = []
    for i in range(0, nModes):
        col = np.genfromtxt(f'{fpath}/P_id{splitId(db, i)}_ref.dat')
        data.append(col)

    data = np.array(data)
    np.savetxt(f'{fpath}/P_id{db}_ref.dat', data.T)

def processSplitMeta(id):
    print(f'Processing id = {id}...')
    args = ['./TransformWorlandTests', '[Poly::P:projector]', '--id', f'{id}', '--dumpData']
    sproc.call(args, stdout=sproc.DEVNULL, stderr=sproc.DEVNULL,)

def cleanupSplitFiles(datapath, refpath, db):
    nModes = countModes(db)

    for i in range(0, nModes):
        os.remove(f'{refpath}/P_id{splitId(db, i)}_meta.dat')
        os.remove(f'{refpath}/P_id{splitId(db, i)}_in.dat')
        os.remove(f'{datapath}/P_id{splitId(db, i)}_ref.dat')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--pool", type = int, default = 1, help = "Worker pool size")
    parser.add_argument("--db", type = int, required = True, help = "DB id")
    parser.add_argument("--operator", default = 'Projector', help = "Operator")
    args = parser.parse_args()

    nCpu = args.pool
    db = args.db
    op = args.operator

    refbase = "_refdata/Transform/Worland"
    database = "_data/Transform/Worland"
    op = "Projector"
    refdir = f'{refbase}/{op}'
    datadir = f'{database}/{op}'

    # Create split metadata
    createSplitMeta(refdir, db)

    # Process split metadata
    pool = Pool(nCpu)
    ids = [f'{splitId(db, i)}' for i in range(0, countModes(db))]
    pool.map(processSplitMeta, ids)
    pool.close()

    # Collect split reference
    combineSplitReference(datadir, db)
    
    # Cleanup the generated split files
    cleanupSplitFiles(datadir, refdir, db)
