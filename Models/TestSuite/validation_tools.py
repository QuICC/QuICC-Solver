import os,sys, getopt
import numpy as np

import struct
import math

import colorcodes as cc
_c = cc.Colorcodes()

def compute_ulp(x):
    """Return the value of the least significant bit of a
    float x, such that the first float bigger than x is x+ulp(x).
    Then, given an expected result x and a tolerance of n ulps,
    the result y should be such that abs(y-x) <= n * ulp(x).
    The results from this function will only make sense on platforms
    where native doubles are represented in IEEE 754 binary64 format.

    From official test_math.py. Python 3.9 has now math.ulp
    """

    x = abs(float(x))
    if math.isnan(x) or math.isinf(x):
        return x

    # Find next float up from x.
    n = struct.unpack('<q', struct.pack('<d', x))[0]
    x_next = struct.unpack('<d', struct.pack('<q', n + 1))[0]
    if math.isinf(x_next):
        # Corner case: x was the largest finite float. Then it's
        # not an exact power of two, so we can take the difference
        # between x and the previous float.
        x_prev = struct.unpack('<d', struct.pack('<q', n - 1))[0]
        return x - x_prev
    else:
        return x_next - x


def processArgv(argv):
    ref_dir = None
    data_dir = None

    usage = 'validate_benchmark.py -d <data_dir> -r <ref_dir>'
    try:
        opts, args = getopt.getopt(argv,"hd:r:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-d"):
            data_dir = arg + "/"
        elif opt in ("-r"):
            ref_dir = arg + "/"

    return (ref_dir, data_dir)

def printResult(condition, msg, ntabs = 1):
    """Pretty print test results"""

    res = np.zeros(2, dtype='i8')
    if condition:
        status = (_c.green + b'passed' + _c.reset).decode()
    else:
        status = (_c.red + b'failed' + _c.reset).decode()
        res[1] += 1
    res[0] += 1
    tabs = "\t"*ntabs
    print(tabs + msg.ljust(50) + f': {status}')

    return res

def printSummary(results, rows, reftol = None):
    """Pretty print validation test summary"""
    fail = 0
    tot = 0
    minTol = [0]*len(rows)
    for r in results:
        tot += r[0]
        fail += r[1]
        i = rows.index(r[3][0])
        minTol[i] = max(minTol[i], r[3][1])

    if reftol is None:
        reftol = [0]*len(minTol)
    print("")
    print("Passing tolerances (+10%):")
    e = b"["
    epc = b"["
    newtol = []
    for v, r in zip(minTol,reftol):
        vScaled = int(round(v*1.1 + 0.5))
        rr = r
        if r == 0:
            cs = b""
            ce = b""
            rr = vScaled
            newtol.append(vScaled)
        elif vScaled > r:
            cs = _c.red
            ce = _c.reset
            newtol.append(vScaled)
        else:
            cs = _c.green
            ce = _c.reset
            newtol.append(r)
        e += cs + (f'{vScaled:.0f}').encode() + ce + b", "
        epc += cs + (f'{100*(vScaled-r)/rr:.0f}%').encode() + ce + b", "
    e = e[:-2] + b"]"
    epc = epc[:-2] + b"]"
    print(e.decode())
    print("")
    print("Tolerance changes:")
    print(epc.decode())
    print("")
    print("New tolerances:")
    print(newtol)

    print("")
    if(fail == 0):
        print(f'All benchmark validation tests passed!')
    else:
        t = 'test'
        if fail > 1:
            t += 's'
        tc = _c.red + (f'{fail} benchmark validation {t} failed').encode()  + _c.reset
        msg = tc + (f' out of {tot}').encode()
        print(msg.decode())
        print("Failed tests")
        for r in results:
            if r[1] > 0:
                print("\t" + r[2])

def tableTest(fname, ref_dir, data_dir, tid, tol = 11, usecols = None, max_rows = None, threshold = -1, percol = False, perrow = False, max_firstcol = 0):

    # Validate nusselt number
    extra = ''
    if usecols is not None:
        extra = f' usecols = {usecols}'
    if max_rows is not None:
        extra = f' max_rows = {max_rows}'
    if extra:
        extra = ' (' + extra + ' )'
    print(f'Validating {fname}{extra}')
    checks = np.zeros(2, dtype='i8')

    cond = True

    if cond:
        cond = os.path.exists(ref_dir + fname)
        checks += printResult(cond, 'Checking if reference exists')

    if cond:
        cond = os.path.exists(data_dir + fname)
        checks += printResult(cond, 'Checking if data exists')

    if cond:
        ref = np.genfromtxt(ref_dir + fname, usecols=usecols, max_rows = max_rows)
        ref = np.atleast_2d(ref)
        data = np.genfromtxt(data_dir + fname, usecols=usecols, max_rows = max_rows)
        data = np.atleast_2d(data)
        cond = (ref.shape == data.shape)
        checks += printResult(cond, f'Checking file size match (shape: {ref.shape})')

    if cond:
        max_ulp = 0
        # compute reference ulp
        if percol and perrow:
            ref_max = np.max(np.abs(ref[:,max_firstcol:]))
            def get_ulp(r, idx):
                if idx[1] >= max_firstcol:
                    ulp = compute_ulp(ref_max)
                else:
                    ulp = compute_ulp(r)
                return ulp
        elif percol:
            col_max = np.max(np.abs(ref), axis = 0)
            def get_ulp(r, idx):
                if idx[1] >= max_firstcol:
                    ulp = compute_ulp(col_max[idx[1]])
                else:
                    ulp = compute_ulp(r)
                return ulp
        elif perrow:
            row_max = np.max(np.abs(ref[:,max_firstcol:]), axis = 1)
            def get_ulp(r, idx):
                if idx[1] >= max_firstcol:
                    ulp = compute_ulp(row_max[idx[0]])
                else:
                    ulp = compute_ulp(r)
                return ulp
        else:
            def get_ulp(r, idx):
                return compute_ulp(r)

        # Compute error on data
        for idx, r in np.ndenumerate(ref):
            if r > threshold:
                d = data[idx]
                diff = np.abs(r-d)
                ulp = diff/get_ulp(r,idx)
                if ulp > max_ulp:
                    max_ulp = ulp
                if ulp > tol:
                    print((r.item(), d.item(), diff, ulp))

        cond = (max_ulp < tol)
        if tol > 1e3:
            details = f'(tol: {tol:.3e}, '
        else:
            details = f'(tol: {tol:.0f}, '
        if max_ulp > 1e3:
            msg = f'Checking error tolerance'
            details += f'max ulp: {max_ulp:.3e})'
        else:
            msg = f'Checking error tolerance'
            details += f'max ulp: {max_ulp:.0f})'
        checks += printResult(cond, msg)
        print(f'\t{details}')

    return (1, int(not cond), f'{fname}{extra}', (tid, max_ulp))

def check_setup(fname, ref_dir, data_dir, trigger, lines_to_check):
    checked = lines_to_check
    dlines = []
    with open(data_dir + fname) as f:
        for line in f:
            if line.find(trigger) != -1:
                checked = 0
            if checked < lines_to_check:
                checked += 1
                dlines.append(line)
    rlines = []
    with open(ref_dir + fname) as f:
        for line in f:
            if line.find(trigger) != -1:
                checked = 0
            if checked < lines_to_check:
                checked += 1
                rlines.append(line)

    print(f'Validating setup from {fname}')
    # Check lines were found
    cond = (len(dlines) == lines_to_check and len(rlines) == lines_to_check)
    printResult(cond, 'Checking setup is present')

    # Check lines match
    cond = True
    for d,r in zip(dlines,rlines):
        cond = cond and (d == r)
    printResult(cond, 'Checking setup match: ')
