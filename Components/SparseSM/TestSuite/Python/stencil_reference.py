import scipy.io as io
import quicc.geometry.spherical.shell_radius_boundary as bc

cases = dict()

id = 0
nN = 10
riro = 0.35
ri = riro/(1.0 - riro)
ro = 1.0/(1.0 - riro)
l = 1
cases[id] = [nN, nN, ri, ro, l]

id = 1
nN = 21
l = 6
cases[id] = [nN, nN, ri, ro, l]

id = 10
nN = 10
riro = 0.47
ri = riro/(1.0 - riro)
ro = 1.0/(1.0 - riro)
l = 3
cases[id] = [nN, nN, ri, ro, l]

id = 11
nN = 21
l = 7
cases[id] = [nN, nN, ri, ro, l]


def writeMeta(fbase, opts):
    """ Write meta data file"""

    fname = fbase + '_meta.dat'
    with open(fname, 'w') as f:
        f.write(str(len(opts)) + '\n')
        for v in opts:
            f.write(str(v) + '\n')

def writeMatrixMarket(fbase, mat):
    """ Write meta data file"""

    fname = fbase + '_ref.dat'
    with open(fname, 'wb') as f:
        io.mmwrite(f, mat)

for i in [0, 1, 10, 11]:
    opts = cases[i][:-1]
    nN, _, ri ,ro = opts
    fid = f'_id{i}'
    # Value
    fbase = f'Value' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_value(nN, ri, ro, 0, {'c':None})
    writeMatrixMarket(fbase, S) 

    # Diff
    fbase = f'D1' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_diff(nN, ri, ro, 0, {'c':None})
    writeMatrixMarket(fbase, S) 

    # R1D1DivR1
    fbase = f'R1D1DivR1' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_rdiffdivr(nN, ri, ro, 0, {'c':None})
    writeMatrixMarket(fbase, S) 

    # ValueD1
    fbase = f'ValueD1' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_value_diff(nN, ri, ro, 0, {'c':None})
    writeMatrixMarket(fbase, S) 

    # ValueD2
    fbase = f'ValueD2' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_value_diff2(nN, ri, ro, 0, {'c':None})
    writeMatrixMarket(fbase, S) 

    # L dependent conditions
    opts = cases[i]
    nN, _, ri ,ro, l = opts
    # InsulatingShell
    fbase = f'InsulatingShell' + fid
    writeMeta(fbase, opts) 
    S = bc.stencil_insulating(nN, ri, ro, 0, {'c':None, 'l':l})
    writeMatrixMarket(fbase, S) 
