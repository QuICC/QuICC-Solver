import difflib

schemes=["WLFl", "WLFm", "SLFl", "SLFm"]
algorithms=["tubular", "single1d", "single2d"]
fbase = "Data_distribution_"
stages=["TRAB1D", "TRAB2D", "TRAB3D"]
nps=[4]

validated = True
for alg in algorithms:
    print(f'Validating {alg}:')
    for sch in schemes:
        print(f'\t{sch}')
        for np in nps:
            print(f'\t\tNp = {np}')
            for s in stages:
                fname = fbase + s + ".vtp"
                with open(f'{alg}/{sch}/Np{np}/{fname}') as f:
                    data = f.readlines()
                with open(f'../../_refdata/Communicators/{alg}/{sch}/Np{np}/{fname}') as f:
                    ref = f.readlines()
                delta = difflib.unified_diff(data, ref, fromfile=f'ref/{fname}', tofile=f'data/{fname}')
                isSame = (len(list(delta)) == 0)
                validated = validated and isSame
                if(isSame):
                    print(f'\t\t\t{fname} are the same')
                else:
                    print(f'\t\t\t{fname} ARE NOT THE SAME')

if validated:
    print("All tests passed")
else:
    print("Some tests failed")
