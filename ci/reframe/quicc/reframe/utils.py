"""@file quicc_library.py"""

import json
import re
import pathlib as pl
from collections import defaultdict
import numpy as np

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

def round_time(val, digits):
    """Round value to digits
    """

    p = (digits - 1) - int(np.floor(np.log10(abs(val))))
    rounded =  np.round(val + 10**-p, p)
    return rounded

def extract_timings(in_file, arch, filter=None):
    """Extract timings from performance report
    """

    # Find start and end of report
    regexp_start = re.compile("^PERFORMANCE REPORT*")
    regexp_end = re.compile("^Log file(s) saved*")

    # Indentify test
    regexp_test = re.compile("^\[(\w+) \%\$testId=(\w+) \%\$splitting=(\w+) .*")

    # Get timings
    if filter == 'Avg':
        regexp_time = re.compile(" *- (\w+Avg): (.*) s \(.*$")
    elif filter == 'Min':
        regexp_time = re.compile(" *- (\w+Min): (.*) s \(.*$")
    elif filter == 'Max':
        regexp_time = re.compile(" *- (\w+Max): (.*) s \(.*$")
    elif filter is None:
        regexp_time = re.compile(" *- (\w+): (.*) s \(.*$")

    db = nested_dict(4, dict)
    is_report = False
    # Scan input file
    with open(in_file) as file:
        for line in file:
            # Check start of performance report
            if not is_report and re.search(regexp_start, line):
                is_report = True
            if is_report and re.search(regexp_end, line):
                is_report = False

            # Process report
            if is_report:
                r = re.search(regexp_test, line)
                # Get test parameters setup
                if r:
                    t = r.group(1)
                    id = r.group(2)
                    split = r.group(3)
                else:
                    # Get timings
                    r = re.search(regexp_time, line)
                    if r:
                        region = r.group(1)
                        time = float(r.group(2))
                        db[t][id][split][arch][region] = round_time(time, 3)

    return db

def write_timings(out_file, db):
    """Write timings dictionary
    """

    # Dump formatted json for dictionary
    with open(out_file, 'w') as file:
        json.dump(db, file, indent=2)

def read_reference_timings(db_file):
    """Read reference timings and build dictionary
    """
    # Read JSON database
    with pl.Path(__file__).parent.joinpath(f'timings/{db_file}') as file:
        db = json.loads(file.read_text())
    return db


def update_timings(new_file, old_file):
    """Update timings if large changes
    """

    # Read new JSON database
    with open(new_file,'r') as file:
        new_db = json.loads(file.read())

    # Read old JSON database
    with open(old_file,'r') as file:
        old_db = json.loads(file.read())

    u_margin = 1.0 + 0.1
    l_margin = 1.0 - 0.25
    # Update entries and add new
    for test,l0 in new_db.items():
        if test not in old_db:
            old_db[test] = l0
        else:
            for id,l1 in l0.items():
                if id not in old_db[test]:
                    old_db[test][id] = l1
                else:
                    for split,l2 in l1.items():
                        if split not in old_db[test][id]:
                            old_db[test][id][split] = l2
                        else:
                            for arch,l3 in l2.items():
                                if arch not in old_db[test][id][split]:
                                    old_db[test][id][split][arch] = l3
                                else:
                                    need_update = False
                                    for region,time in l3.items():
                                        if region not in old_db[test][id][split][arch]:
                                            print((region, old_db[test][id][split][arch]))
                                            old_db[test][id][split][arch][region] = time
                                        else:
                                            old_t = float(old_db[test][id][split][arch][region])
                                            new_t = float(time)
                                            if new_t > old_t*u_margin or new_t < old_t*l_margin:
                                                need_update = True
                                    if need_update:
                                        for region,time in l3.items():
                                            old_db[test][id][split][arch][region] = time
    # Update old timings
    write_timings(old_file, old_db)

def filter_timings(filter, db):
    """Filter timing data"""

    # Crawel through entries and extract requested ones
    timings = dict()
    for name,f in filter.items():
        if f['test'] in db:
            l0 = db[f['test']]
            t = []
            if f['id'] in l0:
                l1 = l0[f['id']]
                for split,l2 in l1.items():
                    if f['arch'] in l2:
                        l3 = l2[f['arch']]
                        if f['region'] in l3:
                            t.append([float(split), float(l3[f['region']])])
            timings[name] = np.array(t)

    return timings

# Custom ticks label function
def fractions(x,pos):
    if x>=1:
        # x is an integer, so just return that
        return '${:.0f}$'.format(x)
    else:
        # find the fractional power of 2
        pow = 2
        for i in range(npow):
            if (x*pow > 1-1e-9):
                break
            pow *= 2
        # this returns a latex formatted fraction
        return '$\\frac{{{:2.0f}}}{{{:2.0f}}}$'.format(x*pow,pow)

def plot_timings(title, filter_cpu, filter_gpu, data_cpu, data_gpu, fname = None):
    """Plot the sweep timings"""

    dpi = 300

    t_cpu = filter_timings(filter_cpu, data_cpu)
    t_gpu = filter_timings(filter_gpu, data_gpu)

    import matplotlib as mpl
    import matplotlib.pylab as plt

    ref = [6.4, 4.8]
    plt.rcParams['figure.figsize'] = [2*ref[0], 2*ref[1]]

    labeled = []
    fig, ax = plt.subplots()
    ax.set_title(title)
    for lname, d in t_cpu.items():
        h, = ax.plot(d[:,0]/filter_cpu[lname].get('cores',1), d[:,1], lw=2, label=lname)
        labeled.append(h)
    for lname, d in t_gpu.items():
        h, = ax.plot(d[:,0]/filter_gpu[lname].get('cores',1), d[:,1], lw=2, label=lname)
        labeled.append(h)
    if 'baseline' in t_cpu:
        d = t_cpu['baseline']
        nodes = d[:,0]/filter_cpu['baseline'].get('cores',1)
        ax.plot(nodes, (2**1)*d[:,1], 'k-', lw=0.5, label='2X speedup')
        ax.text(nodes[0], (2**1)*d[0,1], '0.5x')
        ax.plot(nodes, (2**-1)*d[:,1], 'k-', lw=0.5, label='2X speedup')
        ax.text(nodes[0], (2**-1)*d[0,1], '2x')
        ax.plot(nodes, (2**-2)*d[:,1], 'k-', lw=0.5, label='4X speedup')
        ax.text(nodes[0], (2**-2)*d[0,1], '4x')
        ax.plot(nodes, (2**-3)*d[:,1], 'k-', lw=0.5, label='8X speedup')
        ax.text(nodes[0], (2**-3)*d[0,1], '8x')
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=10)
    ax.set_xlabel('Nodes')
    ax.set_ylabel('Time [s]')
    ax.xaxis.set_ticks(np.logspace(3, 8, num=6, base=2.0))
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(fractions))
    ax.legend(handles=labeled)
    if fname is not None:
        fig.savefig(fname, dpi=dpi)
    else:
        plt.show()

def make_default_plots(file_cpu, file_gpu, id = '109', save = True):
    # Read CPU database, try json first
    try:
        with open(file_cpu) as file:
            d_cpu = json.loads(file.read())
    except:
        d_cpu = extract_timings(file_cpu, 'broadwell')

    # Read GPU database, try json first
    try:
        with open(file_gpu) as file:
            d_gpu = json.loads(file.read())
    except:
        d_gpu = extract_timings(file_gpu, 'p100')

    fname = None

    # ALegendre projectors
    t_base = 'TransformALegendreTests_Poly_P'
    t_ext = '_projector'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'},
            'KokkosCUDA': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyUnitOperatorAvg'},
            'KokkosCUDA with comm': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'}
            }
    title = 'Timings for associated Legendre backward transform'
    if save:
        fname = 'timings_al_backward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)

    # ALegendre integrator
    t_base = 'TransformALegendreTests_Poly_P'
    t_ext = '_integrator'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'},
            'KokkosCUDA': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyUnitOperatorAvg'},
            'KokkosCUDA with comm': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'}
            }
    title = 'Timings for associated Legendre forward transform'
    if save:
        fname = 'timings_al_forward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)

    # Worland projector
    t_base = 'TransformWorlandTests_Poly_P'
    t_ext = '_projector'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'},
            'KokkosCUDA': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyUnitOperatorAvg'},
            'KokkosCUDA with comm': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'}
            }
    title = 'Timings for Jones-Worland backward transform'
    if save:
        fname = 'timings_jw_backward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)

    # Worland integrator
    t_base = 'TransformWorlandTests_Poly_P'
    t_ext = '_integrator'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'applyOperatorsAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'},
            'KokkosCUDA': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyUnitOperatorAvg'},
            'KokkosCUDA with comm': {'test':f'{t_base}_kokkos_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyOperatorsAvg'}
            }
    title = 'Timings for Jones-Worland forward transform'
    if save:
        fname = 'timings_jw_forward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)

    # Mixed Fourier projector
    t_base = 'TransformFourierTests_Mixed_P'
    t_ext = '_projector'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'transformAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'transformAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGpu (min)': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplMin'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformAvg'},
            'viewGpu with comm (min)': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformMin'},
            'VkFFT': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'VkFFT (min)': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplMin'},
            'VkFFT with comm': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformAvg'},
            'VkFFT with comm (min)': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformMin'}
            }
    title = 'Timings for mixed Fourier backward transform'
    if save:
        fname = 'timings_fourier_mixed_backward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)

    # Mixed Fourier integrator
    t_base = 'TransformFourierTests_Mixed_P'
    t_ext = '_integrator'
    test = f'{t_base}_base_t_{t_ext}'
    f_cpu = {
            'baseline': {'test': f'{t_base}_base_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'transformAvg', 'cores':36},
            'viewCpu': {'test': f'{t_base}_viewCpu_t{t_ext}', 'id':id, 'arch':'broadwell', 'region':'transformAvg', 'cores':36}
            }
    f_gpu = {
            'viewGpu': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'viewGp (min)': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplMin'},
            'viewGpu with comm': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformAvg'},
            'viewGpu with comm (min)': {'test':f'{t_base}_viewGpu_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformMin'},
            'VkFFT': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplAvg'},
            'VkFFT (min)': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'applyImplMin'},
            'VkFFT with comm': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformAvg'},
            'VkFFT with comm (min)': {'test':f'{t_base}_viewGpuVkFFT_t{t_ext}', 'id':id, 'arch':'p100', 'region':'transformMin'}
            }
    title = 'Timings for mixed Fourier backward transform'
    if save:
        fname = 'timings_fourier_mixed_forward.pdf'
    plot_timings(title, filter_cpu = f_cpu, filter_gpu = f_gpu, data_cpu = d_cpu, data_gpu = d_gpu, fname = fname)
