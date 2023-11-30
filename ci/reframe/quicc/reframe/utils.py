"""@file quicc_library.py"""

import json
import re
import importlib.resources as imlib_res
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

def read_timings(db_file):
    """Read timings and build dictionary
    """

    # Read JSON database
    with imlib_res.files(__package__).joinpath(f'timings/{db_file}').open() as file:
        db = json.loads(file.read())

    return db

def update_timings(new_file, old_file):
    """Update timings if large changes
    """

    # Read new JSON database
    with open(new_file,'r') as file:
        new_db = json.loads(file.read())

    # Read new JSON database
    old_db = read_timings(old_file)

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
