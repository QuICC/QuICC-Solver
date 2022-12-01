#
# Default models for QuICC model tests
#
from typing import NamedTuple, List, Dict

class config(NamedTuple):
    name: str
    tag: str
    test: bool = False
    tasks: int = 1
    variant: str = ''

    def variantextension(self):
        """Combine variant prepended by underscore or nothing"""
        if not self.variant:
            return ''
        else:
            return '_' + self.variant
    def fulltag(self):
        """Combine tag and variant"""
        if not self.variant:
            return self.tag
        else:
            return self.tag + self.variantextension()
    def fullname(self):
        """Combine model name and tag into full name"""
        return self.name + self.fulltag()
    def shortname(self):
        """Combine model name and tag"""
        return self.name + self.tag

class variant(NamedTuple):
    tag: str = 'default'
    tasks: int = -1

default_variants = {
                    'serial' : variant('serial', 1),
                    'mpi' : variant('tubular', 4),
                    'kk' : variant('serial', 1),
                    'kkgpu' : variant('serial', 1),
                    'perf' : variant('', -1),
                        }

# [model][model_tag][pipe_tag] list(variants)
configurations = {  'BoussinesqSphereDynamo': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant()],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        },
                        'Implicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single2d', 4)],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        }
                    },
                    'BoussinesqShellDynamo': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant()],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        },
                        'Implicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single2d', 4)],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant('none', -1)]
                        }
                    },
                    'BoussinesqSphereRTC': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single1d', 6), variant('single2d', 6), variant('tubular', 6) ],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        },
                        'Implicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single1d', 6), variant('single2d', 6), variant('tubular', 6) ],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        }
                    },
                    'BoussinesqShellRTC': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single1d', 6), variant('single2d', 6), variant('tubular', 6) ],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        },
                        'Implicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant('single2d', 6)],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant('none', -1)]
                        }
                    },
                    'BoussinesqSphereTC': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant() ],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        }
                    },
                    'BoussinesqShellTC': {
                        'Explicit' : {
                            'serial' : [variant()],
                            'mpi' : [variant()],
                            'kk' : [variant()],
                            'kkgpu' : [variant()],
                            'perf' : [variant()]
                        }
                    }
                }

def default_configs(pipeline):
    """return a list of models configurations based on pipeline tag """
    confs = []
    for model in configurations:
        for model_tag in configurations[model]:
            for pipe_tag in configurations[model][model_tag]:
                if pipe_tag == pipeline:
                    for v in configurations[model][model_tag][pipe_tag]:
                        if v.tag == 'none':
                            continue
                        elif v.tag == 'build_only':
                            confs.append(config(model, model_tag,))
                        elif v.tag == 'default':
                            confs.append(config(model, model_tag, True,
                                default_variants[pipe_tag].tasks, default_variants[pipe_tag].tag))
                        else:
                            confs.append(config(model, model_tag, True, v.tasks, v.tag))

    return confs

