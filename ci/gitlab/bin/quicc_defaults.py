#
# Default models for QuICC model tests
#
from typing import NamedTuple

class config(NamedTuple):
    name: str
    tag: str
    test: bool

def defaultConfigs():
    confs = [
        config('BoussinesqSphereDynamo', 'Explicit', True),
        config('BoussinesqShellDynamo', 'Explicit', True),
        config('BoussinesqShellDynamo', 'Implicit', True),
        config('BoussinesqSphereRTC', 'Explicit', True),
        config('BoussinesqSphereRTC', 'Implicit', True),
        config('BoussinesqShellRTC', 'Explicit', True),
        config('BoussinesqShellRTC', 'Implicit', True),
        config('BoussinesqShellTC', 'Explicit', True),
        config('BoussinesqSphereTC', 'Explicit', True),
        ]

    return confs

def defaultModels():
    l = list()
    for c in defaultConfigs():
        l.append(c.name+c.tag)

    return l
