#
# Default models for QuICC model tests
#
from typing import NamedTuple

class config(NamedTuple):
    name: str
    tag: str
    test: bool
    tasks: int = 4
    variant: str = ''
    ignore: list[str] = []

    def fullname(self):
        """Combine model name and tag into full name"""
        return self.name + self.tag + self.variant

def defaultConfigs():
    confs = [
        config('BoussinesqSphereDynamo', 'Explicit', True),
        config('BoussinesqSphereDynamo', 'Implicit', True, 4, '_single2d'),
        config('BoussinesqShellDynamo', 'Explicit', True),
        config('BoussinesqShellDynamo', 'Implicit', True, 4, '_single2d'),
        config('BoussinesqSphereRTC', 'Explicit', True, 6, '_serial', ['mpi']),
        config('BoussinesqSphereRTC', 'Explicit', True, 6, '_single1d', ['serial']),
        config('BoussinesqSphereRTC', 'Explicit', True, 6, '_single2d', ['serial']),
        config('BoussinesqSphereRTC', 'Explicit', True, 6, '_tubular', ['serial']),
        config('BoussinesqSphereRTC', 'Implicit', True, 6, '_serial', ['mpi']),
        config('BoussinesqSphereRTC', 'Implicit', True, 6, '_single1d', ['serial']),
        config('BoussinesqSphereRTC', 'Implicit', True, 6, '_single2d', ['serial']),
        config('BoussinesqSphereRTC', 'Implicit', True, 6, '_tubular', ['serial']),
        config('BoussinesqShellRTC', 'Explicit', True, 6, '_serial', ['mpi']),
        config('BoussinesqShellRTC', 'Explicit', True, 6, '_single1d', ['serial']),
        config('BoussinesqShellRTC', 'Explicit', True, 6, '_single2d', ['serial']),
        config('BoussinesqShellRTC', 'Explicit', True, 6, '_tubular', ['serial']),
        config('BoussinesqShellRTC', 'Implicit', True, 6, '_serial', ['mpi']),
        config('BoussinesqShellRTC', 'Implicit', True, 6, '_single2d', ['serial']),
        config('BoussinesqShellTC', 'Explicit', True),
        config('BoussinesqSphereTC', 'Explicit', True),
        ]

    return confs
