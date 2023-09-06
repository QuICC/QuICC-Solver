#
# Base class to write yml files for QuICC
#
import yaml
import math
import sys

def validate_yaml(yml):
    try:
        yaml.safe_load(yml)
        return yml
    except:
        sys.exit('Failed to validate yml.')


"""Base class, defines methods to write yaml file"""
class base_yaml:
    def __init__(self):
        self.actions = [self.empty]

    def empty(self):
        self.config = {}

    """apply actions to build pipeline"""
    def pipe(self):
        for a in self.actions:
            a()

    """hook to overwrite file name"""
    def set_file_name(self):
        self.file_name = '.quicc_base_pipeline'

    def yaml(self):
        self.pipe()
        config_ret = yaml.dump(self.config, default_flow_style=False,
        sort_keys=False, width=math.inf)
        return config_ret

    def write(self):
        self.set_file_name()
        gitlab_yml = open(self.file_name+'.yml','w')
        gitlab_yml.write('# This file is generated, do not modify!\n')
        yml = validate_yaml(self.yaml())
        gitlab_yml.write(yml)
        gitlab_yml.close()
