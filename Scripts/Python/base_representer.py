# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:20 2017

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys, os, re
import pandas as pd
from matplotlib import pyplot as pp
import xml.etree.ElementTree as ET
import matplotlib

class BaseRepresenter:

    def __init__(self):
        self.data = None

        # use an incrementing indices so that draw can be used for different plots
        self.idx_draw = 2

        matplotlib.rcParams.update({'font.size': 18})


        pass

    # define the name columns static variable
    _name_columns = None
    @property
    def name_columns(self):
        if self._name_columns == None:
            raise NotImplementedError('Name Columns not implemented in the base class')
        return self._name_columns

    def open(self, filename = None):

        if filename == None:
            try:
                argv = sys.argv
                # print argv
                self.filename = argv[1]
            except RuntimeError as e:
                #print(e)
                #print('Supposed usage: python represent_energies.py filename')
                sys.exit()

        else:
            self.filename = filename

        try:

            folder_name = os.path.relpath(".", "..")
            data = pd.read_csv(self.filename, sep='\t', skiprows=3, names=self.name_columns)

        except IOError as e:
            folderpath = os.path.dirname(self.filename)
            data = []
            for folder in os.listdir(folderpath):
                if (os.path.isdir(folderpath+'/'+folder)):
                    try:
                        datatemp = pd.read_csv(folderpath+'/'+folder + '/' + os.path.basename(self.filename), sep='\t', skiprows=3,
                                               names=self.name_columns)
                        #print('here')
                        #print(folder, data)
                        data.append(datatemp)
                        #print(datatemp)
                        #print(data)
                    except IOError as e:

                        #print(folder + '/' + os.path.basename(self.filename))
                        #print(e)
                        pass

            data = pd.concat(data, ignore_index=True)

            data.reindex()
            # print(data.head())
            data.sort_values([r'$t$'], inplace=True)

        self.data = data

    def draw(self):

        pp.tight_layout()

        # check  for the python pattern at the beginning of the file

        if re.match('python.*', sys.argv[0]):
            pass
        else:
            self.idx_draw = 10000


        try:
            fout = sys.argv[self.idx_draw]
            #print(fout)
            pp.savefig(fout)
            self.idx_draw+=1
        except IndexError as e:
            pp.show()

    def search_in_parameter(self, param):

        dirname = os.path.dirname(self.filename)
        #print(dirname)
        cfg_file = open(dirname + '/parameters.cfg', 'r')
        header = cfg_file.readline()
        root = ET.fromstring(header + '<root>' + cfg_file.read() + '</root>')
        leaf = root.find('*/*/'+param)
        return leaf.text

    def get_data(self):

        return self.data
