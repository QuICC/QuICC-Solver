import numpy as np
import pandas as pd
from matplotlib import pyplot as pp
import sys, os.path, time
import xml.etree.ElementTree as ET

from matplotlib.ticker import OldScalarFormatter

from base_representer import BaseRepresenter

class TorqueRepresenter(BaseRepresenter):

    _name_columns = [r'$t$', 'value']
    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def open(self, filename=None):

        if filename == None:
            try:
                argv = sys.argv
                # print argv
                self.filename = argv[1]
            except RuntimeError as e:
                print(e)
                print('Supposed usage: python represent_energies.py filename')
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
                if (os.path.isdir(folderpath + '/' + folder)):
                    try:
                        local_filename = folderpath + '/' + folder + '/' + os.path.basename(self.filename)
                        datatemp = pd.read_csv(local_filename, sep='\t', skiprows=3, names=self.name_columns)
#print(time.ctime(os.path.getmtime(local_filename)))
                        last_modified_time = pd.Timestamp(time.ctime(os.path.getmtime(local_filename)))
                        #print(folder, data)

                        # if older than 2017-06-19 correct the value
                        if last_modified_time < pd.Timestamp('2017-06-19'):

                            # find Ekman number
                            local_param_file = folderpath + '/' + folder + '/parameters.cfg'
                            print(local_param_file)
                            cfg_file = open(local_param_file, 'r')
                            header = cfg_file.readline()
                            root = ET.fromstring(header + '<root>' + cfg_file.read() + '</root>')
                            Ek = float(root[2][0][0].text)
                            #print('Parameter file found', Ek)

                            ri = 0.35 / (1. - 0.35)
                            T = np.sqrt(4 * np.pi / 3) * ri
                            datatemp['value'] = (datatemp['value'] + ri * T * 8 * np.pi / 3. * Ek) * ri

                        data.append(datatemp)

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
        data = self.data
        data[r'$t$'] -= min(data[r'$t$'])
        data = data[abs(data['value']) < 1.0]

        # post-processing on the value, to be removed once the value calculations are correct
        #ri = 0.35 / (1. - 0.35)
        #T = np.sqrt(4 * np.pi / 3) * ri
        #data['value'] = (data['value'] + ri * T * 8 * np.pi / 3. * 1e-5) * ri

        # e# ax = data.plot(x='time', y='value', title=folder_name)
        ax = data.plot(x=r'$t$', y='value', alpha=0.25)

        # set parameters for plotting

        # pp.rcParams['font.size'] = 14

        alpha = 0.05
        data.set_index(r'$t$', inplace=True)
        # forward
        ewm = data['value'].ewm(alpha=alpha, adjust=True)
        m = ewm.agg(['mean', 'std'])
        # backward
        ewm_bwd = data['value'][::-1].ewm(alpha=alpha, adjust=True)
        m_bwd = ewm_bwd.agg(['mean', 'std'])
        print(m[::-1].head())
        m = (m + m_bwd) / 2.
        ax = m['mean'].plot(legend='mean')
        ax.fill_between(m.index, m['mean'] - m['std'], m['mean'] + m['std'], alpha=.1, label='std')

        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        ax.set_ylabel(r'$G$')
        # ax.set_title(folder_name)
        #ax.set_xlabel('t')

        BaseRepresenter.draw(self)



if __name__=="__main__":

    reader = TorqueRepresenter()
    reader.open()
    reader.draw()