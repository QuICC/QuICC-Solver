import numpy as np
from matplotlib import pyplot as pp
import sys, os
import pandas as pd
#from quicc.geometry.spherical import shell_radius
from base_representer import BaseRepresenter
import xml.etree.ElementTree as ET

class SpectraRepresenter(BaseRepresenter):

    _name_columns = ['L spectrum, toroidal', 'M spectrum, toroidal', 'L spectrum, poloidal', 'M spectrum, poloidal']
    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def open(self, filename=None):

        if filename==None:
            try:
                argv = sys.argv
                self.filename = argv[1]
            except RuntimeError as e:
                print(e)
                print('Supposed usage: python represent_spectra.py filename')
                sys.exit()
        else:
            self.filename=filename


        try:
            #folder_name = os.path.relpath("..","../..")
            #folder_name = os.path.relpath(".", "..")
            folderpath = os.path.dirname(self.filename)
            datafull = np.loadtxt(self.filename, skiprows=3)
        except IOError as e:

            folder_name = os.path.relpath(".", "..")
            datafull = []
            for folder in os.listdir(folderpath):
                if (os.path.isdir(folderpath+'/'+folder)):
                    #print(folder)
                    try:
                        #print(folderpath+'/'+folder + '/' + os.path.basename(self.filename))
                        datatemp = pd.DataFrame(np.loadtxt(folderpath+'/'+folder + '/' + os.path.basename(self.filename), skiprows=3))
                        datafull.append(datatemp)
                        #print(datatemp)
                    except IOError as e:
                        #print(e)
                        pass

            datafull = pd.concat(datafull, ignore_index=True)
            datafull.reindex()
            datafull.sort_values([0], inplace=True)

            datafull = datafull.as_matrix()
        self.datafull = datafull

        try:

            dirname = os.path.dirname(self.filename)
            print(dirname+ '/parameters.cfg')
            cfg_file = open(dirname + '/parameters.cfg', 'r')
            header = cfg_file.readline()
            root = ET.fromstring(header + '<root>' + cfg_file.read() + '</root>')
            self.Lmax = int(root[1][0][1].text)+1
            self.Mmax = int(root[1][0][2].text)+1
            print('Parameter file found', self.Lmax, self.Mmax)
        except BaseException as e:
            self.Lmax = self.Mmax = data.shape[0] / 4
            pass

    def draw(self,type='both'):
        pp.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 14
        # set parameters for plotting
        pp.ticklabel_format(style='sci', axis='y')
        try:
            if sys.argv[2]:
                print('Representing the spectrum over time')
                # remove chunks of data
                t = self.datafull[:,0]
                t = t- np.min(t)
                data = self.datafull[:,1:]
                torL = data[:,:self.Lmax]
                polL = data[:,self.Mmax + self.Lmax:self.Mmax + 2 * self.Lmax]
                torM = data[:,self.Lmax:self.Lmax + self.Mmax]
                polM = data[:,self.Mmax+ 2 * self.Lmax:]
                l = np.arange(self.Lmax)
                m = np.arange(self.Mmax)
                fig, (ax) = pp.subplots(2,2)
                ax[0,0].contourf(t, l, np.log(torL.T), 100)
                ax[0,1].contourf(t, l, np.log(polL.T), 100)
                ax[1,0].contourf(t, m, np.log(torM.T), 100)
                ax[1,1].contourf(t, m, np.log(polM.T), 100)
                pp.show()

            else:
                raise ValueError(sys.argv[2], 'doesn\'t specify a correct flag')

        except Exception as e:
            print(e)


            self.draw_snapshot(self.datafull,type=type)

    def draw(self, data, **kwargs):
        # define the routine for plotting datasets
        if kwargs['type'] != 'l':
            pp.loglog(data[:self.Lmax], label='L spectrum, toroidal')
            pp.loglog(data[self.Mmax + self.Lmax:self.Mmax + 2 * self.Lmax],
                      label='L spectrum, poloidal')

        if kwargs['type'] != 'm':
            pp.loglog(np.cumsum(np.ones_like(data[self.Lmax:self.Lmax + self.Mmax])),
                      data[self.Lmax:self.Lmax + self.Mmax], label='M spectrum, toroidal')
            pp.loglog(np.cumsum(np.ones_like(data[self.Mmax + 2 * self.Lmax:2 * self.Mmax + 2 * self.Lmax])),
                      data[self.Mmax + 2 * self.Lmax:2 * self.Mmax + 2 * self.Lmax], label='M spectrum, poloidal')

        pp.xlabel(r'$l\quad m+1$')
        pp.ylabel(r'$E_{kin}$')
        pp.legend(prop={'size': 14})

        BaseRepresenter.draw(self)


    def draw_snapshot(self, datafull, **kwargs):

        # select the last  record of data to be plotted
        data = datafull[-1, 1:]
        #data = datafull[ 1:]

        # check the length of the data vector
        expected_vec_size = (self.Lmax + self.Mmax)*2
        if len(data) == expected_vec_size :
            self.draw(data, kwargs)

        elif len(data) == 2*expected_vec_size:
            self.draw(data[:expected_vec_size], kwargs)
            pp.figure()
            self.draw(data[expected_vec_size:], kwargs)



        """
        index_vector = np.cumsum(np.ones_like(data[self.Lmax:self.Lmax + self.Mmax]))
        pp.loglog(index_vector, index_vector**(-5./3)*3e-4,'--', label=r'$k^{-\frac{5}{3}}$')

        index_vector = index_vector[index_vector > 10.]
        pp.loglog(index_vector, index_vector**(-3.)*3e-2,'--', label=r'$k^{-3}$')

        index_vector = index_vector[index_vector > 40.]
        pp.loglog(index_vector, index_vector**(-5)*3e2, '--', label=r'$k^{-5}$')
        """
        #pp.ylim(ymin=1e-20)

        """
        pp.figure()
        pp.semilogy(data[0:self.Lmax], label='L spectrum, toroidal')
        pp.semilogy(data[self.Lmax:self.Lmax + self.Mmax], label='M spectrum, toroidal')
        pp.semilogy(data[self.Mmax + self.Lmax:self.Mmax+ 2 * self.Lmax], label='L spectrum, poloidal')
        pp.semilogy(data[self.Mmax+ 2 * self.Lmax:], label='M spectrum, poloidal')
        # pp.title(folder_name)
        """

        return pp.gca()

if __name__=="__main__":
    reader = SpectraRepresenter()
    print(reader.name_columns)
    reader.open()
    reader.draw()
