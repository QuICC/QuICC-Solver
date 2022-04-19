# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
from base_representer import BaseRepresenter
from matplotlib.ticker import OldScalarFormatter
from matplotlib import pyplot as pp

class EnergyRepresenter(BaseRepresenter):

    _name_columns = [r'$t$', r'$E_{tot}$', r'$E_{tor}$', r'$E_{pol}$', r'$E^c_{sym}$', r'$E^c_{asym}$', r'$E^{eq}_{sym}$', r'$E^{eq}_{asym}$']

    def __init__(self):
        BaseRepresenter.__init__(self)
        pass

    def draw(self):
        data = self.data
        data[r'$t$'] -= min(data[r'$t$'])
        data = data[data[r'$E_{tot}$'] < 1e2]

        idx = max(data.index)
        string_ratio = str("%.2f" % (data[r'$E_{tor}$'][idx] / data[r'$E_{tot}$'][idx] * 100))
        # print('Final Toroidal to Total energy ratio of: '+string_ratio)

        #ax = data.plot(x=r'$t$', y=[r'$E_{tot}$', r'$E_{tor}$', r'$E_{pol}$', r'$E^c_{sym}$', r'$E^c_{asym}$', r'$E^{eq}_{sym}$', r'$E^{eq}_{asym}$'], logy = True)
        ax = data.plot(x=r'$t$', y=[r'$E_{tot}$', r'$E_{tor}$', r'$E_{pol}$', r'$E^c_{sym}$', r'$E^c_{asym}$'], logy = True)
        #ax = data.plot(x=r'$t$', y=[r'$E_{tot}$', r'$E_{tor}$', r'$E_{pol}$'], logy=True)
        # set parameters for plotting
        #ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 12

        # ax.set_title(folder_name)#+',  toroidal/total energy ratio: '+ string_ratio+'%')
        ax.set_xlabel('t')
        ax.set_ylabel(r'$E_{kin}$')

        BaseRepresenter.draw(self)

        return ax


if __name__=="__main__":

    reader = EnergyRepresenter()
    reader.open()
    reader.draw()







    
    
    
    
