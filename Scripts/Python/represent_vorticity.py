from base_representer import BaseRepresenter
from matplotlib.ticker import OldScalarFormatter
from matplotlib import pyplot as pp
import numpy as np

class VorticityRepresenter(BaseRepresenter):

    _name_columns = [r'$t$', r'$\omega_x$', r'$\omega_y$', r'$\omega_z$']

    def __init__(self):
        BaseRepresenter.__init__()


    def __init__(self):
        BaseRepresenter.__init__(self)
        pass


    def draw(self):

        data = self.data
        #data[r'$t$'] -= min(data[r'$t$'])
        data = data[data[r'$\omega_z$'] < 1e2]

        ax = data.plot(x=r'$t$', y=[r'$\omega_x$', r'$\omega_y$', r'$\omega_z$'])
        # set parameters for plotting
        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'
        # pp.rcParams['font.size'] = 12

        # ax.set_title(folder_name)#+',  toroidal/total energy ratio: '+ string_ratio+'%')
        ax.set_xlabel('t')
        ax.set_ylabel('omega')

        BaseRepresenter.draw(self)


        try:
            fN = float(self.search_in_parameter('omega'))
            print(fN)
        except AttributeError as e:
            print(e)
            fN=0


        data = self.data
        data['IC_phase'] = (fN*data[r'$t$'])%(2*np.pi)
        data['omega_phase'] = np.arctan2(data[r'$\omega_y$'], data[r'$\omega_x$'])
        data['omega_h'] = (data[r'$\omega_x$']**2+data[r'$\omega_y$']**2)**.5
        data[r'$\theta$'] = np.arctan2(data['omega_h'],data[r'$\omega_z$']+1.)
        data[r'$\delta\phi$'] = (data['omega_phase']-data['IC_phase'])%(2*np.pi)
        #data[r'$t$'] -= min(data[r'$t$'])

        ax = data.plot(x = r'$t$', y = [r'$\delta\phi$', r'$\theta$'])

        ax.yaxis.set_major_formatter(OldScalarFormatter())
        pp.rcParams['font.family'] = 'ubuntu'

        BaseRepresenter.draw(self)
    pass

if __name__ =="__main__":
    reader = VorticityRepresenter()
    reader.open()
    reader.draw()