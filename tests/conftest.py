'''
Module storing data classes needed by tests
Also pytest functions intended to be run after tests
'''
import os

import mplhep
import pandas            as pnd
import matplotlib.pyplot as plt
from rx_efficiencies.decay_names import DecayNames as dn

# -----------------------------------
class ScalesData:
    '''
    data class
    '''
    df_def_wp = pnd.DataFrame(columns=['Process', 'Q2', 'Value', 'Error'])
    df_mva_wp = pnd.DataFrame(columns=['mva_cut', 'Q2', 'Value', 'Error'])

    plt.style.use(mplhep.style.LHCb2)
    # ------------------
    @staticmethod
    def collect_def_wp(proc : str, q2bin : str, value : float, error : float) -> None:
        '''
        Picks test outputs and uses it to fill dataframe
        '''
        size                           = len(ScalesData.df_def_wp)
        ScalesData.df_def_wp.loc[size] = [proc, q2bin, value, error]
    # ------------------
    @staticmethod
    def collect_mva_wp(mvawp: str, q2bin : str, value : float, error : float) -> None:
        '''
        Picks test outputs and uses it to fill dataframe
        '''
        size                          = len(ScalesData.df_mva_wp)
        ScalesData.df_mva_wp.loc[size] = [mvawp, q2bin, value, error]
    # ------------------
    @staticmethod
    def plot_scales_def_wp():
        '''
        Plots scales from dataframe with default WP
        '''
        print(ScalesData.df_def_wp)

        ax = None
        for proc, df_proc in ScalesData.df_def_wp.groupby('Process'):
            decay = dn.tex_from_decay(proc)
            ax = df_proc.plot(x='Q2', y='Value', yerr='Error', label=decay, ax=ax, figsize=(13, 10))

        out_dir = 'plots/prec_scales'
        os.makedirs(out_dir, exist_ok=True)

        plt.ylim(0, 1)
        plt.xlabel('')
        plt.ylabel(r'$N_{PRec}/N_{Signal}$')
        plt.savefig(f'{out_dir}/scales_def_wp.png')
        plt.close()
    # ------------------
    @staticmethod
    def plot_scales_mva_wp():
        '''
        Plots scales from dataframe with scanned MVA WP
        '''
        print(ScalesData.df_mva_wp)

        ax = None
        for q2bin, df_q2bin in ScalesData.df_mva_wp.groupby('Q2'):
            ax = df_q2bin.plot(x='mva_cut', y='Value', yerr='Error', label=q2bin, ax=ax, figsize=(13, 10))

        out_dir = 'plots/prec_scales'
        os.makedirs(out_dir, exist_ok=True)

        plt.ylim(0, 1)
        plt.xlabel('')
        plt.ylabel(r'$N_{PRec}/N_{Signal}$')
        plt.savefig(f'{out_dir}/scales_mva_wp.png')
        plt.close()
# -----------------------------------
def pytest_sessionfinish():
    '''
    Runs at the end
    '''
    ScalesData.plot_scales_def_wp()
    ScalesData.plot_scales_mva_wp()
# -----------------------------------
