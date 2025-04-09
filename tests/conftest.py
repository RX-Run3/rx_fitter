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

        plt.figure(figsize=(20,20))
        for proc, df_proc in ScalesData.df_def_wp.groupby('Process'):
            decay = dn.tex_from_decay(proc)
            plt.plot(df_proc.Q2, df_proc.Value, label=decay, color='blue')
            plt.fill_between(df_proc.Q2, df_proc.Value - df_proc.Error, df_proc.Value + df_proc.Error, color='blue', alpha=0.2)

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
        df_mva            = ScalesData.df_mva_wp
        df_mva['mva_cut'] = df_mva.mva_cut.str.replace('mva_', '')
        df_mva['mva_cut'] = df_mva.mva_cut.str.replace('&&'  , '')
        df_mva['mva_cut'] = df_mva.mva_cut.str.replace('('   , '')
        df_mva['mva_cut'] = df_mva.mva_cut.str.replace(')'   , '')

        plt.figure(figsize=(30,20))
        for q2bin, df in df_mva.groupby('Q2'):
            plt.plot(df.mva_cut, df.Value, label=q2bin, color='blue')
            plt.fill_between(df.mva_cut, df.Value - df.Error, df.Value + df.Error, color='blue', alpha=0.2)

        out_dir = 'plots/prec_scales'
        os.makedirs(out_dir, exist_ok=True)

        plt.legend()
        plt.grid()
        plt.ylim(0.0, 0.4)
        plt.xlabel('')
        plt.xticks(rotation=70)
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
