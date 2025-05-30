import numpy as np
from brokenaxes import brokenaxes
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D

import matplotlib
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['font.sans-serif'] = ['Times New Roman']

def plot_time_vs_alpha_saw_temperature(df):
   
    fig, axs = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    axs[0].plot(df['time[min]'], df['Temperature'], color='black', label='Temperature')
    axs[0].set_ylabel('Temperature [K]', fontsize=18, labelpad=10)
    axs[0].grid(True, which='both', linestyle='--', linewidth=0.75)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].set_xlim([0, 240])
    axs[0].tick_params(axis='x', labelsize=16)
    axs[0].tick_params(axis='y', labelsize=16)

    axs[1].errorbar(df['time[min]'], df['alpha[m^2s^-1]'] * 1e5, 
                    yerr=df['alpha_err[m^2s^-1]'] * 1e5, 
                    fmt='o', linestyle='none', capsize=5, color='red', label='Thermal Diffusivity')
    axs[1].set_ylabel('Thermal Diffusivity [$10^{-5}$ m$^2$s$^{-1}$]', fontsize=18, labelpad=10)
    axs[1].grid(True, which='both', linestyle='--', linewidth=0.75)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].set_xlim([0, 240])
    axs[1].tick_params(axis='x', labelsize=16)
    axs[1].tick_params(axis='y', labelsize=16)

    axs[2].errorbar(df['time[min]'], df['f[Hz]'] * df['grating_spacing[µm]'] * 1e-6,
                    yerr=df['f_err[Hz]'] * df['grating_spacing[µm]'] * 1e-6,
                    fmt='o', linestyle='none', capsize=5, color='blue', label='SAW Speed')
    axs[2].set_ylabel('SAW Speed [ms$^{-1}$]', fontsize=18, labelpad=10)
    axs[2].set_xlabel('Time [min]', fontsize=18, labelpad=10)
    axs[2].grid(True, which='both', linestyle='--', linewidth=0.75)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    axs[2].set_xlim([0, 240])
    axs[2].tick_params(axis='x', labelsize=16)
    axs[2].tick_params(axis='y', labelsize=16)

    plt.tight_layout()
    plt.savefig('figures/time_vs_alpha_saw_temperature.png', dpi=600)

def plot_fluence_vs_alpha(df):

    plt.figure(figsize=(12, 8))
    plt.errorbar(df['fluence[ionsm^-2]'] * 1e-17, 
                 df['alpha[m^2s^-1]'] * 1e4, 
                 yerr=df['alpha_err[m^2s^-1]'] * 1e4, 
                 fmt='o', linestyle='none', capsize=5, color='#B23A4A',
                 label='TGS Measured Data', zorder=3)
    
    x = df['fluence[ionsm^-2]'] * 1e-17
    y = df['alpha[m^2s^-1]'] * 1e4
    errors = df['alpha_err[m^2s^-1]'] * 1e4
    
    sort_idx = np.argsort(x)
    x_sorted = x.iloc[sort_idx]
    y_sorted = y.iloc[sort_idx]
    errors_sorted = errors.iloc[sort_idx]
    
    window_size = min(6, len(y_sorted))
    y_ma = y_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    errors_ma = errors_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    
    plt.fill_between(x_sorted, 
                     y_ma - errors_ma,
                     y_ma + errors_ma,
                     color='red', alpha=0.2,
                     label='1σ Confidence Interval',
                     zorder=1)
    
    plt.plot(x_sorted, y_ma, color='red', linewidth=2, 
             label='Moving Average', zorder=2)
    
    plt.xlabel('Fluence [10$^{17}$ ions m$^{-2}$]', fontsize=18, labelpad=10)
    plt.ylabel('Thermal Diffusivity [10$^{-4}$ m$^2$s$^{-1}$]', fontsize=18, labelpad=10)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, which='both', linestyle='--', linewidth=0.75)
    plt.legend(fontsize=16, handles=plt.gca().get_legend_handles_labels()[0][::-1])
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('figures/fluence_vs_alpha.png', dpi=600)

def plot_dpa_vs_alpha(df):

    plt.figure(figsize=(12, 8))
    plt.errorbar(df['dpa'], df['alpha[m^2s^-1]'] * 1e4, 
                yerr=df['alpha_err[m^2s^-1]'] * 1e4, 
                fmt='o', linestyle='none', capsize=5, color='#B23A4A',
                label='TGS Measured Data', zorder=3)
    
    x = df['dpa']
    y = df['alpha[m^2s^-1]'] * 1e4
    errors = df['alpha_err[m^2s^-1]'] * 1e4
    
    sort_idx = np.argsort(x)
    x_sorted = x.iloc[sort_idx]
    y_sorted = y.iloc[sort_idx]
    errors_sorted = errors.iloc[sort_idx]
    
    window_size = min(6, len(y_sorted))
    y_ma = y_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    errors_ma = errors_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    
    plt.fill_between(x_sorted, 
                     y_ma - errors_ma,
                     y_ma + errors_ma,
                     color='red', alpha=0.2,
                     label='1σ Confidence Interval',
                     zorder=1)
    
    plt.plot(x_sorted, y_ma, color='red', linewidth=2, 
             label='Moving Average', zorder=2)
    
    plt.xlabel('Dpa', fontsize=18, labelpad=10)
    plt.ylabel('Thermal Diffusivity [10$^{-4}$ m$^2$s$^{-1}$]', fontsize=18, labelpad=10)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, which='both', linestyle='--', linewidth=0.75)
    plt.legend(fontsize=16, handles=plt.gca().get_legend_handles_labels()[0][::-1])
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('figures/dpa_vs_alpha.png', dpi=600)

def plot_dpa_vs_saw(df):

    plt.figure(figsize=(12, 8))
    plt.errorbar(df['dpa'], df['f[Hz]'] * df['grating_spacing[µm]'] * 1e-6, 
                yerr=df['f_err[Hz]'] * df['grating_spacing[µm]'] * 1e-6, 
                fmt='o', linestyle='none', capsize=5, color='#3A4AB2',
                label='TGS Measured Data', zorder=3)
    
    x = df['dpa']
    y = df['f[Hz]'] * df['grating_spacing[µm]'] * 1e-6
    errors = df['f_err[Hz]'] * df['grating_spacing[µm]'] * 1e-6
    
    sort_idx = np.argsort(x)
    x_sorted = x.iloc[sort_idx]
    y_sorted = y.iloc[sort_idx]
    errors_sorted = errors.iloc[sort_idx]
    
    window_size = min(6, len(y_sorted))
    y_ma = y_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    errors_ma = errors_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    
    plt.fill_between(x_sorted, 
                     y_ma - errors_ma,
                     y_ma + errors_ma,
                     color='blue', alpha=0.2,
                     label='1σ Confidence Interval',
                     zorder=1)
    
    plt.plot(x_sorted, y_ma, color='blue', linewidth=2, 
             label='Moving Average', zorder=2)
    
    plt.xlabel('Dpa', fontsize=18, labelpad=10)
    plt.ylabel('SAW Speed [ms$^{-1}$]', fontsize=18, labelpad=10)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, which='both', linestyle='--', linewidth=0.75)
    plt.legend(fontsize=16, handles=plt.gca().get_legend_handles_labels()[0][::-1])
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('figures/dpa_vs_saw.png', dpi=600)

def plot_fluence_vs_saw(df):

    plt.figure(figsize=(12, 8))
    plt.errorbar(df['fluence[ionsm^-2]'] * 1e-17, 
                df['f[Hz]'] * df['grating_spacing[µm]'] * 1e-6, 
                yerr=df['f_err[Hz]'] * df['grating_spacing[µm]'] * 1e-6, 
                fmt='o', linestyle='none', capsize=5, color='#3A4AB2',
                label='TGS Measured Data', zorder=3)
    
    x = df['fluence[ionsm^-2]'] * 1e-17
    y = df['f[Hz]'] * df['grating_spacing[µm]'] * 1e-6
    errors = df['f_err[Hz]'] * df['grating_spacing[µm]'] * 1e-6
    
    sort_idx = np.argsort(x)
    x_sorted = x.iloc[sort_idx]
    y_sorted = y.iloc[sort_idx]
    errors_sorted = errors.iloc[sort_idx]
    
    window_size = min(6, len(y_sorted))
    y_ma = y_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    errors_ma = errors_sorted.rolling(window=window_size, min_periods=1, center=True).mean()
    
    plt.fill_between(x_sorted, 
                     y_ma - errors_ma,
                     y_ma + errors_ma,
                     color='blue', alpha=0.2,
                     label='1σ Confidence Interval',
                     zorder=1)
    
    plt.plot(x_sorted, y_ma, color='blue', linewidth=2, 
             label='Moving Average', zorder=2)
    
    plt.xlabel('Fluence [10$^{17}$ ions m$^{-2}$]', fontsize=18, labelpad=10)
    plt.ylabel('SAW Speed [ms$^{-1}$]', fontsize=18, labelpad=10)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, which='both', linestyle='--', linewidth=0.75)
    plt.legend(fontsize=16, handles=plt.gca().get_legend_handles_labels()[0][::-1])
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('figures/fluence_vs_saw.png', dpi=600)

def plot_time_vs_temperature(ref_df, tgs_df, ion_df, ref_tgs_time, ion_time):

    fig = plt.figure(figsize=(10, 6))
    bax = brokenaxes(xlims=((0, ref_tgs_time[-1]), (ion_time[0], ion_time[-1])), hspace=.05, wspace=.03)

    bax.plot(ref_tgs_time, ref_df['Temperature'], color='blue', linewidth=3, label='REF')
    bax.plot(ref_tgs_time, tgs_df['Temperature'], color='red', linewidth=3, label='TGS')
    bax.plot(ion_time, ion_df['Temperature'], color='orange', linewidth=3, label='TGS+ION')

    bax.set_xlabel('Time [h]', fontsize=18, labelpad=30)
    bax.set_ylabel('Temperature [K]', fontsize=18, labelpad=40)
    bax.tick_params(axis='both', which='major', labelsize=16)
    bax.legend(loc='upper right', fontsize=18)
    bax.grid(True, which='both', linestyle='--', linewidth=0.75)

    for ax in bax.axs:
        ax.yaxis.set_major_locator(MultipleLocator(30))

    plt.savefig('figures/time_vs_temperature.png', dpi=600)

def plot_time_vs_temperature_zoom(ref_df, tgs_df, ion_df, ref_tgs_time, ion_time, zoom_duration):

    ref_df_zoom = ref_df[-zoom_duration:]
    tgs_df_zoom = tgs_df[-zoom_duration:]
    ion_df_zoom = ion_df

    fig = plt.figure(figsize=(10, 6))
    bax = brokenaxes(xlims=((ref_tgs_time[-zoom_duration], ref_tgs_time[-1]), (ion_time[0], ion_time[-1])), hspace=.05, wspace=.05)

    bax.plot(ref_tgs_time[-zoom_duration:], ref_df_zoom['Temperature'], color='blue', linewidth=2, label='REF')
    bax.plot(ref_tgs_time[-zoom_duration:], tgs_df_zoom['Temperature'], color='red', linewidth=2, label='TGS')
    bax.plot(ion_time, ion_df_zoom['Temperature'], color='orange', linewidth=2, label='TGS+ION')

    bax.set_xlabel('Time [h]', fontsize=18, labelpad=30)
    bax.set_ylabel('Temperature [K]', fontsize=18, labelpad=40)
    bax.tick_params(axis='both', which='major', labelsize=16)
    bax.legend(loc='upper right', fontsize=18)
    bax.grid(True, which='both', linestyle='--', linewidth=0.75)

    for ax in bax.axs:
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))

    plt.savefig('figures/time_vs_temperature_zoom.png', dpi=600)

def plot_depth_vs_dpa(df):
    
    fig = plt.figure(figsize=(10, 6))
    gs = plt.GridSpec(1, 1)
    
    ax_main = plt.subplot(gs[0, :])
    ax_main.fill_between(df['depth[um]'], 
                       df['dpa'] * 1e3, 
                       color='purple', 
                       alpha=0.2)
    probe_depth = 1.12
    pattern = '//'
    ax_main.fill_between([0, probe_depth], 0, 40,
                        color='red', alpha=0.1,
                        hatch=pattern)
    
    ax_main.plot(df['depth[um]'], df['dpa'] * 1e3, 
                color='purple', linewidth=2)
    ax_main.axvline(x=probe_depth, color='red', 
                    linestyle='--', linewidth=2)

    ax_main.text(probe_depth/2, 29, 'TGS Thermal\nProbe Region',
                ha='center', va='top', fontsize=18,
                color='red')
    
    ax_main.tick_params(axis='both', which='major', labelsize='large')
    
    ax_main.set_xlabel('Depth [μm]', fontsize=18)
    ax_main.set_ylabel('Displacements per Atom [mdpa]', fontsize=18)
    ax_main.tick_params(axis='both', which='major', labelsize=16)
    ax_main.set_ylim(0, 35)
    ax_main.set_xlim(0.0, 3)
    ax_main.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/depth_vs_dpa.png', dpi=600)
