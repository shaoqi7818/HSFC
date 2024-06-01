import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy.interpolate import make_interp_spline


def plot_recovered_fraction_20(ax):
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # depth = 20
    fraction_recovered_dbgps = [1.0000, 1.0000, 0.9980, 0.9510, 0.7070, 0.2910, 0.0560, 0.0040, 0.0000, 0.0000]
    fraction_recovered_minhash = [1.0000, 1.0000, 1.0000, 1.0000, 0.9980, 0.9710, 0.8850, 0.6820, 0.4420, 0.2310]
    fraction_recovered_clover = [1.0000, 1.0000, 1.0000, 1.0000, 0.9930, 0.9600, 0.8580, 0.6390, 0.4070, 0.2000]
    fraction_recovered_our = [1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9970, 0.9880, 0.9660, 0.9250]

    error_rates = [x * 100 for x in error_rates]
    fraction_recovered_dbgps = [x * 100 for x in fraction_recovered_dbgps]
    fraction_recovered_minhash = [x * 100 for x in fraction_recovered_minhash]
    fraction_recovered_clover = [x * 100 for x in fraction_recovered_clover]
    fraction_recovered_our = [x * 100 for x in fraction_recovered_our]

    # Create a DataFrame
    data = {
        'Error Rate': error_rates,
        'DBGPS': fraction_recovered_dbgps,
        'Clover': fraction_recovered_clover,
        'Minhash': fraction_recovered_minhash,
        'Our': fraction_recovered_our
    }

    df = pd.DataFrame(data)

    x = df['Error Rate']
    y = df['DBGPS']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='DBGPS', color='#215a59', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Clover']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Clover+Muscle5', color='#ed723f', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Minhash']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Minhash+Muscle5', color='#2F5597', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Our']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='HSFC+Muscle5', color='#c3272b', ls='-', lw=2.5)

    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)

    ax.set_xlabel('Error rates')
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylabel('Recovery rate (×20)')
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylim(-0, 105)
    ax.legend()


def plot_recovered_fraction_30(ax):
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # depth = 30
    fraction_recovered_dbgps = [1.0000, 1.0000, 0.9990, 0.9970, 0.9850, 0.9270, 0.6770, 0.2650, 0.0450, 0.0060]
    fraction_recovered_minhash = [1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9990, 0.9810, 0.9150, 0.7390, 0.5270]
    fraction_recovered_clover = [1.0000, 1.0000, 1.0000, 1.0000, 0.9990, 0.9970, 0.9760, 0.9000, 0.6790, 0.4390]
    fraction_recovered_our = [1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9980, 0.9920]

    error_rates = [x * 100 for x in error_rates]
    fraction_recovered_dbgps = [x * 100 for x in fraction_recovered_dbgps]
    fraction_recovered_minhash = [x * 100 for x in fraction_recovered_minhash]
    fraction_recovered_clover = [x * 100 for x in fraction_recovered_clover]
    fraction_recovered_our = [x * 100 for x in fraction_recovered_our]

    # Create a DataFrame
    data = {
        'Error Rate': error_rates,
        'DBGPS': fraction_recovered_dbgps,
        'Clover': fraction_recovered_clover,
        'Minhash': fraction_recovered_minhash,
        'Our': fraction_recovered_our
    }

    df = pd.DataFrame(data)

    x = df['Error Rate']
    y = df['DBGPS']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='DBGPS', color='#215a59', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Clover']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Clover+Muscle5', color='#ed723f', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Minhash']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Minhash+Muscle5', color='#2F5597', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Our']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='HSFC+Muscle5', color='#c3272b', ls='-', lw=2.5)

    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)

    ax.set_xlabel('Error rates')
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylabel('Recovery rate (×30)')
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylim(-0, 105)
    ax.legend()


def plot_reconstruction_rate_20(ax):
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # depth = 20
    reconstruction_rate_dbgps = [0.9991, 0.9991, 0.9991, 0.9782, 0.8734, 0.6589, 0.4388, 0.3068, 0.2511, 0.2239]
    reconstruction_rate_minhash = [0.9999, 0.9996, 0.9993, 0.9988, 0.9981, 0.9972, 0.9956, 0.9930, 0.9903, 0.9858]
    reconstruction_rate_clover = [1.0000, 0.9998, 0.9994, 0.9988, 0.9977, 0.9954, 0.9929, 0.9884, 0.9846, 0.9786]
    reconstruction_rate_our = [1.0000, 1.0000, 0.9999, 0.9997, 0.9996, 0.9995, 0.9991, 0.9984, 0.9975, 0.9957]

    error_rates = [x * 100 for x in error_rates]
    reconstruction_rate_dbgps = [x * 100 for x in reconstruction_rate_dbgps]
    reconstruction_rate_minhash = [x * 100 for x in reconstruction_rate_minhash]
    reconstruction_rate_clover = [x * 100 for x in reconstruction_rate_clover]
    reconstruction_rate_our = [x * 100 for x in reconstruction_rate_our]

    # Create a DataFrame
    data = {
        'Error Rate': error_rates,
        'DBGPS': reconstruction_rate_dbgps,
        'Clover': reconstruction_rate_clover,
        'Minhash': reconstruction_rate_minhash,
        'Our': reconstruction_rate_our
    }

    df = pd.DataFrame(data)

    x = df['Error Rate']
    y = df['DBGPS']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='DBGPS', color='#215a59', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Clover']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Clover+Muscle5', color='#ed723f', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Minhash']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Minhash+Muscle5', color='#2F5597', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Our']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='HSFC+Muscle5', color='#c3272b', ls='-', lw=2.5)

    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)

    ax.set_xlabel('Error rates')
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylabel('Reconstruction rate (×20)')
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylim(-0, 105)
    ax.legend()


def plot_reconstruction_rate_30(ax):
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # depth = 30
    reconstruction_rate_dbgps = [0.9991, 0.9991, 0.9991, 0.9989, 0.9964, 0.9692, 0.8689, 0.6439, 0.4380, 0.3119]
    reconstruction_rate_minhash = [0.9997, 0.9996, 0.9993, 0.9990, 0.9986, 0.9976, 0.9961, 0.9944, 0.9919, 0.9887]
    reconstruction_rate_clover = [1.0000, 0.9998, 0.9996, 0.9992, 0.9985, 0.9972, 0.9949, 0.9916, 0.9871, 0.9816]
    reconstruction_rate_our = [1.0000, 1.0000, 0.9998, 0.9997, 0.9996, 0.9994, 0.9991, 0.9984, 0.9974, 0.9964]

    error_rates = [x * 100 for x in error_rates]
    reconstruction_rate_dbgps = [x * 100 for x in reconstruction_rate_dbgps]
    reconstruction_rate_minhash = [x * 100 for x in reconstruction_rate_minhash]
    reconstruction_rate_clover = [x * 100 for x in reconstruction_rate_clover]
    reconstruction_rate_our = [x * 100 for x in reconstruction_rate_our]

    # Create a DataFrame
    data = {
        'Error Rate': error_rates,
        'DBGPS': reconstruction_rate_dbgps,
        'Clover': reconstruction_rate_clover,
        'Minhash': reconstruction_rate_minhash,
        'Our': reconstruction_rate_our
    }

    df = pd.DataFrame(data)

    x = df['Error Rate']
    y = df['DBGPS']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='DBGPS', color='#215a59', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Clover']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Clover+Muscle5', color='#ed723f', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Minhash']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='Minhash+Muscle5', color='#2F5597', ls='-', lw=2.5)

    x = df['Error Rate']
    y = df['Our']
    m = make_interp_spline(x, y)
    xs = np.linspace(x.min(), x.max(), 1000)
    ys = m(xs)
    ax.plot(xs, ys, label='HSFC+Muscle5', color='#c3272b', ls='-', lw=2.5)

    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)

    ax.set_xlabel('Error rates')
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylabel('Reconstruction rate (×30)')
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_ylim(-0, 105)
    ax.legend()


if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 20
    plt.rc('legend', fontsize=16)

    fig, axs = plt.subplots(2, 2, figsize=(19, 13))

    plot_recovered_fraction_30(axs[0, 0])
    plot_reconstruction_rate_30(axs[0, 1])
    plot_recovered_fraction_20(axs[1, 0])
    plot_reconstruction_rate_20(axs[1, 1])

    axs[0, 0].text(-0.105, 1.09, 'A', transform=axs[0, 0].transAxes, fontsize=22, fontweight='bold', va='top',
                   ha='right')
    axs[0, 1].text(-0.105, 1.09, 'B', transform=axs[0, 1].transAxes, fontsize=22, fontweight='bold', va='top',
                   ha='right')
    axs[1, 0].text(-0.105, 1.09, 'C', transform=axs[1, 0].transAxes, fontsize=22, fontweight='bold', va='top',
                   ha='right')
    axs[1, 1].text(-0.105, 1.09, 'D', transform=axs[1, 1].transAxes, fontsize=22, fontweight='bold', va='top',
                   ha='right')

    plt.tight_layout()
    plt.savefig("recovery and reconstruction.png", dpi=300)
    # plt.show()
