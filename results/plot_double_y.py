import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import numpy as np


def plot_clustering_metrics(ax):
    labels = ['NMI', 'AMI', 'ARI', 'FMI', 'HOMO', 'COMP']
    values_A = [0.981, 0.947, 0.885, 0.890, 0.997, 0.965]
    values_B = [0.982, 0.952, 0.892, 0.896, 0.997, 0.968]
    values_C = [0.987, 0.966, 0.937, 0.937, 0.994, 0.980]
    values_D = [0.989, 0.971, 0.943, 0.944, 0.994, 0.984]
    x = np.arange(len(labels))

    bars_A = ax.bar(x - 0.3, values_A, width=0.2, label='-Merging-Refinement', color='#376795', alpha=0.85)
    bars_B = ax.bar(x - 0.1, values_B, width=0.2, label='-Merging', color='#72bcd5', alpha=0.85)
    bars_C = ax.bar(x + 0.1, values_C, width=0.2, label='-Refinement', color='#ffd06f', alpha=0.85)
    bars_D = ax.bar(x + 0.3, values_D, width=0.2, label='HSFC', color='#ea6254', alpha=0.85)

    for bars in [bars_A, bars_B, bars_C, bars_D]:
        for bar in bars:
            bar.set_edgecolor('black')
            bar.set_linewidth(0.75)

    ax.set_ylabel('Values')
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.17)

    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    ax.legend(loc='upper right', fancybox=True, ncol=2)


def plot_count_time(ax):
    module = ['-Merging-Refinement', '-Merging', '-Refinement', 'HSFC']
    y1 = [17819, 16048, 15530, 13888]
    y2 = [186.55, 189.50, 192.62, 193.74]

    bars = ax.bar(module, y2, color='#376795', label='Time(s)', width=0.5, alpha=0.8)
    for bar in bars:
        bar.set_edgecolor('black')
        bar.set_linewidth(0.75)
    ax.tick_params(axis='y', labelcolor='#000000')
    ax.set_ylim(95, 205)

    ax_ = ax.twinx()

    ax_.plot(module, y1, color='#ea6254', marker='o', linestyle='-', linewidth=2.5, markersize=10, label='Count')
    ax_.tick_params(axis='y', labelcolor='#000000')
    ax_.set_ylim(9500, 20000)

    ax.tick_params(axis='x', rotation=0)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.legend(loc='upper left', fancybox=True, ncol=4)
    ax_.legend(loc='upper right', fancybox=True, ncol=4)


if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 24
    plt.rc('legend', fontsize=18)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))

    plot_clustering_metrics(ax1)
    plot_count_time(ax2)

    ax1.text(-0.06, 1.06, 'A', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
    ax2.text(-0.06, 1.06, 'B', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')

    plt.tight_layout()

    plt.savefig("ablation", dpi=300)
    # plt.show()
