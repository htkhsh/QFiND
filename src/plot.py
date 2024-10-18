import numpy as np
import matplotlib.pyplot as plt

ls1 = 20
ls2 = 15
lw1 = 2.5
lw2 = 2.5
color1 = 'orangered'
color2 = 'royalblue'
color3 = 'black'

def plot_bcf(t, exact, approx, error):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), gridspec_kw={'hspace': 0.04, 'wspace': 0})
    ax1 = axes[0]
    ax2 = axes[1]
    ax1.plot(t, np.real(approx), linestyle='-', label=r'$\mathrm{Re}\;C(t)$', lw=lw1, color=color1)
    ax1.plot(t, np.imag(approx), linestyle='-', label=r'$\mathrm{Im}\;C(t)$', lw=lw1, color=color2)
    ax1.plot(t, np.real(exact), linestyle='--', label='Reference', lw=lw2, color=color3)
    ax1.plot(t, np.imag(exact), linestyle='--', lw=lw2, color=color3)
    ax2.plot(t, np.real(error), linestyle='-', lw=lw1, color=color1)
    ax2.plot(t, np.imag(error), linestyle='-', lw=lw1, color=color2)
    ax1.legend(loc='upper right', fontsize=ls2)

    tmax = t[np.size(t)-1]
    for ax in axes.flat:
        ax.set_xlim(0, tmax)
        ax.tick_params(labelsize=ls2)
    ax1.tick_params(labelbottom=False)
    ax1.set_ylabel('$C(t)$', fontsize=ls1)
    ax2.set_ylabel('$\delta C(t)$', fontsize=ls1)
    ax2.set_xlabel('$t$ $(\mathrm{fs})$', fontsize=ls1)

    plt.subplots_adjust(left=0.17, right=0.95, top=0.99, bottom=0.11, hspace=0, wspace=0)
    plt.savefig('../bcf.png')
    plt.show()
