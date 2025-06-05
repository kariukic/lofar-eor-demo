# From soltool plot by Florent

import os
import click
import numpy as np
import losoto.h5parm
import matplotlib.pyplot as plt
from multiprocessing import Pool


class GainSol(object):

    def __init__(self, time, freqs, ant, directions, pol, amp, phase):
        self.time = time
        self.freqs = freqs
        self.ant = ant
        self.dir = directions
        self.pol = pol
        self.amp = amp
        self.phase = phase
        self.d = self.amp * np.exp(1j * self.phase)


def open_sol(file_h5):
    # amp: time, freqs, antm, dir, pol
    sol_file = None
    try:
        sol_file = losoto.h5parm.h5parm(file_h5)
        solset = sol_file.getSolsets()[0]
        soltab, soltab_phase = solset.getSoltabs(useCache=True)

        ant = soltab.getAxisValues('ant')
        directions = soltab.getAxisValues('dir')
        time = soltab.getAxisValues('time')
        pol = soltab.getAxisValues('pol')

        freqs = soltab.getAxisValues('freq')

        weight = soltab.getValues(weight=True)[0].astype(bool)
        amp = np.ma.array(soltab.getValues(weight=False)[0], mask=~weight)
        phase = np.ma.array(soltab_phase.getValues(weight=False)[0], mask=~weight)

        if directions is None:
            amp = amp[:, :, :, None, :]
            phase = phase[:, :, :, None, :]
            directions = ['di']
    finally:
        if sol_file is not None:
            sol_file.close()

    return GainSol(time, freqs, ant, directions, pol, amp, phase)


def plot_sol(sol, dir, pol, data_type, filename):
    if data_type == 'Amplitude':
        v = sol.amp[:, :, :, dir, pol]
    elif data_type == 'Phase':
        v = sol.phase[:, :, :, dir, pol]
    else:
        print(f'Error: data type {data_type} unknown')
        return

    vmax = np.nanquantile(v[~v.mask & ~np.isnan(v) & (v != 0)], 0.999)
    vmin = np.nanquantile(v[~v.mask & ~np.isnan(v) & (v != 0)], 0.001)
    extent = [0, len(sol.time), sol.freqs.min() * 1e-6, sol.freqs.max() * 1e-6]

    n = v.shape[2]
    ncols, nrows = int(np.ceil(np.sqrt(n))), int(np.ceil(n / np.ceil(np.sqrt(n))))

    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, sharey=True, figsize=(1 + 2 * ncols, 1 + 1.5 * nrows),
                            sharex=True)

    im = None

    for i, ax in zip(range(v.shape[2]), axs.flatten()):
        if v.shape[0] > 1 and v.shape[1] > 1:
            im = ax.imshow(v[:, :, i].T, aspect='auto', vmax=vmax, vmin=vmin, extent=extent)
        elif v.shape[0] == 1:
            ax.plot(sol.freqs * 1e-6, v[0, :, i].T)
        elif v.shape[1] == 1:
            ax.plot(v[:, 0, i].T)
        ax.text(0.025, 0.975, sol.ant[i], transform=ax.transAxes, fontsize=11, va='top')

    ylabel = ''
    xlabel = ''

    if v.shape[0] > 1 and v.shape[1] > 1 and im is not None:
        cax = fig.add_axes([0.6, 1.04, 0.39, 0.02])
        cax.set_xlabel(data_type)
        fig.colorbar(im, cax=cax, orientation='horizontal')
        xlabel = 'Time index'
        ylabel = 'Frequency [Mhz]'
    elif v.shape[0] == 1:
        xlabel = 'Frequency [MHz]'
        ylabel = data_type
    elif v.shape[1] == 1:
        xlabel = 'Time index'
        ylabel = data_type

    for ax in axs[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axs[-1, :]:
        ax.set_xlabel(xlabel)

    fig.tight_layout(pad=0)
    print(filename)
    fig.savefig(filename, dpi=120, bbox_inches="tight")


@click.group()
def main():
    ''' plot DP3 solutions ...'''


@main.command('plot')
@click.argument('sol_file')
@click.option('--n_cpu', help='Number of CPU to use', type=int, default=4)
def plot(sol_file, n_cpu):
    ''' Plot solutions of the h5 file SOLS '''
    sol = open_sol(sol_file)
    with Pool(n_cpu) as pool:
        for data_type in ['Amplitude', 'Phase']:
            for dir in range(len(sol.dir)):
                for pol in range(len(sol.pol)):
                    dirname = sol.dir[dir].strip('[').strip(']')
                    filename = f'{data_type}_dir{dirname}_pol{sol.pol[pol]}.png'
                    pool.apply_async(plot_sol, [sol, dir, pol, data_type, filename])
        pool.close()
        pool.join()


if __name__ == '__main__':
    main()
