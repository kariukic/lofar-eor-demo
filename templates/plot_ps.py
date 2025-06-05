import click
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pspipe import settings, database

mpl.style.use("default")
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["image.interpolation"] = "none"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.aspect"] = "auto"


class Spectra:
    def __init__(
        self,
        toml_file: str,
        obsid: str,
        do_flag: bool = False,
        fmhz_range: list = [],
        eor_bin: int = 0,
        flags_to_use=None
    ):
        self.toml_file = toml_file
        self.obsid = obsid
        self.flag = do_flag
        self.fmhz_range = fmhz_range
        self.eor_bin = eor_bin
        self.flags_to_use = flags_to_use

        self.revision = database.VisRevision(
            settings.Settings.load_with_defaults(self.toml_file)
        )

        self.data = self.revision.get_data(self.obsid)
        self.data.filter_uvrange(50, 250)

        if self.flag:
            self.data.do_flag()

        if self.flags_to_use:
            self.data.set_flag(self.flags_to_use)

        self.ps_gen = self.data.get_ps_gen(
            window_fct='blackmanharris',
            rmean_freqs=False,
            umin=50,
            umax=250,
            du=8,
            empirical_weighting=True,
            empirical_weighting_polyfit_deg='2',
            primary_beam='lofar_hba',
            kbins_n=8,
            ft_method='lssa',
        )

        self.z_val = self.ps_gen.eor.z

        self.kbins = np.logspace(np.log10(self.ps_gen.kmin), np.log10(0.5), 10)

    def plot(self, plotdir):
        fig = plt.figure(figsize=(8, 8), dpi=300)
        gs = gridspec.GridSpec(3, 3, figure=fig)

        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])

        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])
        ax5 = fig.add_subplot(gs[1, 2])

        ax6 = fig.add_subplot(gs[2, 0])
        ax7 = fig.add_subplot(gs[2, 1:3])

        self.ps_gen.get_ps2d(self.data.i).plot(ax=ax0, title='I')
        self.ps_gen.get_ps2d(self.data.v).plot(ax=ax1, title='V')
        (self.ps_gen.get_ps2d(self.data.i) / self.ps_gen.get_ps2d(self.data.v)).plot(ax=ax2, title='I/V')

        self.ps_gen.get_ps(self.data.i).plot(ax=ax3, title='I')
        self.ps_gen.get_ps(self.data.v).plot(ax=ax4, title='V')
        (self.ps_gen.get_ps(self.data.i) / self.ps_gen.get_ps(self.data.v)).plot(ax=ax5, title='I/V')

        self.ps_gen.get_variance(self.data.i).plot(ax=ax6, label='I')
        self.ps_gen.get_variance(self.data.v).plot(ax=ax6, label='V')
        self.ps_gen.get_variance(self.data.v_dt).plot(ax=ax6, label='V_dt')
        ax6.legend()

        self.ps_gen.get_ps3d(self.kbins, self.data.i).plot(ax=ax7, label='I')
        self.ps_gen.get_ps3d(self.kbins, self.data.v).plot(ax=ax7, label='V')
        self.ps_gen.get_ps3d(self.kbins, self.data.v_dt).plot(ax=ax7, label='V_dt')
        ax7.legend()

        for a in [ax0, ax1, ax2]:
            a.set_ylim(0, 1.2)

        fig.tight_layout()
        plt.savefig(f"{plotdir}/ps_spectra.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        return


@click.group()
def main():
    ''' plot flags occupancy ...'''


@main.command('plot_ps')
@click.argument('config')
@click.option('--obsid', help='pspipe obsID', type=str, default='')
@click.option('--plotdir', help='plots directory', type=str, default='')
def plot_aoqstats(config, obsid, plotdir=None):
    ''' Plot aoquality statistics '''
    spec = Spectra(toml_file=config, obsid=obsid)
    spec.plot(plotdir)


if __name__ == '__main__':
    main()
