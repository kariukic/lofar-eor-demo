import click
import aoquality as ao
import matplotlib.pyplot as plt


class QualityStats(object):

    def __init__(self, qs_file, name):
        self.qs_file = qs_file
        self.name = name
        self.available_stats = ['RFIPercentage', 'Mean', 'Std', 'DStd', 'Variance', 'DVariance', 'SNR']
        self.pols = ['XX', 'XY', 'YX', 'YY']

    def plot_freq_qstats(self):
        aof = ao.AOQualityFrequencyStat(self.qs_file)

        fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 8), sharex=True)
        axes = axes.ravel()
        for k, stat in enumerate(self.available_stats):
            for p in [0, 3]:
                axes[k].plot(aof.freqs/1e6, aof.get_stat(stat)[:, p])
            axes[k].set_title(stat)
            axes[k].set_yscale('log')
        fig.tight_layout()
        plt.savefig("aoq_frequency_stats.pdf", bbox_inches='tight')
        plt.close(fig)
        return

    def plot_time_qstats(self):
        aot = ao.AOQualityTimeStat(self.qs_file)

        fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 8), sharex=True)
        axes = axes.ravel()
        x_ax = range(aot.Std.shape[0])
        for k, stat in enumerate(self.available_stats):
            for p in [0, 3]:
                axes[k].plot(x_ax, aot.get_stat(stat)[:, p])
            axes[k].set_ylabel(stat)
            axes[k].set_yscale('log')
        fig.tight_layout()
        plt.savefig("aoq_time_stats.pdf", bbox_inches='tight')
        plt.close(fig)
        return

    def plot_baseline_qstats(self):
        aob = ao.AOQualityBaselineStat(self.qs_file)

        for k, stat in enumerate(self.available_stats):
            for p in [0, 3]:
                fb1 = aob.plot_baseline_stats(stat, log=False, name='', pol=p)
                fb2 = aob.plot_antennae_stats(stat, log=False, name='', pol=p)
                fb3 = aob.plot_baseline_length_stats(stat, log=False, name='', pol=p)

                fb1.savefig(f"aoq_{self.name}_baseline_{stat}_{self.pols[p]}.pdf", bbox_inches='tight')
                fb2.savefig(f"aoq_{self.name}_antenna_{stat}_{self.pols[p]}.pdf", bbox_inches='tight')
                fb3.savefig(f"aoq_{self.name}_baseline_length_{stat}_{self.pols[p]}.pdf", bbox_inches='tight')
                for ff in [fb1, fb2, fb3]:
                    plt.close(ff)
        return


@click.group()
def main():
    ''' plot flags occupancy ...'''


@main.command('plot_aoq')
@click.argument('qstats_file')
@click.option('--name', help='output name', type=str, default='')
def plot_aoqstats(qstats_file, name=''):
    ''' Plot aoquality statistics '''

    qq = QualityStats(qstats_file, name)

    qq.plot_time_qstats()
    qq.plot_freq_qstats()
    qq.plot_baseline_qstats()


if __name__ == '__main__':
    main()
