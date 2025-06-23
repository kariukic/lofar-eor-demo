import click
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
import numpy as np
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib import gridspec
from glob import glob


class Image(object):
    def __init__(self, fitspath, pix_box=[100, 100]):
        self.fitspath = fitspath
        self.pix_box = pix_box

        with fits.open(self.fitspath) as hdus:
            img_hdu = hdus["PRIMARY"]

            self.data_array = img_hdu.data.squeeze()
            self.header = img_hdu.header
            try:
                self.image_ID = img_hdu.header["OBJECT"]
            except KeyError:
                image_id = int("".join(filter(str.isdigit, fitspath)))
                self.image_ID = image_id
            self.obsdate = img_hdu.header["DATE-OBS"]
            self.image_size = [img_hdu.header["NAXIS1"], img_hdu.header["NAXIS2"]]
            self.xcellsize = np.abs(img_hdu.header["CDELT1"])
            self.ycellsize = np.abs(img_hdu.header["CDELT2"])
            self.beam_major = img_hdu.header["BMAJ"]
            self.beam_minor = img_hdu.header["BMIN"]
            self.beam_parallactic_angle = img_hdu.header["BPA"]
            self.beam_major_px = self.beam_major / self.xcellsize
            self.beam_minor_px = self.beam_minor / self.ycellsize
            self.beam_area = self.beam_major * self.beam_minor * np.pi
            self.beam_npix = self.beam_area / (self.xcellsize * self.ycellsize)
            self.beam_radius_px = np.sqrt(self.beam_major_px**2 + self.beam_minor_px**2)
            self.mean = np.nanmean(self.data_array)
            self.rms = np.sqrt(np.nanmean(self.data_array**2))
            self.std = np.nanstd(self.data_array)
            self.polarization = img_hdu.header["CRVAL4"]
            region = self.data_array[0 : self.pix_box[0], 0 : self.pix_box[1]]
            self.mean_across_box = np.nanmean(region)
            self.std_across_box = np.nanstd(region)
            self.rms_across_box = np.sqrt(np.nanmean(region**2))

    def plot_image(
        self,
        ax=None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        title=None,
        cbar=True,
        cbar_label="mJy/PSF",
        cbar_ax=None,
        **kwargs,
    ):
        """
        Plot the image with WCS coordinates.
        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            If provided, will plot on this axes. Otherwise creates a new figure.
        cmap : str, optional
            Colormap to use (default: 'viridis')
        vmin, vmax : float, optional
            Minimum and maximum values for the colormap
        title : str, optional
            Title for the plot
        cbar : bool, optional
            Whether to add a colorbar (default: True)
        cbar_label : str, optional
            Label for the colorbar (default: "mJy/PSF")
        cbar_ax : matplotlib.axes.Axes, optional
            Axes to draw the colorbar in (for multi-panel plots)
        kwargs : dict
            Additional arguments passed to imshow
        """
        my_wcs = wcs.WCS(self.header, naxis=[wcs.WCSSUB_CELESTIAL])
        # Create new figure if no axes provided
        if ax is None:
            fig = plt.figure(figsize=(5, 4))
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=my_wcs)
            standalone = True
        else:
            # Ensure the existing axes has the correct WCS projection
            if not isinstance(ax, WCSAxes):
                ax.set_projection(my_wcs)
            standalone = False
        # Set data limits
        if vmin is None:
            vmin = np.nanmin(self.data_array)
        if vmax is None:
            vmax = np.nanmax(self.data_array)
        # Plot the data
        plot_data = self.data_array.squeeze()
        im = ax.imshow(
            plot_data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, **kwargs
        )
        # Formatting
        ax.set_xlabel("R.A. [deg]")
        ax.set_ylabel("Dec [deg]")
        ax.grid(ls="dotted", color="black")
        if title:
            ax.set_title(title)
        # Handle colorbar
        if cbar:
            if standalone or cbar_ax is None:
                # For standalone plot or no cbar_ax specified
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label(cbar_label)
                # Format colorbar ticks
                func = lambda x, pos: "{:g}".format(x * 1e3)
                fmt = mpl.ticker.FuncFormatter(func)
                cbar.ax.yaxis.set_major_formatter(fmt)
            else:
                # For multi-panel plot with shared colorbar
                cbar = plt.colorbar(im, cax=cbar_ax)
                cbar.set_label(cbar_label)
                func = lambda x, pos: "{:g}".format(x * 1e3)
                fmt = mpl.ticker.FuncFormatter(func)
                cbar.ax.yaxis.set_major_formatter(fmt)
        if standalone:
            plt.tight_layout()
        return im


def plot_image(data, header, titles, vmin=0.1, vmax=10):
    # Setup figure with gridspec
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.05], wspace=0)

    if not vmin:
        vmin = min([np.nanmin(imdat) for imdat in data])
    if not vmax:
        vmax = max([np.nanmax(imdat) for imdat in data])

    for i, (arr, title) in enumerate(zip(data, titles)):
        ax = fig.add_subplot(
            gs[0, 0],
            projection=wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL]),
        )
        im = ax.imshow(arr, norm=Normalize(vmin=vmin, vmax=vmax))
        ax.set_title(title)

        ax.set_ylabel("Dec")
        ax.set_xlabel("RA")

    cbar_ax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label("mJy/PSF")
    func = lambda x, _: "{:g}".format(x)
    fmt = mpl.ticker.FuncFormatter(func)
    cbar.ax.yaxis.set_major_formatter(fmt)

    fig.tight_layout()
    return fig


def plot_all_images(data, header, titles, shape=(2, 2), vmin=None, vmax=None):
    # Setup figure with gridspec
    fig = plt.figure(figsize=(4 * shape[1], 3 * shape[0]))
    nrows, ncols = shape
    gs = gridspec.GridSpec(
        nrows, ncols + 1, width_ratios=[1] * ncols + [0.05], wspace=0, hspace=0.2
    )

    if not vmin:
        vmin = min([np.nanmin(imdat) for imdat in data])
    if not vmax:
        vmax = max([np.nanmax(imdat) for imdat in data])

    images = []
    for i, (arr, title) in enumerate(zip(data, titles)):
        ax = fig.add_subplot(
            gs[i // ncols, i % ncols],
            projection=wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL]),
        )
        im = ax.imshow(arr, norm=Normalize(vmin=vmin, vmax=vmax))
        # ax.set_title(title)
        images.append(im)

        if i % ncols != 0:
            ax.set_yticks([])
            ax.set_ylabel("")
        else:
            ax.set_ylabel("Dec")
        if i // ncols != nrows - 1:
            ax.set_xticks([])
            ax.set_xlabel("")
        else:
            ax.set_xlabel("R")

    cbar_ax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(images[0], cax=cbar_ax)
    cbar.set_label("mJy/PSF")
    func = lambda x, _: "{:g}".format(x)
    fmt = mpl.ticker.FuncFormatter(func)
    cbar.ax.yaxis.set_major_formatter(fmt)

    fig.tight_layout()
    return fig


@click.group()
def main():
    """plot DP3 solutions ..."""


@main.command("plot")
@click.option(
    "--images",
    "-i",
    multiple=True,
    type=click.STRING,
    help="images",
)
@click.option("--imagelist", help="imagelist", type=str, default="")
@click.option("--filename", help="output filename", type=str, default="images")
@click.option(
    "--vmin",
    help="vmin",
    type=float,
)
@click.option(
    "--vmax",
    help="vmax",
    type=float,
)
def plot(images="", imagelist="", filename="images", vmin=None, vmax=None):
    if imagelist:
        with open(imagelist) as t:
            images = [line.strip() for line in t]
    else:
        assert images

    if len(images) == 1 and ("*" in images[0] or "?" in images[0]):
        images = sorted(glob(images[0]))

    click.echo(f"Processing {len(images)} files")
    click.echo("images")

    data = [Image(fits_file).data_array for fits_file in images]
    titles = [str(idx).zfill(3) for idx in range(len(images))]
    header = Image(images[0]).header
    # fig = plot_all_images(data, header, titles, shape=(7, 10), vmin=vmin, vmax=vmax)
    fig = plot_image(data, header, titles, vmin=None, vmax=None)
    fig.savefig(f"{filename}_images.png", bbox_inches="tight")


if __name__ == "__main__":
    main()
