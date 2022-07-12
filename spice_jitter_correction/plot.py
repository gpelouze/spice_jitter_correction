from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.visualization as viz
import numpy as np


class PlotPointing:

    @staticmethod
    def plot_hdu(hdu, ax, spice_name):
        img = hdu.data
        if img.ndim > 2:
            img = np.nansum(img, axis=1)  # Sum over wavelengths
            img = np.squeeze(img)  # Collapse 1-depth axis (t or X)

        norm = viz.ImageNormalize(
            img,
            interval=viz.AsymmetricPercentileInterval(1, 99),
            stretch=viz.LogStretch(a=9),  # log10
            )
        ax.imshow(
            img,
            origin='lower',
            norm=norm,
            aspect=hdu.header['CDELT2'] / hdu.header['CDELT1'],
            )
        plt.title(f'{spice_name}\n{hdu.name}', fontsize=12)
        plt.xlabel('X [px]')
        plt.ylabel('Y [px]')
        plt.gca().ticklabel_format(scilimits=(-3, 4))

    def plot_hdulist(self, hdulist, spice_name, filename):
        with PdfPages(filename) as pdf:
            for hdu in hdulist:
                if hdu.is_image and hdu.name != 'WCSDVARR':
                    plt.clf()
                    self.plot_hdu(hdu, plt.gca(), spice_name)
                    pdf.savefig()
