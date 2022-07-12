import argparse

from .jitter_correction import correct_jitter


def main():
    p = argparse.ArgumentParser(
        description=('Correct the pointing of SPICE for each slit position, '
                     'using the AOCS data present in the FITS files'),
        )
    p.add_argument(
        'files', nargs='+',
        help='SPICE L2 FITS'
        )
    p.add_argument(
        '-O', '--output-dir', required=True,
        help='output directory (required)'
        )
    p.add_argument(
        '-p', '--plot-results', action='store_true',
        help='generate plots to visualize the results'
        )
    p.add_argument(
        '--sum-wvl', action='store_true',
        help='save wavelength-integrated images'
        )
    args = p.parse_args()

    for filename in args.files:
        correct_jitter(
            filename,
            args.output_dir,
            plot_results=args.plot_results,
            sum_wvl=args.sum_wvl,
            )
