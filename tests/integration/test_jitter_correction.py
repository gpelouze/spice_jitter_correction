"""
Integration test for jitter_correction.py
"""

import copy
import os
import unittest

from astropy.io import fits

import spice_jitter_correction
from spice_jitter_correction.utils import SpiceFilename


class TestSpiceJitterCorrection(unittest.TestCase):

    def setup_filenames(self):
        self.filenames = [
            SpiceFilename('solo_L2_spice-n-ras_20210914T025031_V88_67109159'
                          '-000.fits'),
            ]

    def assert_fits_identical(self, fn_generated, fn_reference):
        diff = fits.FITSDiff(
            fn_generated, fn_reference,
            ignore_keywords=['CHECKSUM'],  # header checksum
            ignore_comments=['CHECKSUM', 'DATASUM'],  # checksums creation time
            )
        self.assertTrue(diff.identical, msg=diff.report())

    def align_file(self, filename, sum_wvl=False):
        filename_L2r = copy.copy(filename)
        if sum_wvl:
            filename_L2r['level'] = 'L2r_quicklook'
        else:
            filename_L2r['level'] = 'L2r'
        filename_reference = filename_L2r.to_str(path='data/L2r_reference')
        filename_generated = filename_L2r.to_str(path='data/L2r_generated')

        # Apply correction
        spice_jitter_correction.correct_jitter(
            filename.to_str(path='data/L2'),
            './data/L2r_generated/',
            sum_wvl=sum_wvl,
            )

        # Verification
        self.assert_fits_identical(filename_reference, filename_generated)

        # Clear results
        if os.path.isfile(filename_generated):
            os.remove(filename_generated)

    def test_integration_full_spectrum(self):
        self.setup_filenames()

        for filename in self.filenames:
            msg = f'file={filename}'
            with self.subTest(msg=msg):
                self.align_file(filename, sum_wvl=False)

    def test_integration_sum_wvl(self):
        self.setup_filenames()

        for filename in self.filenames:
            msg = f'file={filename}'
            with self.subTest(msg=msg):
                self.align_file(filename, sum_wvl=True)


if __name__ == '__main__':
    unittest.main()
