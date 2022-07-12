import os
import re


class SpiceFilename(dict):
    re_spice_L123_filename = re.compile(
        r'''
        solo
        _(?P<level>L[123])
        _spice
            (?P<concat>-concat)?
            -(?P<slit>[wn])
            -(?P<type>(ras|sit|exp))
            (?P<db>-db)?
            (?P<int>-int)?
        _(?P<time>\d{8}T\d{6})
        _(?P<version>V\d{2})
        _(?P<SPIOBSID>\d+)-(?P<RASTERNO>\d+)
        (?P<ext>\.fits)
        ''',
        re.VERBOSE
        )

    def __init__(self, filename):
        """ Represent SPICE filename

        Parameters
        ----------
        filename : str
            eg: /a/b/solo_L2_spice-n-ras_20210914T025031_V03_67109159-000.fits
        """
        super().__init__(self._parse_filename(filename))
        self._original_filename = filename

        self._check_filename_reconstruction()

    def to_str(self, **kw):
        """ Generate filename string

        Parameters
        ----------
        **kw :
            If specified, these keywords are replaced in the generated filename
            If not specified, the default values are used. Supported
            keywords are:
            - level ('L1', 'L2', ...)
            - concat ('-concat' or '')
            - slit ('n' or 'w')
            - type ('ras', 'sit', 'exp')
            - db ('-db' or '')
            - int ('-int' or '')
            - time (YYYYMMDDTHHMMSS)
            - version ('V03', ...)
            - SPIOBSID ('67109159', ...)
            - RASTERNO ('000', ...)
            - ext ('.fits')

        Returns
        -------
        filename : str
            New filename with interpolated keywords.
        """
        if kw:
            props = self.copy()
            props.update(**kw)
        else:
            props = self
        basename = ('solo_{level}_spice{concat}-{slit}-{type}{db}{int}'
                    '_{time}_{version}_{SPIOBSID}-{RASTERNO}{ext}')
        filename = os.path.join(
            props['path'],
            basename.format(**props),
            )
        return filename

    def _parse_filename(self, filename):
        path, basename = os.path.split(filename)
        m = self.re_spice_L123_filename.match(basename)
        if m is None:
            raise ValueError(f'could not parse SPICE filename: {basename}')
        d = m.groupdict()
        for k, v in d.items():
            if v is None:
                d[k] = ''
        d['path'] = path
        return d

    def _check_filename_reconstruction(self):
        filename = self.to_str()
        msg = (f"Filename reconstruction error: "
               f"{self._original_filename} -> {filename}")
        assert filename == self._original_filename, msg

    def __repr__(self):
        return f'<{type(self).__name__} {self}>'

    def __str__(self):
        return self.to_str()
