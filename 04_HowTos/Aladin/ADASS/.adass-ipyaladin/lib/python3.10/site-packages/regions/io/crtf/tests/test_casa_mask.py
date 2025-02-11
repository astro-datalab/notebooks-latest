# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pickle
import tempfile

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

from regions.io.crtf.write import _write_crtf
from regions.shapes.ellipse import EllipseSkyRegion

try:
    from casatasks import tclean
    from casatools import image
    from casatools import measures as me
    from casatools import simulator
    HAS_CASATOOLS = True
except ImportError:
    HAS_CASATOOLS = False


@pytest.mark.skipif(not HAS_CASATOOLS, reason='casatools is required')
def test_casa_masking():
    with tempfile.TemporaryDirectory() as tmpdir:
        # SIMULATE SOME DATA SET

        # Define antennas
        diam = [25, 25, 25, 25, 25]
        xx = [50, 100, 150, 200, 250]
        yy = [2, -5, -20, -50, -100]
        zz = [-0.5, -1.0, -1.5, -2.0, -2.5]

        sm = simulator()
        sm.open(tmpdir + '/SIM.ms')
        # do configuration
        posvla = me.observatory('VLA')
        sm.setconfig(telescopename='VLA', x=xx, y=yy, z=zz, dishdiameter=diam,
                     mount='alt-az', antname='VLA',
                     coordsystem='local', referencelocation=posvla)

        # Initialize the spectral windows
        sm.setspwindow(spwname='CBand', freq='5GHz', deltafreq='50MHz',
                       freqresolution='50MHz', nchannels=1,
                       stokes='RR RL LR LL')

        # Initialize the source and calibrater
        sm.setfield(sourcename='My cal',
                    sourcedirection=['J2000', '00h0m0.0', '+45.0.0.000'],
                    calcode='A')
        sm.setfield(sourcename='My source',
                    sourcedirection=['J2000', '01h0m0.0', '+47.0.0.000'])

        sm.setlimits(shadowlimit=0.001, elevationlimit='8.0deg')
        sm.setauto(autocorrwt=0.0)

        sm.settimes(integrationtime='10s', usehourangle=False,
                    referencetime=me.epoch('utc', 'today'))

        sm.observe('My cal', 'CBand', starttime='720s', stoptime='1020s')
        sm.observe('My source', 'CBand', starttime='1030s', stoptime='1500s')

        sm.close()

        # Create mask to use during clean
        reg = EllipseSkyRegion(SkyCoord(0.0 * u.deg, 45.0 * u.deg,
                                        frame='fk5'),
                               width=1.0 * u.arcmin, height=2.0 * u.arcmin,
                               angle=45 * u.deg)
        _write_crtf([reg], tmpdir + '/SIM.crtf', 'fk5', '.6f', 'deg')

        # Image the dataset
        tclean(vis=tmpdir + '/SIM.ms', imagename=tmpdir + '/SIM', imsize=100,
               cell='5arcsec', niter=1, mask=tmpdir + '/SIM.crtf',
               interactive=False)

        ia = image()
        ia.open(tmpdir + '/SIM.mask')
        mask_array = ia.getregion()
        ia.close()

        with open('data/binary_mask.pkl', 'rb') as f:
            ref_mask = pickle.load(f)

        assert all(mask_array == ref_mask)
