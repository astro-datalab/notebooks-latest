import os

import numpy as np
from astropy.io import fits
from astropy.table import Table

import ipywidgets as widgets
from IPython.display import display, HTML
#
# Bokeh
#
from bokeh.plotting import figure, output_file
from bokeh.io import push_notebook, show, output_notebook
from bokeh.models import Arrow, Band, CustomJS, ColumnDataSource, Label, Legend, Range1d, Slider, Span, VeeHead
from bokeh.layouts import row, column, widgetbox
from bokeh.models.widgets import Div
import bokeh.events
import bokeh.palettes

# from dl import authClient as ac, storeClient as sc
from dl import storeClient as sc

#
# The first bunched set are emission lines from the spZline files.
# The second bunched set are absorption lines.
# See $IDLSPEC2D_DIR/etc/emlines.par
#
# Note: Wavelengths are in air for lambda > 2000, vacuum for lambda < 2000.
#
lines = [
    #
    # Emission lines
    #
    {"name" : "Ly-α",           "lambda" : 1215.67,  "emission": True },
    {"name" : "N V 1240",       "lambda" : 1240.81,  "emission": True },
    {"name" : "C IV 1549",      "lambda" : 1549.48,  "emission": True },
    {"name" : "He II 1640",     "lambda" : 1640.42,  "emission": True },
    {"name" : "C III] 1908",    "lambda" : 1908.734, "emission": True },
    {"name" : "Mg II 2799",     "lambda" : 2799.49,  "emission": True },
    {"name" : "[O II] 3725",    "lambda" : 3726.032, "emission": True },
    {"name" : "[O II] 3727",    "lambda" : 3728.815, "emission": True },
    {"name" : "[Ne III] 3868",  "lambda" : 3868.76,  "emission": True },
    {"name" : "Hζ",             "lambda" : 3889.049, "emission": True },
    {"name" : "[Ne III] 3970",  "lambda" : 3970.00,  "emission": True },
    {"name" : "Hδ",             "lambda" : 4101.734, "emission": True },
    {"name" : "Hγ",             "lambda" : 4340.464, "emission": True },
    {"name" : "[O III] 4363",   "lambda" : 4363.209, "emission": True },
    {"name" : "He II 4685",     "lambda" : 4685.68,  "emission": True },
    {"name" : "Hβ",             "lambda" : 4861.325, "emission": True },
    {"name" : "[O III] 4959",   "lambda" : 4958.911, "emission": True },
    {"name" : "[O III] 5007",   "lambda" : 5006.843, "emission": True },
    {"name" : "He II 5411",     "lambda" : 5411.52,  "emission": True },
    {"name" : "[O I] 5577",     "lambda" : 5577.339, "emission": True },
    {"name" : "[N II] 5755",    "lambda" : 5754.59,  "emission": True },
    {"name" : "He I 5876",      "lambda" : 5875.68,  "emission": True },
    {"name" : "[O I] 6300",     "lambda" : 6300.304, "emission": True },
    {"name" : "[S III] 6312",   "lambda" : 6312.06,  "emission": True },
    {"name" : "[O I] 6363",     "lambda" : 6363.776, "emission": True },
    {"name" : "[N II] 6548",    "lambda" : 6548.05,  "emission": True },
    {"name" : "Hα",             "lambda" : 6562.801, "emission": True },
    {"name" : "[N II] 6583",    "lambda" : 6583.45,  "emission": True },
    {"name" : "[S II] 6716",    "lambda" : 6716.44,  "emission": True },
    {"name" : "[S II] 6730",    "lambda" : 6730.82,  "emission": True },
    {"name" : "[Ar III] 7135",  "lambda" : 7135.790, "emission": True },
    #
    # Absorption lines
    #
    {"name" : "Hζ",             "lambda" : 3889.049, "emission": False },
    {"name" : "K (Ca II 3933)", "lambda" : 3933.7,   "emission": False },
    {"name" : "H (Ca II 3968)", "lambda" : 3968.5,   "emission": False },
    {"name" : "Hε",             "lambda" : 3970.072, "emission": False },
    {"name" : "Hδ",             "lambda" : 4101.734, "emission": False },
    {"name" : "G (Ca I 4307)",  "lambda" : 4307.74,  "emission": False },
    {"name" : "Hγ",             "lambda" : 4340.464, "emission": False },
    {"name" : "Hβ",             "lambda" : 4861.325, "emission": False },
    {"name" : "Mg I 5175",      "lambda" : 5175.0,   "emission": False },
    {"name" : "D2 (Na I 5889)", "lambda" : 5889.95,  "emission": False },
    {"name" : "D1 (Na I 5895)", "lambda" : 5895.92,  "emission": False },
    {"name" : "Hα",             "lambda" : 6562.801, "emission": False },
    ]


def _airtovac(w):
    """Convert air wavelengths to vacuum wavelengths. Don't convert less than 2000 Å.

    Parameters
    ----------
    w : :class:`float`
        Wavelength [Å] of the line in air.

    Returns
    -------
    :class:`float`
        Wavelength [Å] of the line in vacuum.
    """
    if w < 2000.0:
        return w;
    vac = w
    for iter in range(2):
        sigma2 = (1.0e4/vac)*(1.0e4/vac)
        fact = 1.0 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
        vac = w*fact
    return vac

# Mapping of human friendly strings to integers for visual scan results
scan_map = {
    'flag': -1,     #- flag for data expert followup
    'bad':   0,     #- bad target (e.g. low S/N, can't measure z)
    'no':    1,     #- ok data but wrong redshift
    'maybe': 2,     #- redshift might be right
    'yes':   3,     #- redshift definitely is right
}
#- Add reverse lookup (int -> string) to scan_map
scan_names = list(scan_map.keys())
scan_values = [scan_map[name] for name in scan_names]
for _name, _value in zip(scan_names, scan_values):
    scan_map[_value] = _name

#
# Resampling code, borrowed from https://github.com/desihub/desispec/blob/master/py/desispec/interpolation.py
#
def resample_flux(xout, x, flux, ivar=None, extrapolate=False):
    """Returns a flux conserving resampling of an input flux density.
    The total integrated flux is conserved.

    Args:
        - xout: output SORTED vector, not necessarily linearly spaced
        - x: input SORTED vector, not necessarily linearly spaced
        - flux: input flux density dflux/dx sampled at x

    both x and xout must represent the same quantity with the same unit

    Options:
        - ivar: weights for flux; default is unweighted resampling
        - extrapolate: extrapolate using edge values of input array, default is False,
          in which case values outside of input array are set to zero.

    Setting both ivar and extrapolate raises a ValueError because one cannot
    assign an ivar outside of the input data range.

    Returns:
        if ivar is None, returns outflux
        if ivar is not None, returns outflux, outivar
    """

    if ivar is None:
        return _unweighted_resample(xout, x, flux, extrapolate=extrapolate)
    else:
        if extrapolate :
            raise ValueError("Cannot extrapolate ivar. Either set ivar=None and extrapolate=True or the opposite")
        a = _unweighted_resample(xout, x, flux*ivar, extrapolate=False)
        b = _unweighted_resample(xout, x, ivar, extrapolate=False)
        mask = (b>0)
        outflux = np.zeros(a.shape)
        outflux[mask] = a[mask] / b[mask]
        dx = np.gradient(x)
        dxout = np.gradient(xout)
        outivar = _unweighted_resample(xout, x, ivar/dx)*dxout

        return outflux, outivar


def _unweighted_resample(output_x,input_x,input_flux_density, extrapolate=False) :
    """Returns a flux conserving resampling of an input flux density.
    The total integrated flux is conserved.

    Args:
        output_x: SORTED vector, not necessarily linearly spaced
        input_x: SORTED vector, not necessarily linearly spaced
        input_flux_density: input flux density dflux/dx sampled at x

    both must represent the same quantity with the same unit
    input_flux_density =  dflux/dx sampled at input_x

    Options:
        extrapolate: extrapolate using edge values of input array, default is False,
                     in which case values outside of input array are set to zero

    Returns:
        returns output_flux
    """

    # shorter names
    ix=input_x
    iy=input_flux_density
    ox=output_x

    # boundary of output bins
    bins=np.zeros(ox.size+1)
    bins[1:-1]=(ox[:-1]+ox[1:])/2.
    bins[0]=1.5*ox[0]-0.5*ox[1]     # = ox[0]-(ox[1]-ox[0])/2
    bins[-1]=1.5*ox[-1]-0.5*ox[-2]  # = ox[-1]+(ox[-1]-ox[-2])/2

    # make a temporary node array including input nodes and output bin bounds
    # first the boundaries of output bins
    tx=bins.copy()

    # if we do not extrapolate,
    # because the input is a considered a piece-wise linear function, i.e. the sum of triangles f_i(x),
    # we add two points at ixmin = ix[0]-(ix[1]-ix[0]) and  ixmax = ix[-1]+(ix[-1]-ix[-2])
    # with zero flux densities, corresponding to the edges of the first and last triangles.
    # this solves naturally the edge problem.
    if not extrapolate :
        # note we have to keep the array sorted here because we are going to use it for interpolation
        ix = np.append( 2*ix[0]-ix[1] , ix)
        iy = np.append(0.,iy)
        ix = np.append(ix, 2*ix[-1]-ix[-2])
        iy = np.append(iy, 0.)

    # this sets values left and right of input range to first and/or last input values
    # first and last values are=0 if we are not extrapolating
    ty=np.interp(tx,ix,iy)

    #  add input nodes which are inside the node array
    k=np.where((ix>=tx[0])&(ix<=tx[-1]))[0]
    if k.size :
        tx=np.append(tx,ix[k])
        ty=np.append(ty,iy[k])

    # sort this node array
    p = tx.argsort()
    tx=tx[p]
    ty=ty[p]

    # now we do a simple integration in each bin of the piece-wise
    # linear function of the temporary nodes

    # integral of individual trapezes
    trapeze_integrals=(ty[1:]+ty[:-1])*(tx[1:]-tx[:-1])/2.

    # output flux
    # for each bin, we sum the trapeze_integrals that belong to that bin
    # and divide by the bin size

    trapeze_centers=(tx[1:]+tx[:-1])/2.
    binsize = bins[1:]-bins[:-1]

    if np.any(binsize<=0)  :
        raise ValueError("Zero or negative bin size")

    return np.histogram(trapeze_centers, bins=bins, weights=trapeze_integrals)[0] / binsize


class SDSSSpectra(object):
    """Simple container object for SDSS spectra.
    """
    def __init__(self, wave, flux, ivar, sky, plugmap):
        self.bands = ['c']  # Coadded SDSS spectra have one channel
        self.wave = dict(c=wave)
        self.flux = dict(c=flux)
        self.ivar = dict(c=ivar)
        self.extra = dict(c=sky)
        self.fibermap = plugmap

    def wavelength_grid(self, band):
        if band not in self.bands:
            raise KeyError("{} is not a valid band".format(band))
        return self.wave[band]

    def target_ids(self):
        return np.arange(self.num_spectra(), dtype=np.int32)

    def num_spectra(self):
        if self.fibermap is not None:
            return len(self.fibermap)
        return 0

    def num_targets(self):
        if self.fibermap is not None:
            return len(self.fibermap)
        return 0


def _read_templates():
    pass


def load_sdss_spectra(spPlate, spZbest=None):
    '''
    Load spectra and return an Inspector object

    Parameters
    ----------
    spPlate : :class:`str`
        Full path to input spPlate file.
    spZbest : :class:`str`, optional
        Full path to input spZbest file. If not provided,
        the path will be deduced from `spPlate`.

    Returns
    -------
    :class:`Inspector`
        Object loaded with spectra from files.
    '''
    if spZbest is None:
        d, b = os.path.split(spPlate)
        dd, plate = os.path.split(d)
        ddd, run2d = os.path.split(dd)
        if run2d in ('26', '103', '104'):
            spZbest = os.path.join(ddd, run2d, plate, b.replace('spPlate', 'spZbest'))
        else:
            spZbest = os.path.join(ddd, run2d, plate, run2d, b.replace('spPlate', 'spZbest'))
    with fits.open(sc.get(spPlate, mode='fileobj')) as hdulist:
        wave = 10.0**(np.arange(hdulist['PRIMARY'].header['NAXIS1'],
                                dtype=np.float32)*hdulist['PRIMARY'].header['COEFF1'] +
                                hdulist['PRIMARY'].header['COEFF0'])
        flux = hdulist['PRIMARY'].data
        ivar = hdulist['IVAR'].data
        sky = hdulist['SKY'].data
        plugmap = Table(hdulist['PLUGMAP'].data)
    with fits.open(sc.get(spZbest, mode='fileobj')) as hdulist:
        zbest = Table(hdulist[1].data)
        zbest_model = hdulist[2].data

    spectra = SDSSSpectra(wave, flux, ivar, sky, plugmap)

    return Inspector(lines, spectra, zbest, zbest_model)


class Inspector(object):
    """An interface to plotting spectra with Bokeh in a Jupyter notebook"""

    def __init__(self, lines, spectra, zbest, templates=None):
        """Create an Inspector object.

        Parameters
        ----------
        spectra : :class:`desispec.spectra.Spectra` object
        zbest : Table of zbest output from redrock
        """
        self.lines = lines
        self.zbest = zbest
        self.spectra = spectra
        if templates is None:
            self.templates = _read_templates()
        else:
            self.templates = templates
        self.nspec = len(self.zbest)

        # assert np.all(self.spectra.target_ids() == self.zbest['TARGETID'])
        # assert np.all(self.spectra.target_ids() == self.spectra.fibermap['TARGETID'])

        self.data = dict()     #- high resolution
        self.xdata = dict()    #- low resolution
        self.ispec = 0
        self._update_data()
        self._emission = False
        self._absorption = False
        self.print_targets_info()
        self._plotted = False

        #- dictionary for holding results of visual scan
        self.visual_scan = Table(dtype=[
            ('targetid', int),
            ('scanner', 'S16'),
            ('z', float),
            ('spectype', 'S6'),
            ('subtype', 'S6'),
            ('intresult', 'int16'),
            ('result', 'S6')
        ])
        #- Add header keywords for mapping scan names/values
        for value, name in sorted(zip(scan_values, scan_names)):
            key = 'VSCAN{:02d}'.format(value)
            self.visual_scan.meta[key] = name

        output_notebook()

    #- Property accessors for common target properties
    @property
    def z(self):
        """The redshift of the current target."""
        return self.zbest[self.izbest]['Z']

    @property
    def spectype(self):
        """The spectral classification type of the current target."""
        return self.zbest[self.izbest]['SPECTYPE']

    @property
    def targetid(self):
        """The targetid of the current target."""
        return self.zbest[self.izbest]['TARGETID']

    def select(self, targetids, verbose=False):
        '''Filter spectra to only the specified targetids
        '''
        ii = np.in1d(self.zbest['TARGETID'], targetids)
        self.zbest = self.zbest[ii]
        self.spectra = self.spectra.select(targets=targetids)
        self.nspec = len(self.zbest)
        self.ispec = 0
        self._update()
        if verbose:
            self.print_targets_info()

    def print_targets_info(self):
        '''Print information about the targets in this Inspector object'''
        ntargets = self.spectra.num_targets()
        print('{} targets'.format(ntargets), end='')
        fm = self.spectra.fibermap
        nexp = len(fm)
        try:
            nelg = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.ELG)
            nlrg = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.LRG)
            nqso = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.QSO)
            nbgs = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.BGS_ANY)
            nmws = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.MWS_ANY)
            print(' including {} ELG, {} LRG, {} QSO, {} BGS, {} MWS'.format(
                nelg, nlrg, nqso, nbgs, nmws))
        except KeyError:
            pass

    def plot(self):
        '''
        Plot the spectra
        '''

        #- Make notebook use full width of screen
        display(HTML("<style>.container { width:100% !important; }</style>"))

        #-----
        #- Main spectrum plot; use p for shorthand
        #- set dummy y_range that will be updated later
        tools = 'pan,box_zoom,wheel_zoom,crosshair,undo,redo,reset,save'
        self.specplot = p = figure(plot_height=400, plot_width=800,
                        y_range=(-1,1),
                        output_backend="webgl",
                        toolbar_location='above', tools=tools)

        p.toolbar.active_drag = p.tools[0]    #- pan
        p.toolbar.active_scroll = p.tools[2]  #- wheel zoom

        #- Assemble data for the current spectrum
        bands = self.spectra.bands
        colors = dict(b='#1f77b4', r='#d62728', z='maroon', c='#d62728')
        flux_lines = list()
        model_lines = list()
        for channel in bands:
            flux_lines.append(p.line('wave', 'flux',
                source=self.xdata[channel],
                line_color=colors[channel], line_width=1, alpha=1.0))
            flux_lines.append(p.line('wave', 'upper',
                              source=self.xdata[channel],
                              line_color=colors[channel],
                              line_width=1, alpha=0.3))
            flux_lines.append(p.line('wave', 'lower',
                              source=self.xdata[channel],
                              line_color=colors[channel],
                              line_width=1, alpha=0.3))
            model_lines.append(p.line('wave', 'model',
                source=self.xdata[channel],
                line_color='black', line_width=1, alpha=1.0))

        #- Add horizontal line at y=0
        xtmp = [self.xdata[bands[0]].data['wave'][0],
                self.xdata[bands[-1]].data['wave'][-1]]
        ytmp = [0,0]
        p.line(xtmp, ytmp, color='black', alpha=0.5)

        #- main spectrum plot formatting
        p.yaxis.axis_label = 'Flux [10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹]'
        p.xaxis.axis_label = 'Wavelength [Å]'
        p.xaxis.axis_label_text_font_style = 'normal'
        p.yaxis.axis_label_text_font_style = 'normal'
        p.min_border_left = 60
        p.min_border_bottom = 40

        #- Add legend for flux and model lines
        legend = Legend(items=[
            ("flux",  flux_lines),
            ("model", model_lines),
        ])
        p.add_layout(legend, 'center')
        p.legend.click_policy = 'hide'    #- or 'mute'

        #- Unclear why this is needed here, but if it isn't, the toolbar
        #- disappears when it is called later.
        self._update_lines()

        #-----
        #- Zoom plot of wherever the mouse is hovering on main specplot
        #- use pz for shorthand
        self.zoomplot = figure(title=None,
                plot_height=200, plot_width=200,
                y_range=p.y_range, x_range=(5000,5100),
                output_backend="webgl",
                toolbar_location=None, tools=[])
        for channel in bands:
            self.zoomplot.line('wave', 'flux', source=self.data[channel],
                    line_color=colors[channel], line_width=1, line_alpha=1.0)
            self.zoomplot.line('wave', 'upper', source=self.data[channel],
                    line_color=colors[channel], line_width=1, line_alpha=0.3)
            self.zoomplot.line('wave', 'lower', source=self.data[channel],
                    line_color=colors[channel], line_width=1, line_alpha=0.3)
            self.zoomplot.line('wave', 'model', source=self.data[channel],
                    line_color='black', line_width=1, alpha=1.0)

        #- Callback to update zoom window x-range
        def zoom_callback(zoomplot):
            return CustomJS(args=dict(xr=zoomplot.x_range), code="""
                xr.start = cb_obj.x - 100;
                xr.end = cb_obj.x + 100;
            """)

        p.js_on_event(bokeh.events.MouseMove, zoom_callback(self.zoomplot))

        #-----
        #- Imaging cutout of target location
        self.im = figure(plot_width=200, plot_height=200,
                         x_range=(0, 256), y_range=(0, 256),
                         x_axis_location=None, y_axis_location=None,
                         output_backend="webgl",
                         toolbar_location=None, tools=[])
        self.im.min_border_left = 0
        self.im.min_border_right = 0
        self.im.min_border_top = 0
        self.im.min_border_bottom = 0

        #- Unclear why this is needed here, but otherwise the callback
        #- to open the URL upon clicking doesn't work
        self._update_cutout()

        #-----
        #- Text area with targeting info
        self.info_div = Div(text='Hello<br/>There', width=400)

        #-----
        #- Put it all together
        self.plot_handle = show(
            column(
                row(
                    self.specplot,
                    column(self.im, self.zoomplot),
                    ),
                row(widgetbox(self.info_div, width=600),),
                height=550,
                ),
            notebook_handle=True
            )
        # self.plot_handle = show(p, notebook_handle=True)

        #- Update the contents of the plots
        self._plotted = True
        self._update()
        self._add_inspection_buttons()

    def _update(self, ispec=None):
        '''Update the data and plots for target number ispec

        If ispec is None, use self.ispec; otherwise set self.ispec = ispec
        '''
        if ispec is not None:
            self.ispec = ispec

        self._update_data()
        if not self._plotted:
            return

        self._update_xylim()
        self._update_lines()

        zb = self.zbest[self.izbest]
        try:
            title = '{0} z={1:.4f} zwarn={2}'.format(
                zb['SPECTYPE'].strip(), zb['Z'], zb['ZWARN'])
        except KeyError:
            title = '{0} ({1}) z={2:.4f} zwarn={3}'.format(
                zb['CLASS'].strip(), zb['SUBCLASS'].strip(), zb['Z'], zb['ZWARNING'])
        self.specplot.title.text = title

        self._update_info_div()
        self._update_cutout()

        push_notebook(handle=self.plot_handle)

    def _update_data(self, ispec=None):
        '''
        Update the data containers for target number ispec

        If ispec is None, use self.ispec; otherwise set self.ispec = ispec

        updates self.ispec, .izbest, .data, .xdata
        '''
        if ispec is not None:
            self.ispec = ispec

        try:
            targetid = self.spectra.fibermap['TARGETID'][self.ispec]
            self.izbest = np.where(self.zbest['TARGETID']==targetid)[0][0]
        except KeyError:
            # SDSS Spectra are row-aligned
            targetid = self.izbest = self.ispec
        zb = self.zbest[self.izbest]
        sdss_model = False
        try:
            tx = self.templates[(zb['SPECTYPE'], zb['SUBTYPE'])]
            coeff = zb['COEFF'][0:tx.nbasis]
            model = tx.flux.T.dot(coeff).T
        except KeyError:
            model = self.templates[self.ispec, :]
            sdss_model = True
        for channel in self.spectra.bands:
            good = self.spectra.ivar[channel][self.ispec, :] > 0
            wave = self.spectra.wave[channel][good]
            flux = self.spectra.flux[channel][self.ispec, good]
            ivar = self.spectra.ivar[channel][self.ispec, good]
            std = 1.0/np.sqrt(ivar)
            xwave = np.arange(wave[0], wave[-1], 3)
            xflux, xivar = resample_flux(xwave, wave, flux, ivar=ivar, extrapolate=False)
            xstd = 1.0/np.sqrt(xivar)
            if sdss_model:
                xmodel = resample_flux(xwave, wave, model[good])
                rmodel = model[good]
            else:
                xmodel = resample_flux(xwave, tx.wave*(1+zb['Z']), model)
                rmodel = resample_flux(wave, tx.wave*(1+zb['Z']), model)
            #
            # Avoid warnings when the length of columns change by changing all the columns simultaneously.
            #
            data_update = dict(wave=wave, flux=flux, ivar=ivar, upper=(flux + std),
                               lower=(flux - std), model=rmodel)
            xdata_update = dict(wave=xwave, flux=xflux, ivar=xivar, upper=(xflux + xstd),
                                lower=(xflux - xstd), model=xmodel)
            if channel in self.data:
                self.data[channel].data = data_update
                self.xdata[channel].data = xdata_update
            else:
                self.data[channel] = ColumnDataSource(data_update)
                self.xdata[channel] = ColumnDataSource(xdata_update)

    def _update_xylim(self):
        '''Update the spectrum and zoom plots xy limits for current data'''
        ymin = ymax = 0.0
        for channel in self.spectra.bands:
            model = self.data[channel].data['model']
            flux = self.data[channel].data['flux']
            ymax = max(ymax, np.max(model)*1.05)
            ymax = max(ymax, np.percentile(flux, 98))
            ymin = min(ymin, np.percentile(flux, 10))
            ymin = min(0, ymin)

        self.specplot.y_range.start = ymin
        self.specplot.y_range.end = ymax

        self.zoomplot.y_range.start = ymin
        self.zoomplot.y_range.end = ymax
        self.zoomplot.x_range.start = 3727*(1 + self.z) - 100
        self.zoomplot.x_range.end = 3727*(1 + self.z) + 100

    def _update_cutout(self, zoom=13, layer='ls-dr67'):
        """Update image cutout plot.

        Returns URL to full interactive legacysurvey.org/viewer at ra,dec
        for current target
        """

        #- Get ra,dec from new or old format fibermap for current target
        try:
            ra = self.spectra.fibermap[self.ispec]['RA_TARGET']
            dec = self.spectra.fibermap[self.ispec]['DEC_TARGET']
        except KeyError:
            try:
                ra = self.spectra.fibermap[self.ispec]['TARGET_RA']
                dec = self.spectra.fibermap[self.ispec]['TARGET_DEC']
            except KeyError:
                ra = self.spectra.fibermap[self.ispec]['RA']
                dec = self.spectra.fibermap[self.ispec]['DEC']
                layer = 'sdss2'

        #- JPEG cutout URL
        u = "http://legacysurvey.org/viewer/jpeg-cutout?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}".format(ra, dec, zoom, layer)

        #- Full legacysurvey.org viewer URL
        v = "http://legacysurvey.org/viewer/?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}".format(ra, dec, zoom, layer)

        #- Update cutout plot
        img = self.im.image_url([u], 1, 1, 256, 256, anchor='bottom_left')
        radec = 'RA, Dec = {:.4f}, {:.4f}'.format(ra, dec)
        self.im.text(10, 256-30, dict(value=radec),
            text_color='yellow', text_font_size='8pt')
        ### self.im.title.text = radec

        #- Add callback to open legacysurvey.org viewer when clicking cutout
        callback = CustomJS(code="window.open('{}', '_blank');".format(v))
        self.im.js_event_callbacks.clear()
        self.im.js_on_event('tap', callback)

        return v

    def _update_info_div(self):
        '''Update the text div with information about the current target'''
        fibermap = self.spectra.fibermap[self.ispec]
        zb = self.zbest[self.izbest]
        info = list()
        info.append('<table>')
        sdss = False
        try:
            targetid = zb['TARGETID']
            label = 'Target ID'
        except KeyError:
            targetid = zb['SPECOBJID']
            label = 'SpecObjID'
            sdss = True
        info.append('<tr><th>{0}</th><td>{1:d}</td></tr>'.format(label, targetid))
        if sdss:
            info.append('<tr><th>Plate</th><td>{0:d}</td></tr>'.format(plate))
            info.append('<tr><th>MJD</th><td>{0:d}</td></tr>'.format(mjd))
            info.append('<tr><th>Fiber Number</th><td>{0:d}</td></tr>'.format(self.ispec + 1))
            info.append('<tr><th>BOSS_TARGET1</th><td>{0:d}</td></tr>'.format(
                fibermap['BOSS_TARGET1']))
            info.append('<tr><th>BOSS_TARGET2</th><td>{0:d}</td></tr>'.format(
                fibermap['BOSS_TARGET2']))
            info.append('<tr><th>ANCILLARY_TARGET1</th><td>{0:d}</td></tr>'.format(
                fibermap['ANCILLARY_TARGET1']))
            info.append('<tr><th>ANCILLARY_TARGET2</th><td>{0:d}</td></tr>'.format(
                fibermap['ANCILLARY_TARGET2']))
        else:
            info.append('<tr><th>DESI_TARGET</th><td>{0}</td></tr>'.format(
                ' '.join(desi_mask.names(fibermap['DESI_TARGET']))))
            info.append('<tr><th>BGS_TARGET</th><td>{0}</td></tr>'.format(
                ' '.join(bgs_mask.names(fibermap['BGS_TARGET']))))
            info.append('<tr><th>MWS_TARGET</th><td>{0}</td></tr>'.format(
                ' '.join(mws_mask.names(fibermap['MWS_TARGET']))))
        info.append('</table>')

        self.info_div.text = '\n'.join(info)

    #-------------------------------------------------------------------------
    #- Navigation and visual inspection buttons

    def _add_inspection_buttons(self):
        #- Create the button objects
        buttons = list()
        layout = widgets.Layout(width='60px')
        buttons.append(widgets.Button(
            description='prev', tooltip='Go to previous target',
            layout=layout))
        buttons.append(widgets.Button(
            description='flag', tooltip='Flag for more inspection',
            layout=layout, button_style='warning'))
        b = widgets.Button(
            description='bad', tooltip='Bad data (e.g. low-S/N)',
            layout=layout)
        b.style.button_color = 'gold'
        buttons.append(b)
        buttons.append(widgets.Button(
            description='no', tooltip='Redshift is not correct',
            layout=layout, button_style='danger'))
        buttons.append(widgets.Button(
            description='maybe', tooltip='Uncertain if redshift is correct',
            layout=layout, button_style='primary'))
        buttons.append(widgets.Button(
            description='yes', tooltip='Confident that redshift is correct',
            layout=layout, button_style='success'))
        buttons.append(widgets.Button(
            description='next', tooltip='Skip to next target without recording yes/no/maybe',
            layout=layout))
        buttons.append(widgets.Button(description='Emission Lines', button_style='primary',
                                      tooltip='Toggle display of emission lines.'))
        buttons.append(widgets.Button(description='Absorption Lines', button_style='danger',
                                      tooltip='Toggle display of absorption lines.'))

        #- What to do when a button is clicked
        def button_callback(source):
            if source.description == 'prev':
                self.prev()
            elif source.description == 'next':
                self.next()
            elif source.description == 'Emission Lines':
                self.emission()
            elif source.description == 'Absorption Lines':
                self.absorption()
            elif source.description in scan_names:
                try:
                    targetid = self.zbest['TARGETID'][self.izbest]
                except KeyError:
                    targetid = self.zbest['SPECOBJID'][self.izbest]
                z = self.zbest['Z'][self.izbest]
                try:
                    spectype = self.zbest['SPECTYPE'][self.izbest]
                    subtype = self.zbest['SUBTYPE'][self.izbest]
                except KeyError:
                    spectype = self.zbest['CLASS'][self.izbest]
                    subtype = self.zbest['SUBCLASS'][self.izbest]

                #- remove previous result if needed
                if targetid in self.visual_scan['targetid']:
                    ii = np.where(self.visual_scan['targetid'] == targetid)[0]
                    self.visual_scan.remove_rows(ii)

                #- Add new visual scan result
                self.visual_scan.add_row(dict(
                    targetid=targetid,
                    scanner=os.getenv('USER'),
                    z=z,
                    spectype=spectype,
                    subtype=subtype,
                    intresult=scan_map[source.description],
                    result=source.description,
                ))
                self.next()
            else:
                raise ValueError('Unknown button {}'.format(source.description))

        #- Add the callback function to every button
        for b in buttons:
            b.on_click(button_callback)

        #- Display the buttons
        display(widgets.HBox(buttons))

        #- Don't display widget close button; javascript magic code from
        #- https://groups.google.com/forum/#!topic/jupyter/r67iMlSmuEg
        hideclose = "<script>$('.widget-area .prompt .close').hide()</script>"
        display(HTML(hideclose))

    def next(self):
        '''Advance to the next target'''
        if self.ispec+1 < self.nspec:
            self.ispec += 1
        else:
            print('end of targets')
        self._update()

    def prev(self):
        '''Go to the previous target'''
        if self.ispec > 0:
            self.ispec -= 1
        else:
            print('Already at first target')
        self._update()

    #-------------------------------------------------------------------------
    #- Toggling emission and absorption line markers

    def emission(self, toggle=None):
        """Toggle the display of known emission lines.

        Parameters
        ----------
        toggle : :class:`bool`, optional
            ``True`` and ``False`` turn on and off emission lines,
            respectively.  If not set, the state will be set to the
            opposite of the current state.
        """
        if toggle is None:
            self._emission = not self._emission
        else:
            self._emission = bool(toggle)
        self._update_lines()
        push_notebook(handle=self.plot_handle)

    def absorption(self, toggle=None):
        """Toggle the display of known absorption lines.

        Parameters
        ----------
        toggle : :class:`bool`, optional
            ``True`` and ``False`` turn on and off emission lines,
            respectively.  If not set, the state will be set to the
            opposite of the current state.
        """
        if toggle is None:
            self._absorption = not self._absorption
        else:
            self._absorption = bool(toggle)
        self._update_lines()
        push_notebook(handle=self.plot_handle)

    def _update_lines(self, line_size=0.25, line_scale=2.0):
        for i, l in enumerate(self.lines):
            shiftedWave = _airtovac(l['lambda'])*(1.0 + self.z)
            visible = (self._line_in_range(shiftedWave) and
                       ((l['emission'] and self._emission) or
                        (self._absorption and not l['emission'])))
            shiftedWave_y = 0.0
            for channel in self.spectra.bands:
                sign = -1.0
                if l['emission']: sign = 1.0
                y_envelope = self.xdata[channel].data['model'] + sign*line_scale/np.sqrt(self.xdata[channel].data['ivar'])
                if self.xdata[channel].data['wave'].min() < shiftedWave < self.xdata[channel].data['wave'].max():
                    shiftedWave_y = np.interp(shiftedWave,
                                              self.xdata[channel].data['wave'],
                                              y_envelope)
                    break
            if l['emission']:
                lc = 'blue'
                y_start = shiftedWave_y + line_size
                y_end = shiftedWave_y
            else:
                lc = 'red'
                y_start = shiftedWave_y - line_size
                y_end = shiftedWave_y
            if 'span' in l:
                l['source'].data = dict(x_start=[shiftedWave],
                                        y_start=[y_start],
                                        x_end=[shiftedWave],
                                        y_end=[y_end])
                l['span'].visible = visible
                l['label'].x = shiftedWave
                l['label'].y = y_start
                l['label'].visible = visible
            else:
                l['source'] = ColumnDataSource(data=dict(x_start=[shiftedWave],
                                                         y_start=[y_start],
                                                         x_end=[shiftedWave],
                                                         y_end=[y_end]))
                l['span'] = Arrow(end=VeeHead(size=2,
                                              line_color=lc, line_alpha=0.3,
                                              fill_color=lc, fill_alpha=0.3),
                                  line_color=lc, line_width=2, line_alpha=0.3,
                                  x_start='x_start', y_start='y_start',
                                  x_end='x_end', y_end='y_end',
                                  source=l['source'], visible=visible)
                l['label'] = Label(x=shiftedWave, y=y_start,
                                   text=l['name'], text_color=lc, text_alpha=0.5,
                                   visible=visible)
                self.specplot.add_layout(l['span'])
                self.specplot.add_layout(l['label'])

    def _line_in_range(self, l):
        """True if a spectral line is within the range of the plot.

        Parameters
        ----------
        l : :class:`float`
            Wavelength [Å] of the line to be tested.

        Returns
        -------
        :class:`bool`
            ``True`` if the line should be plotted.
        """
        bands = self.spectra.bands
        return self.xdata[bands[0]].data['wave'].min() < l < self.xdata[bands[-1]].data['wave'].max()
