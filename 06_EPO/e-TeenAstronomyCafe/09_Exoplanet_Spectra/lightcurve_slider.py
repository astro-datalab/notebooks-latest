import requests
import numpy as np
import sys
from bokeh.layouts import row, column
from bokeh.models import CustomJS, Slider, Text
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.io import output_notebook
from bokeh.palettes import Spectral
#from bokeh.models import LinearColorMapper
import pdb
import warnings
from json import JSONEncoder
import os
try:
    from astropy.io import fits, ascii
    from astropy.table import Table
except ImportError:
    warnings.warn("Could not find astropy. Data plotter may not work")

if sys.version_info < (3,5):
    warnings.warn("Use a Python 3.5 or later for best results")

axes_font_size = "14pt"

output_notebook()
from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

def initial_imports():
	
	#Next, we download the files needed for the activity.
	url = 'https://raw.githubusercontent.com/davidvargasmora/notebooks-latest/Colab-Badges/06_EPO/e-TeenAstronomyCafe/09_Exoplanet_Spectra'

	r = requests.get(url +'/lc_functions.js', allow_redirects=True, stream=True)
	open('lc_functions.js', 'wb').write(r.content)

	r = requests.get(url +'/scattering_functions.js', allow_redirects=True, stream=True)
	open('scattering_functions.js', 'wb').write(r.content)

	r = requests.get(url +'/transmission_spec_functions.js', allow_redirects=True, stream=True)
	open('transmission_spec_functions.js', 'wb').write(r.content)

	r = requests.get(url +'/data/mystery_lc_1.fits', allow_redirects=True, stream=True)
	open('mystery_lc_1.fits', 'wb').write(r.content)

	r = requests.get(url +'/data/mystery_lc_2.fits', allow_redirects=True, stream=True)
	open('mystery_lc_2.fits', 'wb').write(r.content)

	r = requests.get(url +'/data/mystery_lc_3.fits', allow_redirects=True, stream=True)
	open('mystery_lc_3.fits', 'wb').write(r.content)

	r = requests.get(url +'/data/mystery_lc_4.fits', allow_redirects=True, stream=True)
	open('mystery_lc_4.fits', 'wb').write(r.content)

	r = requests.get(url +'/data/opacity_breakdown_gto_f_hd189733b.fits', allow_redirects=True, stream=True)
	open('opacity_breakdown_gto_f_hd189733b.fits', 'wb').write(r.content)

	return "All files downloaded successfully."


def limb_dark(z,r,u=0.2):
    """ Simple limb darkening law
        Ignores the variations across the planet
    """
    C = 1./ (1. - u/6.)
    
    f = np.zeros_like(z)
    outside_pt = (z >= (1. + r))
    f[outside_pt] = 1.0
    
    inside_pt = (z < 1.0)
    mu = np.sqrt(1. - z[inside_pt]**2)
    f[inside_pt] = C * (1.0 - u * (1. - mu))
    
    intersect_pt = (z >= 1.0) & (z < (1. + r))
    f[intersect_pt] = C * (1. - u)
    
    return f
    

def light_c(t,aOr=6.,b=0.2,r=0.1,p=24.0,u=0.2):
    x = aOr * np.sin(t * np.pi * 2. / p)
    bp = b *  np.cos(t * np.pi * 2. / p)
    z = np.sqrt(x**2 + bp**2)
    Aint = area_intersect(z,r)
    f = 1. - Aint * limb_dark(z,r,u=u)
    return f * 100.
    
def area_intersect(z,r):
    f = np.zeros_like(z)
    outside_pt = (z >= (1. + r))
    f[outside_pt] = 0.0
    
    inside_pt = (z <= (1. - r))
    f[inside_pt] = r**2
    
    intersect_pt = (z > (1.0 -r)) & (z < (1. + r)) 
    if np.sum(intersect_pt) > 0:
        x= (1. - r**2 + z[intersect_pt]**2)/(2. * z[intersect_pt])
        theta1 = np.arccos(x)
        theta2 = np.arccos((z[intersect_pt]-x)/r)
        Aint = theta1 + theta2 * r**2 - np.sqrt(1.0 - x**2) * z[intersect_pt]
        f[intersect_pt] = Aint / np.pi
    return f


def practice_slider():
    """ A simple practice slider for the webpage """
    slider = Slider(start=0, end=10, value=0, step=0.25, title='Current Value',
                    bar_color='black')
    source = ColumnDataSource(data=dict(x=[0.0],y=[0.0],txt=['Move Slider to 5.0'],color=['blue']))
    
    callback = CustomJS(args=dict(source=source,s=slider),
                        code="""
                        const data = source.data;
                        const txt = data['txt']
                        const col = data['color']
                        const s_val = s.value
                        
                        if (s_val == 5.0) {
                            txt[0] = 'Good Job!'
                            col[0] = 'green'
                        } else {
                            txt[0] = 'Move Slider to 5.0'
                            col[0] = 'blue'
                        }
                        source.change.emit();
                        """)
    plot1 = figure(plot_width=350,plot_height=80,x_range=[-1,5],y_range=[-1,2],tools="")
    
    txt = Text(x='x',y='y',text='txt',text_color='color')
    plot1.add_glyph(source,txt)
    
    
    slider.js_on_change('value', callback)
    plot1.toolbar_location = None
    plot1.axis.visible = False
    plot1.xgrid.grid_line_color = None
    plot1.ygrid.grid_line_color = None
    
    layout = column([plot1,slider])
    show(layout)

def lightcurve_slider(free_radius=True,free_impact=False,savePlot=False):
    """
    Lightcurve slider to show lightcurve and projected view
    """
    
    
    x = np.linspace(-1.5,1.5,512) ## time (hours)
    y = light_c(x)#np.zeros_like(x) ## flux
    r = [1.0] ## planet radius
    marker_size = [2.0] ## size of time marker
    time_now = [0.0] ## time of interest
    flux_now = light_c(np.array(time_now)) ## flux of interest
    marker_size = [10.0] ## marker size

    xCircle = [0.0]
    yCircle = [2.0]


    source = ColumnDataSource(data=dict(x=x, y=y))
    planet_dict = dict(r=r,x=xCircle,y=yCircle,time_now=time_now,flux_now=flux_now,marker_size=marker_size)
    source_planet = ColumnDataSource(data=planet_dict)

    plot1 = figure(y_range=(97.5, 100.2), plot_width=400, plot_height=200,tools="")

    plot1.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
    plot1.circle('time_now','flux_now',size='marker_size',source=source_planet,color='green')

    plot1.title.text = 'Lightcurve'
    plot1.xaxis.axis_label = "Time from Central Transit (hours)"
    plot1.yaxis.axis_label = "Brightness (%)"
    plot1.xaxis.axis_label_text_font_size = axes_font_size
    plot1.yaxis.axis_label_text_font_size = axes_font_size

    plot2 = figure(x_range=(-20, 20),y_range=(-20, 20), plot_width=400, plot_height=400,tools="")

    
    ## make a limb darkened star
    r_star = 10.0
    u_linear = 0.2 ## linear limb darkening parameter
    img_res = 256
    x_linear = np.linspace(-r_star * 2, r_star * 2, img_res)
    y_linear = np.linspace(-r_star * 2, r_star * 2, img_res)
    xx_grid, yy_grid = np.meshgrid(x_linear, y_linear)
    rr_grid = np.sqrt(xx_grid**2 + yy_grid**2) ## radius
    in_points = rr_grid < r_star ## only the points inside will be calculated
    mu = np.sqrt(r_star**2 - rr_grid[in_points]**2) ##mu
    f_star = np.zeros_like(rr_grid)
    f_star[in_points] = (1.0 - u_linear * (1.0 - mu)) / (1.0 - u_linear / 6.0)
    plot2.image(image=[f_star], x=-2 * r_star, y=-2 * r_star, dw=4 * r_star, dh=4 * r_star, palette="Inferno256", level="image")


    #plot2.circle([0],[0],radius=10,color='yellow')

    plot2.circle('x','y',radius='r',source=source_planet,color='black',line_color='cyan')
    #plot2.line('x2', 'y2', source=source_polar, line_width=3, line_alpha=0.6)
    plot2.xgrid.visible = False
    plot2.ygrid.visible = False
    plot2.xaxis.axis_label = "X Distance (Earth Radii)"
    plot2.yaxis.axis_label = "Y Distance (Earth Radii)"
    plot2.xaxis.axis_label_text_font_size = axes_font_size
    plot2.yaxis.axis_label_text_font_size = axes_font_size
    plot2.title.text = 'Star View'

    t_slider = Slider(start=-1.5, end=1.5, value=0, step=0.01, title='Time from Central Transit (hours)',
                      bar_color='black')
    r_slider = Slider(start=0.0, end=1.5, value=r[0], step=.01, title="Radius (Earth Radii)",
                      bar_color='black')
    b_slider = Slider(start=0.0, end=1.1, value=0.2, step=0.01, title="Impact Parameter",
                      bar_color='black')
    
    sliderList = [t_slider]
    if free_radius == True:
        sliderList.append(r_slider)
    if free_impact == True:
        sliderList.append(b_slider)
    
    with open ("lc_functions.js", "r") as js_file:
        js_code = js_file.read()

    js_args = dict(source=source, source_planet=source_planet, r=r_slider,t=t_slider,b_imp=b_slider)
    callback = CustomJS(args=js_args,
                        code=js_code)
    #    
    
    for oneSlider in sliderList:
        oneSlider.js_on_change('value', callback)
    
    ## Remove the toolbars
    plot1.toolbar_location = None
    plot2.toolbar_location = None

    layout = row(
        column(plot1,plot2),
        column(sliderList),
    )
    
    if savePlot == True:
        outName = "plots/slider_free_rad_{}_free_b_{}.html".format(free_radius,free_impact)
        output_file(outName, title="Radius Slider", mode='inline')

    show(layout)

w0 = 0.67

def calc_radii(w,wRange,thickness=0.3):
    """
    Simple function that converts an "atmospheric thickness" to a radius spectrum
    """
    rad = 0.8 - 1.0 * thickness * (w - w0) / wRange
    return rad

def scattering_slider(savePlot=False,plots=['planet','spectrum','lightcurve']):
    """
    Slider shows the planet, spectrum and lightcurves
    """
    w = np.array([  0.64   ,  0.61    , 0.57  ,   0.53   , 0.47     , 0.41 ])
    nWave = len(w)
    posx = np.zeros_like(w)
    posy = np.zeros_like(w)
    wRange = w[0] - w[-1]
    #colors_array = np.array([  'red' ,'orange','yellow' ,'green',  'blue',  'violet'])
    if nWave <= 11:
        colors_array = np.flip(Spectral[nWave])
    else:
        raise Exception("Too many wavelengths for palette")
    
    thickness = 0.3 ## "atmospheric thickness"
    
    rad_arr = calc_radii(w,wRange,thickness)
    
    source = ColumnDataSource(data=dict(w=w, rad=rad_arr,posx=posx,posy=posy,colors=colors_array))
    
    plot1 = figure(x_range=(-1.3,1.3),y_range=(-1.3,1.3), plot_width=400, plot_height=400,tools="")
    
    plot1.scatter('posx','posy',radius='rad',source=source, line_width=3,
                  fill_color=None,line_color='colors')
    plot1.circle(0.0,0.0,radius=0.8,color='black')
    
    plot1.title.text = 'Planet View'
    plot1.xaxis.axis_label = "X Size (Earth Radii)"
    plot1.yaxis.axis_label = "Y Size (Earth Radii)"
    plot1.xaxis.axis_label_text_font_size = axes_font_size
    plot1.yaxis.axis_label_text_font_size = axes_font_size

    plot2 = figure(y_range=[0.77,1.15],plot_width=400, plot_height=400,tools="")
    plot2.line('w','rad',source=source)
    plot2.xaxis.axis_label = "Wavelength (microns)"
    plot2.yaxis.axis_label = "Radius (Earth Radii)"
    plot2.xaxis.axis_label_text_font_size = axes_font_size
    plot2.yaxis.axis_label_text_font_size = axes_font_size
    plot2.scatter('w','rad',source=source,line_width=None,fill_color='colors',size=12)
    plot2.title.text = 'Spectrum Plot'
    
    
    time = np.linspace(-1.2,1.2,256)
    # fluxData = np.zeros([len(time),len(rad_arr)])
    # for waveInd in np.arange(nWave):
    #     fluxData[:,waveInd] = light_c(time,r=rad_arr[waveInd])
    #lc_dict = {'t': time,'f': fluxData}
    # plot3.line('t','f',source=source_lc)
    
    lc_dict = {'t': time}
    for waveInd in np.arange(nWave):
        lc_dict['f {}'.format(waveInd)] = light_c(time,r=rad_arr[waveInd]/10.)
    
    source_lc = ColumnDataSource(data=lc_dict)
    
    plot3 = figure(y_range=[98.5,100.1],plot_width=400, plot_height=400,tools="")
    
    for waveInd in np.arange(nWave):
        plot3.line('t','f {}'.format(waveInd),source=source_lc,
                   color=colors_array[waveInd],line_width=3)
    plot3.xaxis.axis_label = "Time (hours)"
    plot3.yaxis.axis_label = "Relative Brightness (%)"
    plot3.xaxis.axis_label_text_font_size = axes_font_size
    plot3.yaxis.axis_label_text_font_size = axes_font_size
    plot3.title.text = 'Lightcurve Plot'
    
    t_slider = Slider(start=0, end=0.3, value=0.3, step=0.01, title='Atmospheric Thickness')
        
    with open ("scattering_functions.js", "r") as js_file:
        js_code = js_file.read()
    
    js_args = dict(source=source, source_lc=source_lc,wRange=wRange,t=t_slider)
    callback = CustomJS(args=js_args,
                        code=js_code)
    
    t_slider.js_on_change('value', callback)
    
    ## Remove the toolbars
    plot1.toolbar_location = None
    plot2.toolbar_location = None
    plot3.toolbar_location = None
    
    leftPlots = []
    rightPlots = [t_slider]
    if 'planet' in plots:
        leftPlots.append(plot1)
    if 'spectrum' in plots:
        leftPlots.append(plot2)
    if 'lightcurve' in plots:
        rightPlots.append(plot3)
    
    layout = row(
        column(leftPlots),
        column(rightPlots),
    )
    
    if savePlot == True:
        output_file("plots/slider_scattering.html", title="Radius Slider", mode='inline')

    show(layout)


def calc_radii(w,wRange,thickness=0.3):
    """
    Simple function that converts an "atmospheric thickness" to a radius spectrum
    """
    rad = 0.8 - 1.0 * thickness * (w - w0) / wRange
    return rad

def transmission_spec_slider(mysteryNum=1,savePlot=False):
    """
    Sliders for the transmission spectrum and their lightcurves
    """
    
    #datName = 'data/mystery_lc_{}.fits'.format(mysteryNum)
    datName = 'mystery_lc_{}.fits'.format(mysteryNum)
    if os.path.exists(datName):
        HDUList = fits.open(datName)
    else:
        raise Exception("Mystery {} not found".format(mysteryNum))
    
    w = HDUList['WAVE'].data
    lcData = HDUList['FLUX'].data
    rad_init = 1.5
    rad_arr = np.ones_like(w) * rad_init
    orig_time = HDUList['TIME'].data
    
    time = np.linspace(-1.0,1.0,256) ## time
    
    nWave = len(w)
    
    if nWave <= 11:
        colors = Spectral[nWave]
    else:
        raise Exception("Need to figure out colors for {} wavelengths".format(nWave))
    
    
    source = ColumnDataSource(data=dict(w=w, rad=rad_arr,colors=colors))
    
    lc_dict = {'t': time}
    lc_data = {'t': orig_time}
    for waveInd in np.arange(nWave):
        offset = 1.5 * waveInd
        lc = light_c(time,aOr=15.,b=0.785,r=rad_arr[waveInd]/10.,p=48.0,u=0.2)
        lc_dict['f {}'.format(waveInd)] = lc - offset
        lc_data['f {}'.format(waveInd)] = lcData[waveInd,:] * 100. - offset
    
    source_lc = ColumnDataSource(data=lc_dict)
    source_data = ColumnDataSource(data=lc_data)
    
    plot1 = figure(y_range=[82,100.1],plot_width=400, plot_height=600,tools='')
    
    for waveInd in np.arange(nWave):
        plot1.scatter('t','f {}'.format(waveInd),source=source_data,
                      color=colors[waveInd],line_color='black')
        plot1.line('t','f {}'.format(waveInd),source=source_lc,
                   line_width=3,color=colors[waveInd])
    
    plot1.xaxis.axis_label = "Time (hours)"
    plot1.yaxis.axis_label = "Relative Brightness (%) - Offset"
    plot1.xaxis.axis_label_text_font_size = axes_font_size
    plot1.yaxis.axis_label_text_font_size = axes_font_size
    plot1.title.text = 'Lightcurve Plot'
    
    
    plot2 = figure(y_range=[1.0,2.0],plot_width=400, plot_height=300,tools='')
    plot2.line('w','rad',source=source,color='black',line_width=3)
    plot2.xaxis.axis_label = "Wavelength (microns)"
    plot2.yaxis.axis_label = "Radius (Earth Radii)"
    plot2.xaxis.axis_label_text_font_size = axes_font_size
    plot2.yaxis.axis_label_text_font_size = axes_font_size
    plot2.square('w','rad',source=source,line_width=None,fill_color='colors',size=16)
    plot2.title.text = 'Spectrum Plot'
    
    slider_list = []
    for waveInd in np.arange(nWave):
        thisTitle = "Radius (Earth Radii) at {:.2f} microns".format(w[waveInd])
        r_slider = Slider(start=1.0, end=2.0, value=rad_init, step=0.01, title=thisTitle,
                          bar_color=colors[waveInd])
        slider_list.append(r_slider)
    
    with open ("transmission_spec_functions.js", "r") as js_file:
        js_code = js_file.read()

    js_args = dict(source=source, source_lc=source_lc,r_slider_list = slider_list)
    callback = CustomJS(args=js_args,
                        code=js_code)

    for one_slider in slider_list:
        one_slider.js_on_change('value', callback)
    
    ## Remove the toolbars
    plot1.toolbar_location = None
    plot2.toolbar_location = None
    
    layout = row(
        column(plot1,plot2),
        column(slider_list),
    )
    
    if savePlot == True:
        output_file("plots/slider_transmission.html", title="Transmission Spectrum Slider", mode='inline')

    show(layout)

convertDict = {'H2O':'Water Vapor','CH4':'Methane','CO2':'Carbon Dioxide','Cloudy':'Cloudy'}

def example_spectra(atmospheres=['H2O','CH4','CO2','No Atmosphere'],savePlot=False):

    ## dat = Table.read('data/opacity_breakdown_gto_f_hd189733b.fits')
    ## Changed the route to facilitate working with Colab.

    dat = Table.read('opacity_breakdown_gto_f_hd189733b.fits')
    dat['No Atmosphere'] = 0.0179
    
    plotList = []
    for ind, atmosphere in enumerate(atmospheres):
        if len(atmospheres) > 1:
            plot_width=180
            plot_height=180
        else:
            plot_width=500
            plot_height=400
        
        plot1 = figure(plot_width=plot_width,plot_height=plot_height,tools='',
                      x_range=[2.3,5.1],y_range=[1.3,2.0])
        rad = np.sqrt(dat[atmosphere]) * 10.
        rad = (rad - np.mean(rad)) * 60. + np.mean(rad) ## exaggerate to see better
        plot1.line(dat['Wave'],rad,line_width=4)
        
        if atmosphere in list(convertDict.keys()):
            thisLabel = convertDict[atmosphere]
        else:
            thisLabel = atmosphere
        
        plot1.title.text = thisLabel
        plot1.xaxis.axis_label = 'Wavelength (microns)'
        plot1.yaxis.axis_label = 'Size (Earth Radii)'
        plot1.toolbar_location = None
        
        plotList.append(plot1)
    
    layout = row(plotList)
    
    if savePlot == True:
        output_file("plots/template_spectra.html",title="Template Spectra",mode="inline")
    
    show(layout)


if __name__ == "__main__":
    lightcurve_slider()
