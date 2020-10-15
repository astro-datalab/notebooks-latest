import numpy as np

from bokeh.layouts import row, column
from bokeh.models import CustomJS, Slider
from bokeh.plotting import figure, output_file, show, ColumnDataSource

x = np.linspace(0, 10, 500)
y = np.sin(x + 0.5)
x2 = np.array([0.0,np.cos(0.5)])
y2 = np.array([0.0,np.sin(0.5)])

source = ColumnDataSource(data=dict(x=x, y=y))
source_polar = ColumnDataSource(data=dict(x2=x2,y2=y2))

plot1 = figure(y_range=(-10, 10), plot_width=400, plot_height=400)

plot1.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

plot2 = figure(x_range=(-10, 10),y_range=(-10, 10), plot_width=400, plot_height=400)
plot2.line('x2', 'y2', source=source_polar, line_width=3, line_alpha=0.6)

amp_slider = Slider(start=0.1, end=10, value=1, step=.1, title="Amplitude")
freq_slider = Slider(start=0.1, end=10, value=1, step=.1, title="Frequency")
phase_slider = Slider(start=0.0, end=6.4, value=0.5, step=.1, title="Phase")
offset_slider = Slider(start=-5, end=5, value=0, step=.1, title="Offset")

callback = CustomJS(args=dict(source=source, source_polar=source_polar, amp=amp_slider, freq=freq_slider, phase=phase_slider, offset=offset_slider),
                    code="""
    const data = source.data;
    const data_polar = source_polar.data;
    const A = amp.value;
    const k = freq.value;
    const phi = phase.value;
    const B = offset.value;
    const x = data['x']
    const y = data['y']
    const x2 = data_polar['x2']
    const y2 = data_polar['y2']
    
    for (var i = 0; i < x.length; i++) {
        y[i] = B + A*Math.sin(k*x[i]+phi);
    }
    x2[1] = A*Math.cos(phi)
    y2[1] = A*Math.sin(phi) + B
    y2[0] = B
    
    source_polar.change.emit();
    source.change.emit();
""")


amp_slider.js_on_change('value', callback)
freq_slider.js_on_change('value', callback)
phase_slider.js_on_change('value', callback)
offset_slider.js_on_change('value', callback)

layout = row(
    plot1,plot2,
    column(amp_slider, freq_slider, phase_slider, offset_slider),
)

output_file("slider.html", title="slider.py example", mode='inline')

show(layout)
