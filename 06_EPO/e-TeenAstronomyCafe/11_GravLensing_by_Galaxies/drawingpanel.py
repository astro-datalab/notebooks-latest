# drawingpanel.py
# =====================================================================
# Simplified Python drawing window class
# to accompany Building Python Programs textbook and associated materials.
# This library is based on Python's Tkinter GUI system.
# see also: http://effbot.org/tkinterbook/canvas.htm
# =====================================================================
# authors:
# Marty Stepp, Stanford University
# Allison Obourn, University of Arizona
# Stuart Reges, University of Washington
# =====================================================================
# To generate the documentation for this library, run this command:
# epydoc --html --no-imports --no-private --no-sourcecode --no-frames -o doc DrawingPanel.py
# =====================================================================
# History:
# 2017/10/15
# - added draw/fill_polygon that accept tuples for points
# - added is/set_resizable
# 2017/09/04
# - added key/mouse event listener functionality
# - added draw/fill_arc
# 2017/09/03
# - added draw_polyline
# - fixed draw/fill_polygon
# - improved PyDoc documentation of parameters
# 2017/09/01
# - separated out Color class
# - refactored fill_* functions
# 2017/08/30
# - center(), set_location(), set_size() functions added
# - DrawingPanel appears centered on screen by default
# 2017/01/27
# - draw_xxx functions added
# - get_pixel_color added
# - save added
# 2009/10/21
# - patch for Python 3 Tkinter as suggested by Steve Geluso
# 2009/01/01
# - initial version
# 
# !! Note to author: also update 'version' in About text !!
# =====================================================================
# COMPATIBILITY NOTE: This library generally should work with Python 2 or 3+.
# If you discover functionality that fails on Python 2, please notify the author.
# =====================================================================
# TODO: support other image formats (BMP, PNG?) for draw and/or save
# - PIL (vs Pillow?)
# - http://effbot.org/tkinterbook/photoimage.htm
# - https://stackoverflow.com/questions/14050281/how-to-check-if-a-python-module-exists-without-importing-it
# TODO: find_font?
# - https://stackoverflow.com/questions/44478807/getting-the-default-font-in-tkinter
# TODO: global exception handler for TclError?
# - https://stackoverflow.com/questions/6598053/python-global-exception-handling
# TODO: make Color into Enum?
# - https://docs.python.org/3/library/enum.html
# See also:

import atexit
import os
import sys
import time
import types

# python3.0 uses "tkinter," python2.x uses "Tkinter."
# detect version info and import appropriate graphics
# code added by UW CSE TA Steve Geluso
if (sys.version_info >= (3, 0)):
    from tkinter import Canvas, Label, Menu, PhotoImage, Tk, TclError
    from tkinter import BOTTOM, LEFT, RIGHT, TOP, W, NW, SW
    from tkinter import filedialog
    from tkinter import messagebox
else:
    from Tkinter import Tk, Canvas, Menu, PhotoImage, TclError
    import tkFileDialog
    import tkMessageBox

class DrawingPanel:
    '''
    A DrawingPanel object represents a simple interface for creating
    graphical windows in Python for drawing shapes, lines, images, and colors.
    The DrawingPanel also supports getting and setting the RGB values of individual
    pixels, which makes it useful for learning about various looping and 2D list-
    processing algorithms.

    See PyDoc comments for individual functions below for more information.

    This library is based on Python's Tkinter GUI system.
    In particular, the drawing surface used by DrawingPanel is a Tkinter Canvas object,
    and many of our functions are wrappers around Tkinter Canvas functions.
    
    @author: Marty Stepp, Stanford University
    @author: Allison Obourn, University of Arizona
    @author: Stuart Reges, University of Washington

    @version: 2017/10/15
    
    @copyright: Marty Stepp, Allison Obourn, Stuart Reges; for individual, educational, non-commercial, private use only; all other rights reserved.

    @see: http://effbot.org/tkinterbook/canvas.htm
    '''

    @staticmethod
    def __adjust_coordinate(coord):
        '''
        Converts between a single 0-based client coordinate and a 1-based Tkinter canvas coordinate.
        '''
        # return coord + 1
        return coord
    
    @staticmethod
    def __adjust_coordinates(*coords):
        '''
        Converts between a list/tuple of 0-based client coordinates and one of 1-based Tkinter canvas coordinates.
        '''
        ADJUST_AMOUNT = 0
        if coords is None or len(coords) == 0:
            raise ValueError
        elif isinstance(coords[0], tuple):
            coords_to_use = coords[0]
        else:
            coords_to_use = coords
        if ADJUST_AMOUNT == 0:
            return coords_to_use
        else:
            return tuple([coord + ADJUST_AMOUNT for coord in coords_to_use])
    
    @staticmethod
    def __adjust_mouse_coordinate(coord):
        '''
        Converts between a single 1-based client mouse coordinate and a 0-based code coordinate.
        '''
        # return coord - 1
        return coord
    
    def __init__(self, width=500, height=500, background="white", **options):
        '''
        Constructs a panel of a given width, height, and optional background color.
        
        Parameters:
        @param width: (optional) width of the DrawingPanel central canvas area in pixels (default 500)
        @param height: (optional) height of the DrawingPanel central canvas area in pixels (default 500)
        @param background: (optional) background color of the DrawingPanel canvas area (default "white")
        @param options: (optional) keyword arguments to pass directly to Tk window constructor (see Tkinter documentation)
        '''

        # create window
        self.x = 0
        self.y = 0
        self.stroke = 1
        self.resizable = True
        self.window = Tk(**options)
        self.window.withdraw()
        self.set_size(width, height)
        
        # fill in window/object properties
        self.fill = "black"
        self.outline = "black"
        self.window.title("DrawingPanel")
        self.font = None
        self.images = []     # list of all drawn images
        self.image_map = dict()
        
        # create window graphical canvas
        self.canvas = Canvas(self.window, width = width, height = height, bd = 0, highlightthickness = 0, relief="ridge")
        self.canvas["bg"] = background
        self.canvas.pack(side = TOP)
        self.canvas.bind("<Motion>", lambda event: self.__callback_mouse_moved(event))
        self.window.bind("<Key>", lambda event: self.__callback_key_pressed(event))
        self.window.bind("<Escape>", lambda event: self.close())
        self.window.bind("<Control-s>", lambda event: self.save_as_dialog())
        
        # create bottom status label
        self.status_label = Label(self.window, text=" ", anchor=SW, justify=LEFT)
        self.status_label.pack(side = BOTTOM, anchor=SW)
        
        # create top menu bar
        menubar = Menu(self.window)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Save As...", command=self.save_as_dialog)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.window.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About...", command=self.about)
        menubar.add_cascade(label="Help", menu=helpmenu)
        self.window.config(menu=menubar)
        
        # finish up
        self.window.update_idletasks()
        self.window.update()
        self.center()
        self.window.deiconify()
        
        # hack - runs mainloop on exit if not interactive
        if not hasattr(sys, 'ps1'):
            self.__install_mainloop_hack()

    def about(self):
        '''
        Shows an About message.
        '''
        msg = '''DrawingPanel

Graphical library class to support Building Java Programs textbook

written by: Marty Stepp, Stanford University; Allison Obourn, University of Arizona; Stuart Reges, University of Washington

Version: 2017/09/04

please visit our web site at:

http://www.buildingpythonprograms.com/
'''
        self.message(msg, "About DrawingPanel")

    def __callback_key_pressed(self, event):
        '''
        Private function that is called when keys are pressed.
        Not to be called by students.
        @see: http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm
        @see: http://www.tcl.tk/man/tcl8.4/TkCmd/keysyms.htm
        '''
        # notes: 
        # event.char contains printable character such as "a";
        # event.keysym contains string for key such as "Escape";
        # event.keycode contains ASCII code such as 97
        pass

    def __callback_mouse_moved(self, event):
        '''
        Private function that is called when the mouse is moved.
        Not to be called by students.
        '''
        x = DrawingPanel.__adjust_mouse_coordinate(int(event.x))   # TODO: adjust?
        y = DrawingPanel.__adjust_mouse_coordinate(int(event.y))
        new_text = "(x=%d, y=%d)" % (x, y)
        rgb = self.get_pixel_color_rgb(x, y)
        new_text += ", r=%d, g=%d, b=%d" % rgb
        self.set_status(new_text)
    
    def add_key_listener(self, function, event_type="press"):
        '''
        Attaches the given function to be called when keyboard events occur in this DrawingPanel window.
        @param function: The function to be called; should accept an event as a parameter.
        @param event_type: (optional) What type of mouse event to respond to. The default is key presses,
        but you can pass any of "press" or "release".
        @see: http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm
        '''
        if event_type.lower() == "press":
            event_type = "<Key>"
        elif event_type.lower() == "release":
            event_type = "<Key>"   # TODO
        self.window.bind(event_type, function)
    
    def add_mouse_listener(self, function, event_type="click"):
        '''
        Attaches the given function to be called when mouse events occur in this DrawingPanel window.
        @param function: The function to be called; should accept an event as a parameter.
        @param event_type: (optional) What type of mouse event to respond to. The default is clicks, but you can pass
        any of "click", "doubleclick", "drag", "enter", "exit", "middleclick", "move", "press", "release", "rightclick", or "wheel".
        @see: http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm
        '''
        if event_type.lower() == "click":
            event_type = "<Button-1>"
        elif event_type.lower() == "doubleclick":
            event_type = "<Double-Button-1>"
        elif event_type.lower() == "drag":
            event_type = "<B1-Motion>"
        elif event_type.lower() == "enter":
            event_type = "<Enter>"
        elif event_type.lower() == "exit":
            event_type = "<Leave>"
        elif event_type.lower() == "middleclick":
            event_type = "<Button-2>"
        elif event_type.lower() == "move":
            event_type = "<Motion>"
        elif event_type.lower() == "press":
            event_type = "<Button-1>"
        elif event_type.lower() == "release":
            event_type = "<ButtonRelease>"
        elif event_type.lower() == "rightclick":
            event_type = "<Button-3>"
        elif event_type.lower() == "wheel":
            event_type = "<MouseWheel>"
            self.window.bind("<4>", function)   # wheel up
            self.window.bind("<5>", function)   # wheel down
            return
        self.window.bind(event_type, function)
    
    def center(self):
        '''
        Centers the DrawingPanel window's location on the screen.
        If the DrawingPanel's size is larger than the screen size in either dimension, does not move it.
        '''
        screenwidth = self.window.winfo_screenwidth()
        screenheight = self.window.winfo_screenheight()
        mywidth = self.window.winfo_width()
        myheight = self.window.winfo_height()
        myx = (screenwidth - mywidth) // 2
        myy = (screenheight - myheight) // 2
        if myx >= 0 and myy >= 0:
            self.set_location(myx, myy)

    def clear(self):
        '''
        Erases all shapes from the panel and fills it with its background color.
        '''
        self.canvas.create_rectangle(0, 0, self.width + 2, self.height + 2, \
                outline=self.canvas["bg"], fill=self.canvas["bg"])
        self.images = []
        self.image_map = dict()
        self.set_status("")
    
    def close(self):
        '''
        Closes the DrawingPanel window so it is no longer visible.
        Once it is closed, it cannot be drawn on or shown on screen again during the current run of your program.
        '''
        self.window.destroy()

    def draw_arc(self, x, y, w, h, start_angle, extent, color="", **options):
        '''
        Draws an outlined arc (a partial oval).
        The resulting arc begins at 'start_angle' and extends for 'extent' degrees.
        Angles are interpreted such that 0 degrees is at the 3 o'clock position.
        A positive value indicates a counter-clockwise rotation while a negative value indicates a clockwise rotation.
        
        If you want to draw a pie-slice style arc rather than just the outer line, pass an option of style="pieslice".
        
        @param x: left x coordinate of bounding box of arc's oval
        @param y: top y coordinate of bounding box of arc's oval
        @param w: width of bounding box of arc's oval
        @param h: height of bounding box of arc's oval
        @param start_angle: starting angle in degrees, where 0 is at "3 o'clock" position
        @param extent: degrees that the arc should extend from its start, where positive numbers are counter-clockwise and negative values are clockwise
        @param color: (optional) color to use for outline of shape (if absent, uses panel's default color as set via set_color)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_arc-method
        '''
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        color = self.__pick_color_with_default(color, self.outline, "outline", options)
        Color.normalize_all_colors(options)
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        if not "start" in options:
            options["start"] = start_angle
        if not "extent" in options:
            options["extent"] = extent
        if not "style" in options:
            options["style"] = "arc"
        bounding_box = (x, y, x + w, y + h)
        self.canvas.create_arc(bounding_box, **options)
    
    def draw_image(self, filename, x, y, w=0, h=0, **options):
        '''
        Draws an image from a file with its top-left corner at the given (x, y) pixel.
        Due to tkinter limitations, can recognize PNG and GIF but not JPG. Sorry.  :-(

        @param filename: path to image file
        @type filename: str
        @param x: top-left x-coordinate at which to place the image
        @param y: top-left y-coordinate at which to place the image
        @param w: (optional) width at which to draw the image (not yet supported)
        @param h: (optional) height at which to draw the image (not yet supported)
        @param options: (optional) keyword arguments to pass directly to Tk create_image function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#canvas.Canvas.create_image-method
        '''
        x = DrawingPanel.__adjust_coordinate(x)
        y = DrawingPanel.__adjust_coordinate(y)
        
        image = PhotoImage(file=filename)
        self.images.append(image)
        
        # image is centered in Python Tkinter; change it so that x/y = top/left corner
        if w > 0 and h > 0:
            # TODO: resize
            w = image.width()
            h = image.height()
        else:
            # use existing size
            w = image.width()
            h = image.height()
        
        result_id = self.canvas.create_image(x + w / 2, y + h / 2, image=image, **options)
        self.image_map[result_id] = image
    
    def draw_line(self, x1, y1, x2, y2, color="", **options):
        '''
        Draws a line between the two given points (x1, y1) and (x2, y2).

        @param x1: starting point's x-coordinate
        @param y1: starting point's y-coordinate
        @param x2: ending point's x-coordinate
        @param y2: ending point's y-coordinate
        @param color: (optional) color with which to draw line (if absent, uses panel's default color as set via set_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_line function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_line-method
        '''
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        x1, y1, x2, y2 = DrawingPanel.__adjust_coordinates(x1, y1, x2, y2)
        color = self.__pick_color_with_default(color, self.outline, "fill", options)
        self.canvas.create_line(x1, y1, x2, y2, **options)
    
    def draw_oval(self, x, y, w, h, color="", **options):
        '''
        Draws an oval the given size with its top-left corner at the given (x, y) pixel.
        You can pass the outline color as a last positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default outline color.
        The polygon will not be filled unless a 'fill' named parameter is passed.
        
        @param x: left x coordinate of bounding box of oval
        @param y: top y coordinate of bounding box of oval
        @param w: width of bounding box of oval
        @param h: height of bounding box of oval
        @param color: (optional) color to use for outline of shape (if absent, uses panel's default color as set via set_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_oval function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_oval-method
        '''
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        color = self.__pick_color_with_default(color, self.outline, "outline", options)
        Color.normalize_all_colors(options)
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        self.canvas.create_oval(x, y, x + w, y + h, **options)

    def draw_pixel(self, x, y, color="", **options):
        '''
        Changes the color at the given (x, y) pixel.
        Equivalent to drawing a 1x1 rectangle.
        
        @param x: x coordinate
        @param y: y coordinate
        @param color: (optional) color to use for pixel (if absent, uses panel's default color as set via set_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_rectangle function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_rectangle-method
        '''
        if not "width" in options:
            options["width"] = 1
        color = self.__pick_color_with_default(color, self.outline, "outline", options)
        Color.normalize_all_colors(options)
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        self.canvas.create_rectangle(x, y, x, y, **options)

    def draw_polyline(self, *coords, **options):
        '''
        Draws a group of lines between the given coordinates passed as individual x/y coordinate parameters or (x,y) tuples.
        
        You can pass the outline color as a last positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default outline color.
        The difference between draw_polyline and draw_polygon is that the polygon will connect its last line point
        back to the start point; the polygon can also be filled, while a polyline cannot.
        
        example: panel.draw_polyline(x1, y1, x2, y2, x3, y3)
        example: panel.draw_polyline(x1, y1, x2, y2, x3, y3, x4, y4, "red")
        example: panel.draw_polyline(p1, p2, p3)
        
        @param coords: x, y coordinate pairs of arbitrary number
        @param options: (optional) keyword arguments to pass directly to Tk create_line function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_line-method
        '''
        if len(coords) == 0:
            return
        if isinstance(coords[0], list):   # allow passing list or tuple
            coords = coords[0]
        if isinstance(coords[0], tuple):
            # unpack tuples: [(10, 20), (30, 40), (50, 60)] => [10, 20, 30, 40, 50, 60]
            coords = [xy for t in coords for xy in t]
        # if arg 'is' a valid color use as color and remove from coords
        if len(coords) % 2 != 0:
            if Color.is_valid_color(coords[-1]):
                options["fill"] = Color.to_hex(coords[-1])
                coords = coords[:-1]
        if len(coords) % 2 != 0:
            raise ValueError("must pass even number of coordinates")
        adjusted_coords = DrawingPanel.__adjust_coordinates(coords)
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)

        color = self.__pick_color_with_default(None, self.outline, "fill", options)
        for i in range(0, len(coords) - 2, 2):
            x1 = coords[i]
            y1 = coords[i + 1]
            x2 = coords[i + 2]
            y2 = coords[i + 3]
            self.draw_line(x1, y1, x2, y2, color, **options)
    
    def draw_polygon(self, *coords, **options):
        '''
        Draws an outlined polygon between the given coordinates passed as individual x/y coordinate parameters or (x,y) tuples.
        You can pass the outline color as a last positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default outline color.
        The polygon will not be filled unless a 'fill' named parameter is passed.

        example: panel.draw_polygon(x1, y1, x2, y2, x3, y3)
        example: panel.draw_polygon(x1, y1, x2, y2, x3, y3, x4, y4, "red")
        example: panel.draw_polygon(p1, p2, p3)

        @param coords: x, y coordinate pairs of arbitrary number
        @param options: (optional) keyword arguments to pass directly to Tk create_polygon and/or create_line function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_polygon-method
        '''
        if len(coords) == 0:
            return
        if isinstance(coords[0], list):   # allow passing list or tuple
            coords = coords[0]
        if isinstance(coords[0], tuple):
            # unpack tuples: [(10, 20), (30, 40), (50, 60)] => [10, 20, 30, 40, 50, 60]
            coords = [xy for t in coords for xy in t]
        # if arg 'is' a valid color use as color and remove from coords
        outline_color = self.outline
        if len(coords) % 2 != 0:
            if Color.is_valid_color(coords[-1]):
                outline_color = Color.to_hex(coords[-1])
                coords = coords[:-1]
        if len(coords) % 2 != 0:
            raise ValueError("must pass even number of coordinates")
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        if not "fill" in options:
            # no fill; draw as poly line to avoid unwanted polygon fill;
            # append start point as end point to close the polygon
            if coords[0] != coords[-2] or coords[1] != coords[-1]:
                coords += (coords[0], coords[1])
            options["fill"] = outline_color
            self.draw_polyline(coords, **options)
            return

        adjusted_coords = DrawingPanel.__adjust_coordinates(coords)

        self.__pick_color_with_default(None, outline_color, "outline", options)
        Color.normalize_all_colors(options)
        self.__pick_color_with_default(None, self.get_background(), "fill", options)   # default transparent fill
        self.canvas.create_polygon(adjusted_coords, **options)

    def draw_rect(self, x, y, w, h, color="", **options):
        '''
        Draws a rectangle of the given size with its top-left corner at the given (x, y) pixel.
        Equivalent to draw_rectangle.
        
        @param x: left x coordinate of bounding box of rectangle
        @param y: top y coordinate of bounding box of rectangle
        @param w: width of bounding box of rectangle
        @param h: height of bounding box of rectangle
        @param color: (optional) color to use for outline of shape (if absent, uses panel's default color as set via set_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_rectangle function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_rectangle-method
        '''
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        color = self.__pick_color_with_default(color, self.outline, "outline", options)
        Color.normalize_all_colors(options)
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        self.canvas.create_rectangle(x, y, x + w, y + h, **options)
    
    def draw_rectangle(self, x, y, w, h, color="", **options):
        '''
        Draws a rectangle of the given size with its top-left corner at the given (x, y) pixel.
        Equivalent to draw_rect.
        
        @param x: left x coordinate of bounding box of rectangle
        @param y: top y coordinate of bounding box of rectangle
        @param w: width of bounding box of rectangle
        @param h: height of bounding box of rectangle
        @param color: (optional) color to use for outline of shape (if absent, uses panel's default color as set via set_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_rectangle function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_rectangle-method
        '''
        self.draw_rect(x, y, w, h, color, **options)

    def draw_string(self, text, x, y, color="", font="", **options):
        '''
        Draws a text string with its top-left corner at the given (x, y) pixel.

        @param x: left x coordinate of start of text
        @param y: top y coordinate of start of text
        @param color: (optional) color to use for outline of shape (if absent, uses panel's default color as set via set_color)
        @param font: (optional) font face to use to render text, in "family size style" format, such as "Times New Roman 16 bold italic"
        @param options: (optional) keyword arguments to pass directly to Tk create_text function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_text-method
        '''
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        color = self.__pick_color_with_default(color, self.outline, "fill", options)
        if font is None or font == "":
            font = self.font
        if not "justify" in options:
            options["justify"] = LEFT
        if not "anchor" in options:
            options["anchor"] = NW
        self.canvas.create_text(x, y, text=text, font=font, **options)

    def fill_arc(self, x, y, w, h, start_angle, extent, color="", **options):
        '''
        Draws a filled arc (a partial oval).
        The resulting arc begins at 'start_angle' and extends for 'extent' degrees.
        Angles are interpreted such that 0 degrees is at the 3 o'clock position.
        A positive value indicates a counter-clockwise rotation while a negative value indicates a clockwise rotation.
        You can pass the outline color as a named 'outline' parameter, else the fill color will be used as the outline color.
        
        @param x: left x coordinate of bounding box of arc's oval
        @param y: top y coordinate of bounding box of arc's oval
        @param w: width of bounding box of arc's oval
        @param h: height of bounding box of arc's oval
        @param start_angle: starting angle in degrees, where 0 is at "3 o'clock" position
        @param extent: degrees that the arc should extend from its start, where positive numbers are counter-clockwise and negative values are clockwise
        @param color: (optional) color to use for fill of shape (if absent, uses panel's default fill color as set via set_fill_color)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_arc-method
        '''
        if not "width" in options and self.stroke != 1: options["width"] = str(self.stroke)
        color = self.__pick_color_with_default(color, self.fill, "fill", options)
        self.__pick_color_with_default(None, color, "outline", options)
        Color.normalize_all_colors(options)
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        if not "start" in options:
            options["start"] = start_angle
        if not "extent" in options:
            options["extent"] = extent
        if not "style" in options:
            options["style"] = "pieslice"
        bounding_box = (x, y, x + w, y + h)
        self.canvas.create_arc(bounding_box, **options)
    
    def fill_oval(self, x, y, w, h, color="", **options):
        '''
        Draws a filled oval of the given size with its top-left corner at the given (x, y) pixel.
        You can pass the fill color as a fifth positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default fill color.
        You can pass the outline color as a named 'outline' parameter, else the fill color will be used as the outline color.
        
        @param x: left x coordinate of bounding box of oval
        @param y: top y coordinate of bounding box of oval
        @param w: width of bounding box of oval
        @param h: height of bounding box of oval
        @param color: (optional) color to use for fill of shape (if absent, uses panel's default color as set via set_fill_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_oval function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_oval-method
        '''
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        color = self.__pick_color_with_default(color, self.fill, "fill", options)
        if not "width" in options and "outline" in options and self.stroke != 1: options["width"] = str(self.stroke)

        # if outline is present as named parameter in options dict, leave it;
        # else set equal to fill color
        self.__pick_color_with_default(None, color, "outline", options)
        
        # draw the filled oval
        self.canvas.create_oval(x, y, x + w, y + h, **options)

    def fill_polygon(self, *coords, **options):
        '''
        Draws a filled polygon between the given coordinates passed as individual x/y coordinate parameters or (x,y) tuples.
        You can pass the 'fill' as a last positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default fill color.
        You can pass the outline color as a named 'outline' parameter, else the fill color will be used as the outline color.
        
        example: panel.fill_polygon(x1, y1, x2, y2, x3, y3)
        example: panel.fill_polygon(x1, y1, x2, y2, x3, y3, x4, y4, "red")
        example: panel.fill_polygon(p1, p2, p3)

        @param coords: x, y coordinate pairs of arbitrary number
        @param options: (optional) keyword arguments to pass directly to Tk create_polygon function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_polygon-method
        '''
        if len(coords) == 0:
            return
        if isinstance(coords[0], list):   # allow passing list or tuple
            coords = coords[0]
        if isinstance(coords[0], tuple):
            # unpack tuples: [(10, 20), (30, 40), (50, 60)] => [10, 20, 30, 40, 50, 60]
            coords = [xy for t in coords for xy in t]
        if not "width" in options and "outline" in options and self.stroke != 1: options["width"] = str(self.stroke)
        
        # if arg 'is' a valid color use as color and remove from coords
        fill_color = ""
        if len(coords) % 2 != 0:
            if Color.is_valid_color(coords[-1]):
                fill_color = options["fill"] = Color.to_hex(coords[-1])
                coords = coords[:-1]
        
        if len(coords) % 2 != 0:
            raise ValueError("must pass even number of coordinates")
        adjusted_coords = DrawingPanel.__adjust_coordinates(coords)

        self.__pick_color_with_default(None, fill_color, "fill", options)
        Color.normalize_all_colors(options)
        self.canvas.create_polygon(adjusted_coords, options)

    def fill_rect(self, x, y, w, h, color="", **options):
        '''
        Draws a filled rectangle of the given size with its top-left corner at the given (x, y) pixel.
        You can pass the fill color as a fifth positional parameter, or it can be a keyword parameter (in options),
        or if none is passed at all, we will fall back to the panel's default fill color.
        You can pass the outline color as a named 'outline' parameter, else the fill color will be used as the outline color.
        Equivalent to fill_rectangle.
        
        @param x: left x coordinate of bounding box of rectangle
        @param y: top y coordinate of bounding box of rectangle
        @param w: width of bounding box of rectangle
        @param h: height of bounding box of rectangle
        @param color: (optional) color to use to fill shape (if absent, uses panel's default color as set via set_fill_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_rectangle function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_rectangle-method
        '''
        x, y = DrawingPanel.__adjust_coordinates(x, y)
        color = self.__pick_color_with_default(color, self.fill, "fill", options)
        if not "width" in options and "outline" in options and self.stroke != 1: options["width"] = str(self.stroke)

        # if outline is present as named parameter in options dict, leave it;
        # else set equal to fill color
        self.__pick_color_with_default(None, color, "outline", options)
        
        # draw the filled rectangle
        self.canvas.create_rectangle(x, y, x + w, y + h, **options)

    def fill_rectangle(self, x, y, w, h, fill="", color="", **options):
        '''
        Draws a filled rectangle of the given size with its top-left corner at the given (x, y) pixel.
        Equivalent to fill_rect.
        
        @param x: left x coordinate of bounding box of rectangle
        @param y: top y coordinate of bounding box of rectangle
        @param w: width of bounding box of rectangle
        @param h: height of bounding box of rectangle
        @param color: (optional) color to use to fill shape (if absent, uses panel's default color as set via set_fill_color)
        @param options: (optional) keyword arguments to pass directly to Tk create_rectangle function (see Tkinter documentation)
        @see: http://effbot.org/tkinterbook/canvas.htm#Tkinter.Canvas.create_rectangle-method
        '''
        self.fill_rect(x, y, w, h, fill, color, **options)

    def get_background(self):
        '''
        Returns the panel's background color. Default is white.
        @return: color as hex string such as "#ff00cc"
        @rtype: str
        '''
        return self.canvas["bg"]

    def get_background_rgb(self):
        '''
        Returns the panel's background color as an RGB tuple. Default is white.
        @return: color as a 3-element rgb tuple such as (255, 0, 192)
        @rtype: tuple
        '''
        return Color.to_rgb(self.canvas["bg"])

    def get_color(self):
        '''
        Returns the default color used to outline drawn shapes.
        This color will be used unless you explicitly pass a color to the drawing function to draw a particular shape.
        Equivalent to get_foreground.
        @rtype: str
        '''
        return self.outline

    def get_fill_color(self):
        '''
        Returns the default color used to fill in shapes.
        This color will be used unless you explicitly pass a color to the filling function to draw a particular shape.
        @rtype: str
        '''
        return self.fill

    def get_font(self):
        '''
        Returns the font used for drawn strings.
        Pass a string in "FONTNAME SIZE STYLE" format, with the last two tokens optional,
        such as "Helvetica" or "Consolas 14" or "Times New Roman 16 bold".
        @return: font as a "family size style" string such as "Times New Roman 16 bold italic"
        @rtype: str
        @see: http://effbot.org/tkinterbook/tkinter-widget-styling.htm
        '''
        return self.font

    def get_foreground(self):
        '''
        Returns the default color used to outline drawn shapes.
        This color will be used unless you explicitly pass a color to the drawing function to draw a particular shape.
        Equivalent to get_color.
        @rtype: str
        '''
        return self.outline

    def get_height(self):
        '''
        Returns the height of the panel's canvas area.
        @return: height in pixels as an int
        @rtype: int
        '''
        return self.height

    def get_location(self):
        '''
        Returns a 2-element tuple containing the panel's window's top-left x/y corner.
        @rtype: tuple
        '''
        return (self.x, self.y)

    def get_outline_color(self):
        '''
        Returns the default color used to outline drawn shapes.
        This color will be used unless you explicitly pass a color to the drawing function to draw a particular shape.
        Equivalent to get_color.
        @rtype: str
        '''
        return self.outline

    def get_pixel_color(self, x, y):
        '''
        Returns the RGB color of the given (x, y) pixel on the canvas as a string such as "#ff00ff".
        NOTE: This function currently fails when there are text strings drawn onto the canvas.
        It interprets every pixel in the bounding box of the text string to be the text string's color.
        
        @param x: x coordinate of pixel
        @param y: y coordinate of pixel
        @rtype: str
        '''
        x = DrawingPanel.__adjust_coordinate(x)
        y = DrawingPanel.__adjust_coordinate(y)
        ids = self.canvas.find_overlapping(x, y, x, y)
        result = None
        if len(ids) > 0:
            index = ids[-1]
            if index in self.image_map:
                # if PhotoImage, get(x, y) to get RGB
                img = self.image_map[index]
                coords = self.canvas.coords(index)
                x = int(x - int(coords[0]) + img.width() / 2)
                y = int(y - int(coords[1]) + img.height() / 2)
                if x >= 0 and x < img.width() and y >= 0 and y < img.height():
                    rgb = img.get(x, y)
                    result = Color.to_hex(rgb[0], rgb[1], rgb[2])
            else:
                # get pixel from an ordinary shape
                color = self.canvas.itemcget(index, "fill")
                if color == "":
                    color = self.canvas.itemcget(index, "outline")
                if color == "":
                    color = self.canvas.itemcget(index, "color")
                if color != "":
                    result = color
        if result is None:
            result = self.canvas["bg"]
        temp = self.canvas.winfo_rgb(result)   # returns a 3-tuple of ints from 0-65535 representing R,G,B
        if temp:
            r = int(temp[0] / 256)
            g = int(temp[1] / 256)
            b = int(temp[2] / 256)
            result = Color.to_hex(r, g, b)   # convert to 0-255
        return result

    def get_pixel_color_rgb(self, x, y):
        '''
        Returns the RGB color of the given (x, y) pixel on the canvas as a 3-tuple of three integers from 0-255
        representing the Red, Green, and Blue components.
        
        @param x: x coordinate of pixel
        @param y: y coordinate of pixel
        @rtype: str
        '''
        return Color.to_rgb(self.get_pixel_color(x, y))

    def get_size(self):
        '''
        Returns a 2-element tuple of the width and height of the panel's canvas area.
        @return: 2-element tuple such as (400, 300)
        @rtype: tuple
        '''
        return (self.width, self.height)

    def get_stroke(self):
        '''
        Returns the width that will be used to draw lines on future shapes.
        Equivalent to get_stroke_width.
        @return: stroke width (default 1)
        @rtype: int
        '''
        return self.stroke

    def get_stroke_width(self):
        '''
        Returns the width that will be used to draw lines on future shapes.
        Equivalent to get_stroke.
        @return: stroke width (default 1)
        @rtype: int
        '''
        return self.stroke

    def get_width(self):
        '''
        Returns the width of the panel's canvas area.
        @rtype: int
        '''
        return self.width

    def get_window(self):
        '''
        Returns the underlying Tkinter window being displayed on the screen.
        @rtype: Tk
        '''
        return self.window

    def get_window_size(self):
        '''
        Returns a 2-element tuple of the width and height of the panel's entire window including borders and status bar.
        @return: 2-element tuple such as (410, 318)
        @rtype: tuple
        '''
        return (self.window.winfo_width(), self.window.winfo_height())

    def get_x(self):
        '''
        Returns the panel's window's left x pixel position.
        @rtype: int
        '''
        return self.x

    def get_y(self):
        '''
        Returns the panel's window's top y pixel position.
        @rtype: int
        '''
        return self.y

    def hide(self):
        '''
        Sets the DrawingPanel window to be not visible on the screen.
        '''
        self.set_visible(False)

    def __install_mainloop_hack(self):
        '''
        Private function that sets up a 'main loop' hack.
        Not to be called by students.
        '''
        # for everything but idle
        atexit.register(self.window.mainloop)

        # hack just for idle:
        # flush_stdout is called immediately after code execution - intercept
        # this call, and use it to call mainloop
        try:
            import idlelib.run
            def mainloop_wrap(orig_func):
                def newfunc(*a, **kw):
                    self.window.mainloop()
                    idlelib.run.flush_stdout = orig_func
                    return orig_func(*a, **kw)
                return newfunc
            idlelib.run.flush_stdout = mainloop_wrap(idlelib.run.flush_stdout)
        except ImportError:
            pass

    def is_resizable(self):
        '''
        Returns whether the DrawingPanel window is able to be resized; initially True.
        @return: True if the DrawingPanel window is able to be resized, otherwise False
        @rtype: bool
        '''
        return self.resizable

    def is_visible(self):
        '''
        Returns whether the DrawingPanel window is visible on the screen.
        If the DrawingPanel has been closed or hidden, returns False.
        @return: True if the DrawingPanel window is visible on the screen, otherwise False
        @rtype: bool
        '''
        try:
            return "normal" == self.window.state()
        except TclError:
            return False

    def message(self, text, title="Message"):
        '''
        Pops up a message box displaying the given text.
        '''
        if (sys.version_info >= (3, 0)):
            messagebox.showinfo(title, text)
        else:
            tkMessageBox.showinfo(title, text)

    def message_error(self, text, title="Error"):
        '''
        Pops up an error message box displaying the given text.
        '''
        if (sys.version_info >= (3, 0)):
            messagebox.showerror(title, text)
        else:
            tkMessageBox.showerror(title, text)

    def message_yes_no(self, text, title="Question"):
        '''
        Pops up a question message box displaying the given text with Yes/No buttons.
        If user presses Yes, returns True; if No, returns False.
        '''
        if (sys.version_info >= (3, 0)):
            return messagebox.askyesno(title, text)
        else:
            return tkMessageBox.askyesno(title, text)

    def __pick_color_with_default(self, color, default, key, options):
        '''
        Private helper function to help us allow for explicitly passed or implicit default colors.
        Chooses a color based on the given value, default, and what option key to store it as.
        If the given color is None/empty and default is present, uses that default.
        If a non-empty color or default is found, places this into options[key] for each key in keys.
        Not to be called by students.
        '''
        if color is None or color == "":
            found = False
            if key in options:
                color = options[key]
                found = True
            if not found:
                color = default
        if not (color is None or color == ""):
            color = Color.to_hex(color)
            options[key] = color
        return color

    def save(self, filename):
        '''
        Saves the drawing panel to the given output file.
        Formats supported: .ppm, .ps
        NOTE: This function currently fails when there are text strings drawn onto the canvas.
        It interprets every pixel in the bounding box of the text string to be the text string's color.
        
        @param filename: path to file to save into (must be .ps or .ppm format)
        '''
        # TODO: use PIL to save if lib is present
        
        if filename.lower().endswith(".ps"):
            self.canvas.postscript(file=filename, colormode="color")
        elif filename.lower().endswith(".ppm"):
            # save as PPM file format
            # https://en.wikipedia.org/wiki/Netpbm_format
            f = open(filename, "w")
            w = self.get_width()
            h = self.get_height()
            f.write("P3\n")                         # ppm header
            f.write(str(w) + " " + str(h) + "\n")   # width and height
            f.write("255\n")                        # max rgb value
            for y in range(h):
                # if y % 10 == 0:
                #     print("  - saved y=" + str(y))
                for x in range(w):
                    px = self.get_pixel_color(x, y)
                    rgb = Color.to_rgb(px)
                    f.write("%3d %3d %3d " % rgb)
                f.write("\n")
            f.close()
        # elif filename.lower().endswith(".bmp"):
            # TODO: save as BMP file format
            # https://en.wikipedia.org/wiki/BMP_file_format
            # f = open(filename, "wb")
            # pass
        else:
            raise ValueError("Unable to save file in this format. Supported formats: *.ppm, *.ps")

    def save_as_dialog(self):
        '''
        Pops up a file choosing dialog to allow the user to save the drawing panel's
        canvas output to a file.
        '''
        # dialog settings
        title = "Save file"
        file_types = (("PPM files", "*.ppm"), \
            ("PS files", "*.ps"), \
            ("PNG files", "*.png"), \
            ("JPG files", "*.jpg,*.jpeg"), \
            ("GIF files", "*.gif"), \
            ("All files", "*.*"))
        current_dir = os.getcwd()
        
        try:
            if (sys.version_info >= (3, 0)):
                file_path = filedialog.asksaveasfilename(initialdir = current_dir, title = title, filetypes = file_types)
            else:
                file_path = tkFileDialog.asksaveasfilename(initialdir = current_dir, title = title, filetypes = file_types)
        except TclError:
            return
        if file_path is None or file_path == "" or not isinstance(file_path, str):
            return
        
        filename = file_path
        slash = file_path.rfind("/")
        if slash < 0:
            slash = file_path.rfind("\\")
        if slash >= 0:
            filename = file_path[slash + 1 : ]

        # if file already exists, ask to overwrite
        # if os.path.isfile(file_path) and not self.message_yes_no(filename + " already exists. Overwrite?", "Overwrite file?"):
        #     return
        
        try:
            self.save(file_path)
            self.set_status("Saved to " + filename + ".")
        except ValueError as ve:
            self.message_error(str(ve), "Error saving file")

    def set_background(self, *args):
        '''
        Sets the panel's background color to be the given color.
        color -- the color to set, as a string such as "yellow" or "black"
        @param args: either a color string such as "yellow" or "#ff00cc", or a 3-element RGB tuple such as (255, 0, 192)
        '''
        self.canvas["bg"] = Color.to_hex(*args)

    def set_color(self, *args):
        '''
        Sets the color used to outline future drawn shapes.
        The color passed can be a color string like "purple" or a hex string like "#ff00ff" or an RGB tuple like (255, 0, 255).
        @param args: either a color string such as "yellow" or "#ff00cc", or a 3-element RGB tuple such as (255, 0, 192)
        @see: http://wiki.tcl.tk/16166
        '''
        color = Color.to_hex(*args)
        self.outline = color
    
    def set_fill_color(self, *args):
        '''
        Sets the color used to fill in future drawn shapes.
        The color passed can be a color string like "purple" or a hex string like "#ff00ff" or an RGB tuple like (255, 0, 255).
        @param args: either a color string such as "yellow" or "#ff00cc", or a 3-element RGB tuple such as (255, 0, 192)
        '''
        color = Color.to_hex(*args)
        self.fill = color
    
    def set_font(self, font):
        '''
        Sets the font used for future drawn strings.
        Pass a string in "FONTNAME SIZE STYLE" format, with the last two tokens optional,
        such as "Helvetica" or "Consolas 14" or "Times New Roman 16 bold".
        @param font: a font string in "FONTNAME SIZE STYLE" format, with the last two tokens optional, such as "Helvetica" or "Consolas 14" or "Times New Roman 16 bold".
        @see: http://effbot.org/tkinterbook/tkinter-widget-styling.htm
        '''
        self.font = font
    
    def set_foreground(self, *args):
        '''
        Sets the color used to outline and fill future drawn shapes.
        The color passed can be a color string like "purple" or a hex string like "#ff00ff" or an RGB tuple like (255, 0, 255).
        Equivalent to set_color.
        @param args: either a color string such as "yellow" or "#ff00cc", or a 3-element RGB tuple such as (255, 0, 192)
        '''
        self.set_color(*args)
    
    def set_height(self, height):
        '''
        Sets the panel's canvas height to the given value.
        @param height: desired height of window in pixels
        '''
        self.set_size(self.width, height)

    def set_location(self, x, y):
        '''
        Sets the window's location such that its top-left corner is at the given x/y pixel position.
        @param x: desired x coordinate of top-left corner of window
        @param y: desired y coordinate of top-left corner of window
        '''
        ww, wh = self.get_window_size()
        self.x = x
        self.y = y
        self.window.geometry("%dx%d+%d+%d" % (ww, wh, x, y))

    def set_outline_color(self, *args):
        '''
        Sets the color used to outline future drawn shapes.
        Equivalent to set_color.
        @param args: either a color string such as "yellow" or "#ff00cc", or a 3-element RGB tuple such as (255, 0, 192)
        '''
        color = Color.to_hex(*args)
        self.outline = color

    def set_resizable(self, resizable):
        '''
        Sets whether the DrawingPanel window is able to be resized; initially True.
        @param resizable: True if window should be resizable; False if not
        '''
        self.resizable = resizable
        self.window.resizable(resizable, resizable)
    
    def set_size(self, width, height):
        '''
        Sets the window's size such that its canvas will be the given width and height.
        @param width: desired canvas width in pixels
        @param height: desired canvas height in pixels
        '''
        # TODO: fix
        self.width = width
        self.height = height
        self.window.geometry("%dx%d+%d+%d" % (self.width, self.height, self.x, self.y))
        if hasattr(self, "canvas"):
            self.canvas.config(width=width, height=height)
            self.canvas.pack(side = TOP)
        if hasattr(self, "status_label"):
            self.status_label.pack(side = LEFT, expand = True)
        self.window.update()

    def set_status(self, text):
        '''
        Sets the text to show in the window's southern status bar area.
        '''
        self.status_label.config(text=text, anchor=W, justify=LEFT)
        self.window.update()

    def set_stroke(self, width):
        '''
        Sets the width that will be used to draw lines on future shapes.
        Equivalent to set_stroke_width.
        '''
        self.stroke = width

    def set_stroke_width(self, width):
        '''
        Sets the width that will be used to draw lines on future shapes.
        Equivalent to set_stroke.
        '''
        self.stroke = width

    def set_title(self, text):
        '''
        Sets the title to display at the top of the window.
        '''
        self.window.title(text)

    def set_visible(self, visible=True):
        '''
        Sets whether the DrawingPanel window should be visible on the screen.
        '''
        try:
            if visible and not self.is_visible():
                self.window.update()
                self.window.deiconify()
            elif not visible and self.is_visible():
                self.window.withdraw()
        except TclError:
            pass

    def set_width(self, width):
        '''
        Sets the panel's canvas width to the given value.
        '''
        self.set_size(width, self.height)

    def set_x(self, x):
        '''
        Sets the panel's window's left x pixel position to the given value.
        @param x: desired x position of top-left corner of window
        '''
        self.set_location(x, self.y)

    def set_y(self, y):
        '''
        Sets the panel's window's top y pixel position to the given value.
        @param y: desired y position of top-left corner of window
        '''
        self.set_location(self.x, y)

    def show(self):
        '''
        Sets the DrawingPanel window to be visible on the screen.
        '''
        self.set_visible(True)

    def sleep(self, ms):
        '''
        Causes the DrawingPanel to pause for the given number of milliseconds.
        Useful for creating simple animations.
        @param ms: number of milliseconds to pause
        '''
        if ms <= 0:
            return
        try:
            self.window.update()
            time.sleep(ms / 1000.0)   # time.sleep takes a number of seconds (can be float)
            self.window.update()
        except Exception:
            pass   # swallow exception

class Color:
    '''
    A class to deal with colors.
    A color can be specified in four ways by the client:
    1) as a Color constant like Color.RED or Color.LIGHT_SKY_BLUE
    2) as a name like "red" or "light sky blue"
    3) as a six-character hex string like "#ff00cc"
    4) as a three-integer tuple like (255, 0, 198) where each int is in range 0-255
    '''
    def __init__(self, *args):
        '''
        Constructs a color.
        '''
        self.hex_string = Color.to_name(*args)
    
    def __str__(self):
        '''
        Returns a string representation of this color.
        '''
        return self.hex_string
    
    @staticmethod
    def is_valid_color(color):
        '''
        Returns True if the color passed is a valid color in any form, such as
        a Color object, a hex string such as "#ff00cc", a color name such as "red" or "light sky blue",
        or RGB color expressed as a 3-element tuple or list.
        '''
        return color is Color or Color.is_hex(color) or Color.is_name(color) or Color.is_rgb(color)

    @staticmethod
    def is_hex(color):
        '''
        Returns True if the color passed is a valid hex color string such as "#ff00cc".
        '''
        if isinstance(color, str):
            color = color.lower()
            if color[0] != "#":
                return False
            for i in range(1, len(color)):
                if not color[i] in "0123456789abcdef":
                    return False
            return True
        else:
            return False

    @staticmethod
    def is_name(color):
        '''
        Returns True if the color passed is a valid color name such as "red" or "light sky blue".
        '''
        return isinstance(color, str) and color.lower() in Color.NAME_TO_RGB_LC

    @staticmethod
    def is_rgb(color):
        '''
        Returns True if the color passed is a valid RGB color expressed as a 3-element tuple or list.
        '''
        if color is None:
            return False
        elif isinstance(color, tuple) or isinstance(color, list):
            if len(color) != 3:
                return False
            return 0 <= color[0] <= 255 and 0 <= color[1] <= 255 and 0 <= color[2] <= 255
        else:
            return False

    @staticmethod
    def to_hex(*args):
        '''
        Converts a color into a hex color string such as "#ff00ff".
        Accepts a Color in any format, such as a Color object, a color name, a hex string, or tuple of three RGB integers from 0-255.
        If an invalid color is passed, raises a ValueError.
        '''
        if len(args) == 0:
            raise ValueError
        if args[0] is None:
            return None
        elif isinstance(args[0], Color):
            # Color objects are always stored/printed as hex already
            return str(args[0])
        elif len(args) == 3:
            if Color.is_rgb(args):
                # RGB tuple
                r, g, b = args
                return "#%02x%02x%02x" % (int(r), int(g), int(b))
        elif isinstance(args[0], str):
            if Color.is_hex(args[0]):
                return args[0]   # already a hex color
            elif args[0].lower() in Color.NAME_TO_HEX_LC:
                # check for color name like "red"
                return Color.NAME_TO_HEX_LC[args[0].lower()]
            else:
                raise ValueError("invalid color: " + str(args[0]))
        elif isinstance(args[0], tuple) and len(args[0]) == 3:
            return "#%02x%02x%02x" % args[0]
        raise ValueError("invalid color: " + str(args))
    
    @staticmethod
    def to_name(*args):
        '''
        Converts a color into a color name string such as "red" or "light sky blue".
        Accepts a Color in any format, such as a Color object, a color name, a hex string, or tuple of three RGB integers from 0-255.
        If an invalid color is passed, or one that does not correspond to a color name, returns the hex string version of the color.
        '''
        hex_string = Color.to_hex(*args)
        if hex_string in Color.HEX_TO_NAME:
            return Color.HEX_TO_NAME[hex_string]
        else:
            return hex_string
    
    @staticmethod
    def to_rgb(*args):
        '''
        Converts a color into an RGB tuple such as (255, 0, 128).
        Accepts a Color in any format, such as a Color object, a color name, a hex string, or tuple of three RGB integers from 0-255.
        If an invalid color is passed, or one that does not correspond to a color name, raises a ValueError.
        '''
        hex_string = Color.to_hex(*args)
        hex_string = str(hex_string).replace("#", "")
        r = int(hex_string[0:2], 16)
        g = int(hex_string[2:4], 16)
        b = int(hex_string[4:6], 16)
        return (r, g, b)
    
    @staticmethod
    def normalize_all_colors(options):
        '''
        Makes sure that any color values in the given dict are in normalized hex string form,
        not tuples or color names.
        Not to be called by students.
        '''
        keys = ["color", "fill", "outline"]
        for key in keys:
            if key in options:
                options[key] = Color.to_hex(options[key])

    # begin looong constants for colors and names
    ALICE_BLUE = "alice blue"
    ALICEBLUE = "AliceBlue"
    ANTIQUE_WHITE = "antique white"
    ANTIQUEWHITE = "AntiqueWhite"
    ANTIQUEWHITE1 = "AntiqueWhite1"
    ANTIQUEWHITE2 = "AntiqueWhite2"
    ANTIQUEWHITE3 = "AntiqueWhite3"
    ANTIQUEWHITE4 = "AntiqueWhite4"
    AQUAMARINE = "aquamarine"
    AQUAMARINE1 = "aquamarine1"
    AQUAMARINE2 = "aquamarine2"
    AQUAMARINE3 = "aquamarine3"
    AQUAMARINE4 = "aquamarine4"
    AZURE = "azure"
    AZURE1 = "azure1"
    AZURE2 = "azure2"
    AZURE3 = "azure3"
    AZURE4 = "azure4"
    BEIGE = "beige"
    BISQUE = "bisque"
    BISQUE1 = "bisque1"
    BISQUE2 = "bisque2"
    BISQUE3 = "bisque3"
    BISQUE4 = "bisque4"
    BLACK = "black"
    BLANCHED_ALMOND = "blanched almond"
    BLANCHEDALMOND = "BlanchedAlmond"
    BLUE = "blue"
    BLUE_VIOLET = "blue violet"
    BLUE1 = "blue1"
    BLUE2 = "blue2"
    BLUE3 = "blue3"
    BLUE4 = "blue4"
    BLUEVIOLET = "BlueViolet"
    BROWN = "brown"
    BROWN1 = "brown1"
    BROWN2 = "brown2"
    BROWN3 = "brown3"
    BROWN4 = "brown4"
    BURLYWOOD = "burlywood"
    BURLYWOOD1 = "burlywood1"
    BURLYWOOD2 = "burlywood2"
    BURLYWOOD3 = "burlywood3"
    BURLYWOOD4 = "burlywood4"
    CADET_BLUE = "cadet blue"
    CADETBLUE = "CadetBlue"
    CADETBLUE1 = "CadetBlue1"
    CADETBLUE2 = "CadetBlue2"
    CADETBLUE3 = "CadetBlue3"
    CADETBLUE4 = "CadetBlue4"
    CHARTREUSE = "chartreuse"
    CHARTREUSE1 = "chartreuse1"
    CHARTREUSE2 = "chartreuse2"
    CHARTREUSE3 = "chartreuse3"
    CHARTREUSE4 = "chartreuse4"
    CHOCOLATE = "chocolate"
    CHOCOLATE1 = "chocolate1"
    CHOCOLATE2 = "chocolate2"
    CHOCOLATE3 = "chocolate3"
    CHOCOLATE4 = "chocolate4"
    CORAL = "coral"
    CORAL1 = "coral1"
    CORAL2 = "coral2"
    CORAL3 = "coral3"
    CORAL4 = "coral4"
    CORNFLOWER_BLUE = "cornflower blue"
    CORNFLOWERBLUE = "CornflowerBlue"
    CORNSILK = "cornsilk"
    CORNSILK1 = "cornsilk1"
    CORNSILK2 = "cornsilk2"
    CORNSILK3 = "cornsilk3"
    CORNSILK4 = "cornsilk4"
    CYAN = "cyan"
    CYAN1 = "cyan1"
    CYAN2 = "cyan2"
    CYAN3 = "cyan3"
    CYAN4 = "cyan4"
    DARK_BLUE = "dark blue"
    DARK_CYAN = "dark cyan"
    DARK_GOLDENROD = "dark goldenrod"
    DARK_GRAY = "dark gray"
    DARK_GREEN = "dark green"
    DARK_GREY = "dark grey"
    DARK_KHAKI = "dark khaki"
    DARK_MAGENTA = "dark magenta"
    DARK_OLIVE_GREEN = "dark olive green"
    DARK_ORANGE = "dark orange"
    DARK_ORCHID = "dark orchid"
    DARK_RED = "dark red"
    DARK_SALMON = "dark salmon"
    DARK_SEA_GREEN = "dark sea green"
    DARK_SLATE_BLUE = "dark slate blue"
    DARK_SLATE_GRAY = "dark slate gray"
    DARK_SLATE_GREY = "dark slate grey"
    DARK_TURQUOISE = "dark turquoise"
    DARK_VIOLET = "dark violet"
    DARKBLUE = "DarkBlue"
    DARKCYAN = "DarkCyan"
    DARKGOLDENROD = "DarkGoldenrod"
    DARKGOLDENROD1 = "DarkGoldenrod1"
    DARKGOLDENROD2 = "DarkGoldenrod2"
    DARKGOLDENROD3 = "DarkGoldenrod3"
    DARKGOLDENROD4 = "DarkGoldenrod4"
    DARKGRAY = "DarkGray"
    DARKGREEN = "DarkGreen"
    DARKGREY = "DarkGrey"
    DARKKHAKI = "DarkKhaki"
    DARKMAGENTA = "DarkMagenta"
    DARKOLIVEGREEN = "DarkOliveGreen"
    DARKOLIVEGREEN1 = "DarkOliveGreen1"
    DARKOLIVEGREEN2 = "DarkOliveGreen2"
    DARKOLIVEGREEN3 = "DarkOliveGreen3"
    DARKOLIVEGREEN4 = "DarkOliveGreen4"
    DARKORANGE = "DarkOrange"
    DARKORANGE1 = "DarkOrange1"
    DARKORANGE2 = "DarkOrange2"
    DARKORANGE3 = "DarkOrange3"
    DARKORANGE4 = "DarkOrange4"
    DARKORCHID = "DarkOrchid"
    DARKORCHID1 = "DarkOrchid1"
    DARKORCHID2 = "DarkOrchid2"
    DARKORCHID3 = "DarkOrchid3"
    DARKORCHID4 = "DarkOrchid4"
    DARKRED = "DarkRed"
    DARKSALMON = "DarkSalmon"
    DARKSEAGREEN = "DarkSeaGreen"
    DARKSEAGREEN1 = "DarkSeaGreen1"
    DARKSEAGREEN2 = "DarkSeaGreen2"
    DARKSEAGREEN3 = "DarkSeaGreen3"
    DARKSEAGREEN4 = "DarkSeaGreen4"
    DARKSLATEBLUE = "DarkSlateBlue"
    DARKSLATEGRAY = "DarkSlateGray"
    DARKSLATEGRAY1 = "DarkSlateGray1"
    DARKSLATEGRAY2 = "DarkSlateGray2"
    DARKSLATEGRAY3 = "DarkSlateGray3"
    DARKSLATEGRAY4 = "DarkSlateGray4"
    DARKSLATEGREY = "DarkSlateGrey"
    DARKTURQUOISE = "DarkTurquoise"
    DARKVIOLET = "DarkViolet"
    DEEP_PINK = "deep pink"
    DEEP_SKY_BLUE = "deep sky blue"
    DEEPPINK = "DeepPink"
    DEEPPINK1 = "DeepPink1"
    DEEPPINK2 = "DeepPink2"
    DEEPPINK3 = "DeepPink3"
    DEEPPINK4 = "DeepPink4"
    DEEPSKYBLUE = "DeepSkyBlue"
    DEEPSKYBLUE1 = "DeepSkyBlue1"
    DEEPSKYBLUE2 = "DeepSkyBlue2"
    DEEPSKYBLUE3 = "DeepSkyBlue3"
    DEEPSKYBLUE4 = "DeepSkyBlue4"
    DIM_GRAY = "dim gray"
    DIM_GREY = "dim grey"
    DIMGRAY = "DimGray"
    DIMGREY = "DimGrey"
    DODGER_BLUE = "dodger blue"
    DODGERBLUE = "DodgerBlue"
    DODGERBLUE1 = "DodgerBlue1"
    DODGERBLUE2 = "DodgerBlue2"
    DODGERBLUE3 = "DodgerBlue3"
    DODGERBLUE4 = "DodgerBlue4"
    FIREBRICK = "firebrick"
    FIREBRICK1 = "firebrick1"
    FIREBRICK2 = "firebrick2"
    FIREBRICK3 = "firebrick3"
    FIREBRICK4 = "firebrick4"
    FLORAL_WHITE = "floral white"
    FLORALWHITE = "FloralWhite"
    FOREST_GREEN = "forest green"
    FORESTGREEN = "ForestGreen"
    GAINSBORO = "gainsboro"
    GHOST_WHITE = "ghost white"
    GHOSTWHITE = "GhostWhite"
    GOLD = "gold"
    GOLD1 = "gold1"
    GOLD2 = "gold2"
    GOLD3 = "gold3"
    GOLD4 = "gold4"
    GOLDENROD = "goldenrod"
    GOLDENROD1 = "goldenrod1"
    GOLDENROD2 = "goldenrod2"
    GOLDENROD3 = "goldenrod3"
    GOLDENROD4 = "goldenrod4"
    GRAY = "gray"
    GRAY0 = "gray0"
    GRAY1 = "gray1"
    GRAY2 = "gray2"
    GRAY3 = "gray3"
    GRAY4 = "gray4"
    GRAY5 = "gray5"
    GRAY6 = "gray6"
    GRAY7 = "gray7"
    GRAY8 = "gray8"
    GRAY9 = "gray9"
    GRAY10 = "gray10"
    GRAY11 = "gray11"
    GRAY12 = "gray12"
    GRAY13 = "gray13"
    GRAY14 = "gray14"
    GRAY15 = "gray15"
    GRAY16 = "gray16"
    GRAY17 = "gray17"
    GRAY18 = "gray18"
    GRAY19 = "gray19"
    GRAY20 = "gray20"
    GRAY21 = "gray21"
    GRAY22 = "gray22"
    GRAY23 = "gray23"
    GRAY24 = "gray24"
    GRAY25 = "gray25"
    GRAY26 = "gray26"
    GRAY27 = "gray27"
    GRAY28 = "gray28"
    GRAY29 = "gray29"
    GRAY30 = "gray30"
    GRAY31 = "gray31"
    GRAY32 = "gray32"
    GRAY33 = "gray33"
    GRAY34 = "gray34"
    GRAY35 = "gray35"
    GRAY36 = "gray36"
    GRAY37 = "gray37"
    GRAY38 = "gray38"
    GRAY39 = "gray39"
    GRAY40 = "gray40"
    GRAY41 = "gray41"
    GRAY42 = "gray42"
    GRAY43 = "gray43"
    GRAY44 = "gray44"
    GRAY45 = "gray45"
    GRAY46 = "gray46"
    GRAY47 = "gray47"
    GRAY48 = "gray48"
    GRAY49 = "gray49"
    GRAY50 = "gray50"
    GRAY51 = "gray51"
    GRAY52 = "gray52"
    GRAY53 = "gray53"
    GRAY54 = "gray54"
    GRAY55 = "gray55"
    GRAY56 = "gray56"
    GRAY57 = "gray57"
    GRAY58 = "gray58"
    GRAY59 = "gray59"
    GRAY60 = "gray60"
    GRAY61 = "gray61"
    GRAY62 = "gray62"
    GRAY63 = "gray63"
    GRAY64 = "gray64"
    GRAY65 = "gray65"
    GRAY66 = "gray66"
    GRAY67 = "gray67"
    GRAY68 = "gray68"
    GRAY69 = "gray69"
    GRAY70 = "gray70"
    GRAY71 = "gray71"
    GRAY72 = "gray72"
    GRAY73 = "gray73"
    GRAY74 = "gray74"
    GRAY75 = "gray75"
    GRAY76 = "gray76"
    GRAY77 = "gray77"
    GRAY78 = "gray78"
    GRAY79 = "gray79"
    GRAY80 = "gray80"
    GRAY81 = "gray81"
    GRAY82 = "gray82"
    GRAY83 = "gray83"
    GRAY84 = "gray84"
    GRAY85 = "gray85"
    GRAY86 = "gray86"
    GRAY87 = "gray87"
    GRAY88 = "gray88"
    GRAY89 = "gray89"
    GRAY90 = "gray90"
    GRAY91 = "gray91"
    GRAY92 = "gray92"
    GRAY93 = "gray93"
    GRAY94 = "gray94"
    GRAY95 = "gray95"
    GRAY96 = "gray96"
    GRAY97 = "gray97"
    GRAY98 = "gray98"
    GRAY99 = "gray99"
    GRAY100 = "gray100"
    GREEN = "green"
    GREEN_YELLOW = "green yellow"
    GREEN1 = "green1"
    GREEN2 = "green2"
    GREEN3 = "green3"
    GREEN4 = "green4"
    GREENYELLOW = "GreenYellow"
    GREY = "grey"
    GREY0 = "grey0"
    GREY1 = "grey1"
    GREY2 = "grey2"
    GREY3 = "grey3"
    GREY4 = "grey4"
    GREY5 = "grey5"
    GREY6 = "grey6"
    GREY7 = "grey7"
    GREY8 = "grey8"
    GREY9 = "grey9"
    GREY10 = "grey10"
    GREY11 = "grey11"
    GREY12 = "grey12"
    GREY13 = "grey13"
    GREY14 = "grey14"
    GREY15 = "grey15"
    GREY16 = "grey16"
    GREY17 = "grey17"
    GREY18 = "grey18"
    GREY19 = "grey19"
    GREY20 = "grey20"
    GREY21 = "grey21"
    GREY22 = "grey22"
    GREY23 = "grey23"
    GREY24 = "grey24"
    GREY25 = "grey25"
    GREY26 = "grey26"
    GREY27 = "grey27"
    GREY28 = "grey28"
    GREY29 = "grey29"
    GREY30 = "grey30"
    GREY31 = "grey31"
    GREY32 = "grey32"
    GREY33 = "grey33"
    GREY34 = "grey34"
    GREY35 = "grey35"
    GREY36 = "grey36"
    GREY37 = "grey37"
    GREY38 = "grey38"
    GREY39 = "grey39"
    GREY40 = "grey40"
    GREY41 = "grey41"
    GREY42 = "grey42"
    GREY43 = "grey43"
    GREY44 = "grey44"
    GREY45 = "grey45"
    GREY46 = "grey46"
    GREY47 = "grey47"
    GREY48 = "grey48"
    GREY49 = "grey49"
    GREY50 = "grey50"
    GREY51 = "grey51"
    GREY52 = "grey52"
    GREY53 = "grey53"
    GREY54 = "grey54"
    GREY55 = "grey55"
    GREY56 = "grey56"
    GREY57 = "grey57"
    GREY58 = "grey58"
    GREY59 = "grey59"
    GREY60 = "grey60"
    GREY61 = "grey61"
    GREY62 = "grey62"
    GREY63 = "grey63"
    GREY64 = "grey64"
    GREY65 = "grey65"
    GREY66 = "grey66"
    GREY67 = "grey67"
    GREY68 = "grey68"
    GREY69 = "grey69"
    GREY70 = "grey70"
    GREY71 = "grey71"
    GREY72 = "grey72"
    GREY73 = "grey73"
    GREY74 = "grey74"
    GREY75 = "grey75"
    GREY76 = "grey76"
    GREY77 = "grey77"
    GREY78 = "grey78"
    GREY79 = "grey79"
    GREY80 = "grey80"
    GREY81 = "grey81"
    GREY82 = "grey82"
    GREY83 = "grey83"
    GREY84 = "grey84"
    GREY85 = "grey85"
    GREY86 = "grey86"
    GREY87 = "grey87"
    GREY88 = "grey88"
    GREY89 = "grey89"
    GREY90 = "grey90"
    GREY91 = "grey91"
    GREY92 = "grey92"
    GREY93 = "grey93"
    GREY94 = "grey94"
    GREY95 = "grey95"
    GREY96 = "grey96"
    GREY97 = "grey97"
    GREY98 = "grey98"
    GREY99 = "grey99"
    GREY100 = "grey100"
    HONEYDEW = "honeydew"
    HONEYDEW1 = "honeydew1"
    HONEYDEW2 = "honeydew2"
    HONEYDEW3 = "honeydew3"
    HONEYDEW4 = "honeydew4"
    HOT_PINK = "hot pink"
    HOTPINK = "HotPink"
    HOTPINK1 = "HotPink1"
    HOTPINK2 = "HotPink2"
    HOTPINK3 = "HotPink3"
    HOTPINK4 = "HotPink4"
    INDIAN_RED = "indian red"
    INDIANRED = "IndianRed"
    INDIANRED1 = "IndianRed1"
    INDIANRED2 = "IndianRed2"
    INDIANRED3 = "IndianRed3"
    INDIANRED4 = "IndianRed4"
    IVORY = "ivory"
    IVORY1 = "ivory1"
    IVORY2 = "ivory2"
    IVORY3 = "ivory3"
    IVORY4 = "ivory4"
    KHAKI = "khaki"
    KHAKI1 = "khaki1"
    KHAKI2 = "khaki2"
    KHAKI3 = "khaki3"
    KHAKI4 = "khaki4"
    LAVENDER = "lavender"
    LAVENDER_BLUSH = "lavender blush"
    LAVENDERBLUSH = "LavenderBlush"
    LAVENDERBLUSH1 = "LavenderBlush1"
    LAVENDERBLUSH2 = "LavenderBlush2"
    LAVENDERBLUSH3 = "LavenderBlush3"
    LAVENDERBLUSH4 = "LavenderBlush4"
    LAWN_GREEN = "lawn green"
    LAWNGREEN = "LawnGreen"
    LEMON_CHIFFON = "lemon chiffon"
    LEMONCHIFFON = "LemonChiffon"
    LEMONCHIFFON1 = "LemonChiffon1"
    LEMONCHIFFON2 = "LemonChiffon2"
    LEMONCHIFFON3 = "LemonChiffon3"
    LEMONCHIFFON4 = "LemonChiffon4"
    LIGHT_BLUE = "light blue"
    LIGHT_CORAL = "light coral"
    LIGHT_CYAN = "light cyan"
    LIGHT_GOLDENROD = "light goldenrod"
    LIGHT_GOLDENROD_YELLOW = "light goldenrod yellow"
    LIGHT_GRAY = "light gray"
    LIGHT_GREEN = "light green"
    LIGHT_GREY = "light grey"
    LIGHT_PINK = "light pink"
    LIGHT_SALMON = "light salmon"
    LIGHT_SEA_GREEN = "light sea green"
    LIGHT_SKY_BLUE = "light sky blue"
    LIGHT_SLATE_BLUE = "light slate blue"
    LIGHT_SLATE_GRAY = "light slate gray"
    LIGHT_SLATE_GREY = "light slate grey"
    LIGHT_STEEL_BLUE = "light steel blue"
    LIGHT_YELLOW = "light yellow"
    LIGHTBLUE = "LightBlue"
    LIGHTBLUE1 = "LightBlue1"
    LIGHTBLUE2 = "LightBlue2"
    LIGHTBLUE3 = "LightBlue3"
    LIGHTBLUE4 = "LightBlue4"
    LIGHTCORAL = "LightCoral"
    LIGHTCYAN = "LightCyan"
    LIGHTCYAN1 = "LightCyan1"
    LIGHTCYAN2 = "LightCyan2"
    LIGHTCYAN3 = "LightCyan3"
    LIGHTCYAN4 = "LightCyan4"
    LIGHTGOLDENROD = "LightGoldenrod"
    LIGHTGOLDENROD1 = "LightGoldenrod1"
    LIGHTGOLDENROD2 = "LightGoldenrod2"
    LIGHTGOLDENROD3 = "LightGoldenrod3"
    LIGHTGOLDENROD4 = "LightGoldenrod4"
    LIGHTGOLDENRODYELLOW = "LightGoldenrodYellow"
    LIGHTGRAY = "LightGray"
    LIGHTGREEN = "LightGreen"
    LIGHTGREY = "LightGrey"
    LIGHTPINK = "LightPink"
    LIGHTPINK1 = "LightPink1"
    LIGHTPINK2 = "LightPink2"
    LIGHTPINK3 = "LightPink3"
    LIGHTPINK4 = "LightPink4"
    LIGHTSALMON = "LightSalmon"
    LIGHTSALMON1 = "LightSalmon1"
    LIGHTSALMON2 = "LightSalmon2"
    LIGHTSALMON3 = "LightSalmon3"
    LIGHTSALMON4 = "LightSalmon4"
    LIGHTSEAGREEN = "LightSeaGreen"
    LIGHTSKYBLUE = "LightSkyBlue"
    LIGHTSKYBLUE1 = "LightSkyBlue1"
    LIGHTSKYBLUE2 = "LightSkyBlue2"
    LIGHTSKYBLUE3 = "LightSkyBlue3"
    LIGHTSKYBLUE4 = "LightSkyBlue4"
    LIGHTSLATEBLUE = "LightSlateBlue"
    LIGHTSLATEGRAY = "LightSlateGray"
    LIGHTSLATEGREY = "LightSlateGrey"
    LIGHTSTEELBLUE = "LightSteelBlue"
    LIGHTSTEELBLUE1 = "LightSteelBlue1"
    LIGHTSTEELBLUE2 = "LightSteelBlue2"
    LIGHTSTEELBLUE3 = "LightSteelBlue3"
    LIGHTSTEELBLUE4 = "LightSteelBlue4"
    LIGHTYELLOW = "LightYellow"
    LIGHTYELLOW1 = "LightYellow1"
    LIGHTYELLOW2 = "LightYellow2"
    LIGHTYELLOW3 = "LightYellow3"
    LIGHTYELLOW4 = "LightYellow4"
    LIME_GREEN = "lime green"
    LIMEGREEN = "LimeGreen"
    LINEN = "linen"
    MAGENTA = "magenta"
    MAGENTA1 = "magenta1"
    MAGENTA2 = "magenta2"
    MAGENTA3 = "magenta3"
    MAGENTA4 = "magenta4"
    MAROON = "maroon"
    MAROON1 = "maroon1"
    MAROON2 = "maroon2"
    MAROON3 = "maroon3"
    MAROON4 = "maroon4"
    MEDIUM_AQUAMARINE = "medium aquamarine"
    MEDIUM_BLUE = "medium blue"
    MEDIUM_ORCHID = "medium orchid"
    MEDIUM_PURPLE = "medium purple"
    MEDIUM_SEA_GREEN = "medium sea green"
    MEDIUM_SLATE_BLUE = "medium slate blue"
    MEDIUM_SPRING_GREEN = "medium spring green"
    MEDIUM_TURQUOISE = "medium turquoise"
    MEDIUM_VIOLET_RED = "medium violet red"
    MEDIUMAQUAMARINE = "MediumAquamarine"
    MEDIUMBLUE = "MediumBlue"
    MEDIUMORCHID = "MediumOrchid"
    MEDIUMORCHID1 = "MediumOrchid1"
    MEDIUMORCHID2 = "MediumOrchid2"
    MEDIUMORCHID3 = "MediumOrchid3"
    MEDIUMORCHID4 = "MediumOrchid4"
    MEDIUMPURPLE = "MediumPurple"
    MEDIUMPURPLE1 = "MediumPurple1"
    MEDIUMPURPLE2 = "MediumPurple2"
    MEDIUMPURPLE3 = "MediumPurple3"
    MEDIUMPURPLE4 = "MediumPurple4"
    MEDIUMSEAGREEN = "MediumSeaGreen"
    MEDIUMSLATEBLUE = "MediumSlateBlue"
    MEDIUMSPRINGGREEN = "MediumSpringGreen"
    MEDIUMTURQUOISE = "MediumTurquoise"
    MEDIUMVIOLETRED = "MediumVioletRed"
    MIDNIGHT_BLUE = "midnight blue"
    MIDNIGHTBLUE = "MidnightBlue"
    MINT_CREAM = "mint cream"
    MINTCREAM = "MintCream"
    MISTY_ROSE = "misty rose"
    MISTYROSE = "MistyRose"
    MISTYROSE1 = "MistyRose1"
    MISTYROSE2 = "MistyRose2"
    MISTYROSE3 = "MistyRose3"
    MISTYROSE4 = "MistyRose4"
    MOCCASIN = "moccasin"
    NAVAJO_WHITE = "navajo white"
    NAVAJOWHITE = "NavajoWhite"
    NAVAJOWHITE1 = "NavajoWhite1"
    NAVAJOWHITE2 = "NavajoWhite2"
    NAVAJOWHITE3 = "NavajoWhite3"
    NAVAJOWHITE4 = "NavajoWhite4"
    NAVY = "navy"
    NAVY_BLUE = "navy blue"
    NAVYBLUE = "NavyBlue"
    OLD_LACE = "old lace"
    OLDLACE = "OldLace"
    OLIVE_DRAB = "olive drab"
    OLIVEDRAB = "OliveDrab"
    OLIVEDRAB1 = "OliveDrab1"
    OLIVEDRAB2 = "OliveDrab2"
    OLIVEDRAB3 = "OliveDrab3"
    OLIVEDRAB4 = "OliveDrab4"
    ORANGE = "orange"
    ORANGE_RED = "orange red"
    ORANGE1 = "orange1"
    ORANGE2 = "orange2"
    ORANGE3 = "orange3"
    ORANGE4 = "orange4"
    ORANGERED = "OrangeRed"
    ORANGERED1 = "OrangeRed1"
    ORANGERED2 = "OrangeRed2"
    ORANGERED3 = "OrangeRed3"
    ORANGERED4 = "OrangeRed4"
    ORCHID = "orchid"
    ORCHID1 = "orchid1"
    ORCHID2 = "orchid2"
    ORCHID3 = "orchid3"
    ORCHID4 = "orchid4"
    PALE_GOLDENROD = "pale goldenrod"
    PALE_GREEN = "pale green"
    PALE_TURQUOISE = "pale turquoise"
    PALE_VIOLET_RED = "pale violet red"
    PALEGOLDENROD = "PaleGoldenrod"
    PALEGREEN = "PaleGreen"
    PALEGREEN1 = "PaleGreen1"
    PALEGREEN2 = "PaleGreen2"
    PALEGREEN3 = "PaleGreen3"
    PALEGREEN4 = "PaleGreen4"
    PALETURQUOISE = "PaleTurquoise"
    PALETURQUOISE1 = "PaleTurquoise1"
    PALETURQUOISE2 = "PaleTurquoise2"
    PALETURQUOISE3 = "PaleTurquoise3"
    PALETURQUOISE4 = "PaleTurquoise4"
    PALEVIOLETRED = "PaleVioletRed"
    PALEVIOLETRED1 = "PaleVioletRed1"
    PALEVIOLETRED2 = "PaleVioletRed2"
    PALEVIOLETRED3 = "PaleVioletRed3"
    PALEVIOLETRED4 = "PaleVioletRed4"
    PAPAYA_WHIP = "papaya whip"
    PAPAYAWHIP = "PapayaWhip"
    PEACH_PUFF = "peach puff"
    PEACHPUFF = "PeachPuff"
    PEACHPUFF1 = "PeachPuff1"
    PEACHPUFF2 = "PeachPuff2"
    PEACHPUFF3 = "PeachPuff3"
    PEACHPUFF4 = "PeachPuff4"
    PERU = "peru"
    PINK = "pink"
    PINK1 = "pink1"
    PINK2 = "pink2"
    PINK3 = "pink3"
    PINK4 = "pink4"
    PLUM = "plum"
    PLUM1 = "plum1"
    PLUM2 = "plum2"
    PLUM3 = "plum3"
    PLUM4 = "plum4"
    POWDER_BLUE = "powder blue"
    POWDERBLUE = "PowderBlue"
    PURPLE = "purple"
    PURPLE1 = "purple1"
    PURPLE2 = "purple2"
    PURPLE3 = "purple3"
    PURPLE4 = "purple4"
    RED = "red"
    RED1 = "red1"
    RED2 = "red2"
    RED3 = "red3"
    RED4 = "red4"
    ROSY_BROWN = "rosy brown"
    ROSYBROWN = "RosyBrown"
    ROSYBROWN1 = "RosyBrown1"
    ROSYBROWN2 = "RosyBrown2"
    ROSYBROWN3 = "RosyBrown3"
    ROSYBROWN4 = "RosyBrown4"
    ROYAL_BLUE = "royal blue"
    ROYALBLUE = "RoyalBlue"
    ROYALBLUE1 = "RoyalBlue1"
    ROYALBLUE2 = "RoyalBlue2"
    ROYALBLUE3 = "RoyalBlue3"
    ROYALBLUE4 = "RoyalBlue4"
    SADDLE_BROWN = "saddle brown"
    SADDLEBROWN = "SaddleBrown"
    SALMON = "salmon"
    SALMON1 = "salmon1"
    SALMON2 = "salmon2"
    SALMON3 = "salmon3"
    SALMON4 = "salmon4"
    SANDY_BROWN = "sandy brown"
    SANDYBROWN = "SandyBrown"
    SEA_GREEN = "sea green"
    SEAGREEN = "SeaGreen"
    SEAGREEN1 = "SeaGreen1"
    SEAGREEN2 = "SeaGreen2"
    SEAGREEN3 = "SeaGreen3"
    SEAGREEN4 = "SeaGreen4"
    SEASHELL = "seashell"
    SEASHELL1 = "seashell1"
    SEASHELL2 = "seashell2"
    SEASHELL3 = "seashell3"
    SEASHELL4 = "seashell4"
    SIENNA = "sienna"
    SIENNA1 = "sienna1"
    SIENNA2 = "sienna2"
    SIENNA3 = "sienna3"
    SIENNA4 = "sienna4"
    SKY_BLUE = "sky blue"
    SKYBLUE = "SkyBlue"
    SKYBLUE1 = "SkyBlue1"
    SKYBLUE2 = "SkyBlue2"
    SKYBLUE3 = "SkyBlue3"
    SKYBLUE4 = "SkyBlue4"
    SLATE_BLUE = "slate blue"
    SLATE_GRAY = "slate gray"
    SLATE_GREY = "slate grey"
    SLATEBLUE = "SlateBlue"
    SLATEBLUE1 = "SlateBlue1"
    SLATEBLUE2 = "SlateBlue2"
    SLATEBLUE3 = "SlateBlue3"
    SLATEBLUE4 = "SlateBlue4"
    SLATEGRAY = "SlateGray"
    SLATEGRAY1 = "SlateGray1"
    SLATEGRAY2 = "SlateGray2"
    SLATEGRAY3 = "SlateGray3"
    SLATEGRAY4 = "SlateGray4"
    SLATEGREY = "SlateGrey"
    SNOW = "snow"
    SNOW1 = "snow1"
    SNOW2 = "snow2"
    SNOW3 = "snow3"
    SNOW4 = "snow4"
    SPRING_GREEN = "spring green"
    SPRINGGREEN = "SpringGreen"
    SPRINGGREEN1 = "SpringGreen1"
    SPRINGGREEN2 = "SpringGreen2"
    SPRINGGREEN3 = "SpringGreen3"
    SPRINGGREEN4 = "SpringGreen4"
    STEEL_BLUE = "steel blue"
    STEELBLUE = "SteelBlue"
    STEELBLUE1 = "SteelBlue1"
    STEELBLUE2 = "SteelBlue2"
    STEELBLUE3 = "SteelBlue3"
    STEELBLUE4 = "SteelBlue4"
    TAN = "tan"
    TAN1 = "tan1"
    TAN2 = "tan2"
    TAN3 = "tan3"
    TAN4 = "tan4"
    THISTLE = "thistle"
    THISTLE1 = "thistle1"
    THISTLE2 = "thistle2"
    THISTLE3 = "thistle3"
    THISTLE4 = "thistle4"
    TOMATO = "tomato"
    TOMATO1 = "tomato1"
    TOMATO2 = "tomato2"
    TOMATO3 = "tomato3"
    TOMATO4 = "tomato4"
    TURQUOISE = "turquoise"
    TURQUOISE1 = "turquoise1"
    TURQUOISE2 = "turquoise2"
    TURQUOISE3 = "turquoise3"
    TURQUOISE4 = "turquoise4"
    VIOLET = "violet"
    VIOLET_RED = "violet red"
    VIOLETRED = "VioletRed"
    VIOLETRED1 = "VioletRed1"
    VIOLETRED2 = "VioletRed2"
    VIOLETRED3 = "VioletRed3"
    VIOLETRED4 = "VioletRed4"
    WHEAT = "wheat"
    WHEAT1 = "wheat1"
    WHEAT2 = "wheat2"
    WHEAT3 = "wheat3"
    WHEAT4 = "wheat4"
    WHITE = "white"
    WHITE_SMOKE = "white smoke"
    WHITESMOKE = "WhiteSmoke"
    YELLOW = "yellow"
    YELLOW_GREEN = "yellow green"
    YELLOW1 = "yellow1"
    YELLOW2 = "yellow2"
    YELLOW3 = "yellow3"
    YELLOW4 = "yellow4"
    YELLOWGREEN = "YellowGreen"
    
    # A map from Tcl/Tk color names like "red" to the corresponding RGB values as tuples.
    # See also:
    # https://www.tcl.tk/man/tcl8.4/TkCmd/colors.htm
    # http://wiki.tcl.tk/1424
    NAME_TO_RGB = {\
        "alice blue": (240, 248, 255), \
        "AliceBlue": (240, 248, 255), \
        "antique white": (250, 235, 215), \
        "AntiqueWhite": (250, 235, 215), \
        "AntiqueWhite1": (255, 239, 219), \
        "AntiqueWhite2": (238, 223, 204), \
        "AntiqueWhite3": (205, 192, 176), \
        "AntiqueWhite4": (139, 131, 120), \
        "aquamarine": (127, 255, 212), \
        "aquamarine1": (127, 255, 212), \
        "aquamarine2": (118, 238, 198), \
        "aquamarine3": (102, 205, 170), \
        "aquamarine4": (69, 139, 116), \
        "azure": (240, 255, 255), \
        "azure1": (240, 255, 255), \
        "azure2": (224, 238, 238), \
        "azure3": (193, 205, 205), \
        "azure4": (131, 139, 139), \
        "beige": (245, 245, 220), \
        "bisque": (255, 228, 196), \
        "bisque1": (255, 228, 196), \
        "bisque2": (238, 213, 183), \
        "bisque3": (205, 183, 158), \
        "bisque4": (139, 125, 107), \
        "black": (0, 0, 0), \
        "blanched almond": (255, 235, 205), \
        "BlanchedAlmond": (255, 235, 205), \
        "blue": (0, 0, 255), \
        "blue violet": (138, 43, 226), \
        "blue1": (0, 0, 255), \
        "blue2": (0, 0, 238), \
        "blue3": (0, 0, 205), \
        "blue4": (0, 0, 139), \
        "BlueViolet": (138, 43, 226), \
        "brown": (165, 42, 42), \
        "brown1": (255, 64, 64), \
        "brown2": (238, 59, 59), \
        "brown3": (205, 51, 51), \
        "brown4": (139, 35, 35), \
        "burlywood": (222, 184, 135), \
        "burlywood1": (255, 211, 155), \
        "burlywood2": (238, 197, 145), \
        "burlywood3": (205, 170, 125), \
        "burlywood4": (139, 115, 85), \
        "cadet blue": (95, 158, 160), \
        "CadetBlue": (95, 158, 160), \
        "CadetBlue1": (152, 245, 255), \
        "CadetBlue2": (142, 229, 238), \
        "CadetBlue3": (122, 197, 205), \
        "CadetBlue4": (83, 134, 139), \
        "chartreuse": (127, 255, 0), \
        "chartreuse1": (127, 255, 0), \
        "chartreuse2": (118, 238, 0), \
        "chartreuse3": (102, 205, 0), \
        "chartreuse4": (69, 139, 0), \
        "chocolate": (210, 105, 30), \
        "chocolate1": (255, 127, 36), \
        "chocolate2": (238, 118, 33), \
        "chocolate3": (205, 102, 29), \
        "chocolate4": (139, 69, 19), \
        "coral": (255, 127, 80), \
        "coral1": (255, 114, 86), \
        "coral2": (238, 106, 80), \
        "coral3": (205, 91, 69), \
        "coral4": (139, 62, 47), \
        "cornflower blue": (100, 149, 237), \
        "CornflowerBlue": (100, 149, 237), \
        "cornsilk": (255, 248, 220), \
        "cornsilk1": (255, 248, 220), \
        "cornsilk2": (238, 232, 205), \
        "cornsilk3": (205, 200, 177), \
        "cornsilk4": (139, 136, 120), \
        "cyan": (0, 255, 255), \
        "cyan1": (0, 255, 255), \
        "cyan2": (0, 238, 238), \
        "cyan3": (0, 205, 205), \
        "cyan4": (0, 139, 139), \
        "dark blue": (0, 0, 139), \
        "dark cyan": (0, 139, 139), \
        "dark goldenrod": (184, 134, 11), \
        "dark gray": (169, 169, 169), \
        "dark green": (0, 100, 0), \
        "dark grey": (169, 169, 169), \
        "dark khaki": (189, 183, 107), \
        "dark magenta": (139, 0, 139), \
        "dark olive green": (85, 107, 47), \
        "dark orange": (255, 140, 0), \
        "dark orchid": (153, 50, 204), \
        "dark red": (139, 0, 0), \
        "dark salmon": (233, 150, 122), \
        "dark sea green": (143, 188, 143), \
        "dark slate blue": (72, 61, 139), \
        "dark slate gray": (47, 79, 79), \
        "dark slate grey": (47, 79, 79), \
        "dark turquoise": (0, 206, 209), \
        "dark violet": (148, 0, 211), \
        "DarkBlue": (0, 0, 139), \
        "DarkCyan": (0, 139, 139), \
        "DarkGoldenrod": (184, 134, 11), \
        "DarkGoldenrod1": (255, 185, 15), \
        "DarkGoldenrod2": (238, 173, 14), \
        "DarkGoldenrod3": (205, 149, 12), \
        "DarkGoldenrod4": (139, 101, 8), \
        "DarkGray": (169, 169, 169), \
        "DarkGreen": (0, 100, 0), \
        "DarkGrey": (169, 169, 169), \
        "DarkKhaki": (189, 183, 107), \
        "DarkMagenta": (139, 0, 139), \
        "DarkOliveGreen": (85, 107, 47), \
        "DarkOliveGreen1": (202, 255, 112), \
        "DarkOliveGreen2": (188, 238, 104), \
        "DarkOliveGreen3": (162, 205, 90), \
        "DarkOliveGreen4": (110, 139, 61), \
        "DarkOrange": (255, 140, 0), \
        "DarkOrange1": (255, 127, 0), \
        "DarkOrange2": (238, 118, 0), \
        "DarkOrange3": (205, 102, 0), \
        "DarkOrange4": (139, 69, 0), \
        "DarkOrchid": (153, 50, 204), \
        "DarkOrchid1": (191, 62, 255), \
        "DarkOrchid2": (178, 58, 238), \
        "DarkOrchid3": (154, 50, 205), \
        "DarkOrchid4": (104, 34, 139), \
        "DarkRed": (139, 0, 0), \
        "DarkSalmon": (233, 150, 122), \
        "DarkSeaGreen": (143, 188, 143), \
        "DarkSeaGreen1": (193, 255, 193), \
        "DarkSeaGreen2": (180, 238, 180), \
        "DarkSeaGreen3": (155, 205, 155), \
        "DarkSeaGreen4": (105, 139, 105), \
        "DarkSlateBlue": (72, 61, 139), \
        "DarkSlateGray": (47, 79, 79), \
        "DarkSlateGray1": (151, 255, 255), \
        "DarkSlateGray2": (141, 238, 238), \
        "DarkSlateGray3": (121, 205, 205), \
        "DarkSlateGray4": (82, 139, 139), \
        "DarkSlateGrey": (47, 79, 79), \
        "DarkTurquoise": (0, 206, 209), \
        "DarkViolet": (148, 0, 211), \
        "deep pink": (255, 20, 147), \
        "deep sky blue": (0, 191, 255), \
        "DeepPink": (255, 20, 147), \
        "DeepPink1": (255, 20, 147), \
        "DeepPink2": (238, 18, 137), \
        "DeepPink3": (205, 16, 118), \
        "DeepPink4": (139, 10, 80), \
        "DeepSkyBlue": (0, 191, 255), \
        "DeepSkyBlue1": (0, 191, 255), \
        "DeepSkyBlue2": (0, 178, 238), \
        "DeepSkyBlue3": (0, 154, 205), \
        "DeepSkyBlue4": (0, 104, 139), \
        "dim gray": (105, 105, 105), \
        "dim grey": (105, 105, 105), \
        "DimGray": (105, 105, 105), \
        "DimGrey": (105, 105, 105), \
        "dodger blue": (30, 144, 255), \
        "DodgerBlue": (30, 144, 255), \
        "DodgerBlue1": (30, 144, 255), \
        "DodgerBlue2": (28, 134, 238), \
        "DodgerBlue3": (24, 116, 205), \
        "DodgerBlue4": (16, 78, 139), \
        "firebrick": (178, 34, 34), \
        "firebrick1": (255, 48, 48), \
        "firebrick2": (238, 44, 44), \
        "firebrick3": (205, 38, 38), \
        "firebrick4": (139, 26, 26), \
        "floral white": (255, 250, 240), \
        "FloralWhite": (255, 250, 240), \
        "forest green": (34, 139, 34), \
        "ForestGreen": (34, 139, 34), \
        "gainsboro": (220, 220, 220), \
        "ghost white": (248, 248, 255), \
        "GhostWhite": (248, 248, 255), \
        "gold": (255, 215, 0), \
        "gold1": (255, 215, 0), \
        "gold2": (238, 201, 0), \
        "gold3": (205, 173, 0), \
        "gold4": (139, 117, 0), \
        "goldenrod": (218, 165, 32), \
        "goldenrod1": (255, 193, 37), \
        "goldenrod2": (238, 180, 34), \
        "goldenrod3": (205, 155, 29), \
        "goldenrod4": (139, 105, 20), \
        "gray": (190, 190, 190), \
        "gray0": (0, 0, 0), \
        "gray1": (3, 3, 3), \
        "gray2": (5, 5, 5), \
        "gray3": (8, 8, 8), \
        "gray4": (10, 10, 10), \
        "gray5": (13, 13, 13), \
        "gray6": (15, 15, 15), \
        "gray7": (18, 18, 18), \
        "gray8": (20, 20, 20), \
        "gray9": (23, 23, 23), \
        "gray10": (26, 26, 26), \
        "gray11": (28, 28, 28), \
        "gray12": (31, 31, 31), \
        "gray13": (33, 33, 33), \
        "gray14": (36, 36, 36), \
        "gray15": (38, 38, 38), \
        "gray16": (41, 41, 41), \
        "gray17": (43, 43, 43), \
        "gray18": (46, 46, 46), \
        "gray19": (48, 48, 48), \
        "gray20": (51, 51, 51), \
        "gray21": (54, 54, 54), \
        "gray22": (56, 56, 56), \
        "gray23": (59, 59, 59), \
        "gray24": (61, 61, 61), \
        "gray25": (64, 64, 64), \
        "gray26": (66, 66, 66), \
        "gray27": (69, 69, 69), \
        "gray28": (71, 71, 71), \
        "gray29": (74, 74, 74), \
        "gray30": (77, 77, 77), \
        "gray31": (79, 79, 79), \
        "gray32": (82, 82, 82), \
        "gray33": (84, 84, 84), \
        "gray34": (87, 87, 87), \
        "gray35": (89, 89, 89), \
        "gray36": (92, 92, 92), \
        "gray37": (94, 94, 94), \
        "gray38": (97, 97, 97), \
        "gray39": (99, 99, 99), \
        "gray40": (102, 102, 102), \
        "gray41": (105, 105, 105), \
        "gray42": (107, 107, 107), \
        "gray43": (110, 110, 110), \
        "gray44": (112, 112, 112), \
        "gray45": (115, 115, 115), \
        "gray46": (117, 117, 117), \
        "gray47": (120, 120, 120), \
        "gray48": (122, 122, 122), \
        "gray49": (125, 125, 125), \
        "gray50": (127, 127, 127), \
        "gray51": (130, 130, 130), \
        "gray52": (133, 133, 133), \
        "gray53": (135, 135, 135), \
        "gray54": (138, 138, 138), \
        "gray55": (140, 140, 140), \
        "gray56": (143, 143, 143), \
        "gray57": (145, 145, 145), \
        "gray58": (148, 148, 148), \
        "gray59": (150, 150, 150), \
        "gray60": (153, 153, 153), \
        "gray61": (156, 156, 156), \
        "gray62": (158, 158, 158), \
        "gray63": (161, 161, 161), \
        "gray64": (163, 163, 163), \
        "gray65": (166, 166, 166), \
        "gray66": (168, 168, 168), \
        "gray67": (171, 171, 171), \
        "gray68": (173, 173, 173), \
        "gray69": (176, 176, 176), \
        "gray70": (179, 179, 179), \
        "gray71": (181, 181, 181), \
        "gray72": (184, 184, 184), \
        "gray73": (186, 186, 186), \
        "gray74": (189, 189, 189), \
        "gray75": (191, 191, 191), \
        "gray76": (194, 194, 194), \
        "gray77": (196, 196, 196), \
        "gray78": (199, 199, 199), \
        "gray79": (201, 201, 201), \
        "gray80": (204, 204, 204), \
        "gray81": (207, 207, 207), \
        "gray82": (209, 209, 209), \
        "gray83": (212, 212, 212), \
        "gray84": (214, 214, 214), \
        "gray85": (217, 217, 217), \
        "gray86": (219, 219, 219), \
        "gray87": (222, 222, 222), \
        "gray88": (224, 224, 224), \
        "gray89": (227, 227, 227), \
        "gray90": (229, 229, 229), \
        "gray91": (232, 232, 232), \
        "gray92": (235, 235, 235), \
        "gray93": (237, 237, 237), \
        "gray94": (240, 240, 240), \
        "gray95": (242, 242, 242), \
        "gray96": (245, 245, 245), \
        "gray97": (247, 247, 247), \
        "gray98": (250, 250, 250), \
        "gray99": (252, 252, 252), \
        "gray100": (255, 255, 255), \
        "green": (0, 255, 0), \
        "green yellow": (173, 255, 47), \
        "green1": (0, 255, 0), \
        "green2": (0, 238, 0), \
        "green3": (0, 205, 0), \
        "green4": (0, 139, 0), \
        "GreenYellow": (173, 255, 47), \
        "grey": (190, 190, 190), \
        "grey0": (0, 0, 0), \
        "grey1": (3, 3, 3), \
        "grey2": (5, 5, 5), \
        "grey3": (8, 8, 8), \
        "grey4": (10, 10, 10), \
        "grey5": (13, 13, 13), \
        "grey6": (15, 15, 15), \
        "grey7": (18, 18, 18), \
        "grey8": (20, 20, 20), \
        "grey9": (23, 23, 23), \
        "grey10": (26, 26, 26), \
        "grey11": (28, 28, 28), \
        "grey12": (31, 31, 31), \
        "grey13": (33, 33, 33), \
        "grey14": (36, 36, 36), \
        "grey15": (38, 38, 38), \
        "grey16": (41, 41, 41), \
        "grey17": (43, 43, 43), \
        "grey18": (46, 46, 46), \
        "grey19": (48, 48, 48), \
        "grey20": (51, 51, 51), \
        "grey21": (54, 54, 54), \
        "grey22": (56, 56, 56), \
        "grey23": (59, 59, 59), \
        "grey24": (61, 61, 61), \
        "grey25": (64, 64, 64), \
        "grey26": (66, 66, 66), \
        "grey27": (69, 69, 69), \
        "grey28": (71, 71, 71), \
        "grey29": (74, 74, 74), \
        "grey30": (77, 77, 77), \
        "grey31": (79, 79, 79), \
        "grey32": (82, 82, 82), \
        "grey33": (84, 84, 84), \
        "grey34": (87, 87, 87), \
        "grey35": (89, 89, 89), \
        "grey36": (92, 92, 92), \
        "grey37": (94, 94, 94), \
        "grey38": (97, 97, 97), \
        "grey39": (99, 99, 99), \
        "grey40": (102, 102, 102), \
        "grey41": (105, 105, 105), \
        "grey42": (107, 107, 107), \
        "grey43": (110, 110, 110), \
        "grey44": (112, 112, 112), \
        "grey45": (115, 115, 115), \
        "grey46": (117, 117, 117), \
        "grey47": (120, 120, 120), \
        "grey48": (122, 122, 122), \
        "grey49": (125, 125, 125), \
        "grey50": (127, 127, 127), \
        "grey51": (130, 130, 130), \
        "grey52": (133, 133, 133), \
        "grey53": (135, 135, 135), \
        "grey54": (138, 138, 138), \
        "grey55": (140, 140, 140), \
        "grey56": (143, 143, 143), \
        "grey57": (145, 145, 145), \
        "grey58": (148, 148, 148), \
        "grey59": (150, 150, 150), \
        "grey60": (153, 153, 153), \
        "grey61": (156, 156, 156), \
        "grey62": (158, 158, 158), \
        "grey63": (161, 161, 161), \
        "grey64": (163, 163, 163), \
        "grey65": (166, 166, 166), \
        "grey66": (168, 168, 168), \
        "grey67": (171, 171, 171), \
        "grey68": (173, 173, 173), \
        "grey69": (176, 176, 176), \
        "grey70": (179, 179, 179), \
        "grey71": (181, 181, 181), \
        "grey72": (184, 184, 184), \
        "grey73": (186, 186, 186), \
        "grey74": (189, 189, 189), \
        "grey75": (191, 191, 191), \
        "grey76": (194, 194, 194), \
        "grey77": (196, 196, 196), \
        "grey78": (199, 199, 199), \
        "grey79": (201, 201, 201), \
        "grey80": (204, 204, 204), \
        "grey81": (207, 207, 207), \
        "grey82": (209, 209, 209), \
        "grey83": (212, 212, 212), \
        "grey84": (214, 214, 214), \
        "grey85": (217, 217, 217), \
        "grey86": (219, 219, 219), \
        "grey87": (222, 222, 222), \
        "grey88": (224, 224, 224), \
        "grey89": (227, 227, 227), \
        "grey90": (229, 229, 229), \
        "grey91": (232, 232, 232), \
        "grey92": (235, 235, 235), \
        "grey93": (237, 237, 237), \
        "grey94": (240, 240, 240), \
        "grey95": (242, 242, 242), \
        "grey96": (245, 245, 245), \
        "grey97": (247, 247, 247), \
        "grey98": (250, 250, 250), \
        "grey99": (252, 252, 252), \
        "grey100": (255, 255, 255), \
        "honeydew": (240, 255, 240), \
        "honeydew1": (240, 255, 240), \
        "honeydew2": (224, 238, 224), \
        "honeydew3": (193, 205, 193), \
        "honeydew4": (131, 139, 131), \
        "hot pink": (255, 105, 180), \
        "HotPink": (255, 105, 180), \
        "HotPink1": (255, 110, 180), \
        "HotPink2": (238, 106, 167), \
        "HotPink3": (205, 96, 144), \
        "HotPink4": (139, 58, 98), \
        "indian red": (205, 92, 92), \
        "IndianRed": (205, 92, 92), \
        "IndianRed1": (255, 106, 106), \
        "IndianRed2": (238, 99, 99), \
        "IndianRed3": (205, 85, 85), \
        "IndianRed4": (139, 58, 58), \
        "ivory": (255, 255, 240), \
        "ivory1": (255, 255, 240), \
        "ivory2": (238, 238, 224), \
        "ivory3": (205, 205, 193), \
        "ivory4": (139, 139, 131), \
        "khaki": (240, 230, 140), \
        "khaki1": (255, 246, 143), \
        "khaki2": (238, 230, 133), \
        "khaki3": (205, 198, 115), \
        "khaki4": (139, 134, 78), \
        "lavender": (230, 230, 250), \
        "lavender blush": (255, 240, 245), \
        "LavenderBlush": (255, 240, 245), \
        "LavenderBlush1": (255, 240, 245), \
        "LavenderBlush2": (238, 224, 229), \
        "LavenderBlush3": (205, 193, 197), \
        "LavenderBlush4": (139, 131, 134), \
        "lawn green": (124, 252, 0), \
        "LawnGreen": (124, 252, 0), \
        "lemon chiffon": (255, 250, 205), \
        "LemonChiffon": (255, 250, 205), \
        "LemonChiffon1": (255, 250, 205), \
        "LemonChiffon2": (238, 233, 191), \
        "LemonChiffon3": (205, 201, 165), \
        "LemonChiffon4": (139, 137, 112), \
        "light blue": (173, 216, 230), \
        "light coral": (240, 128, 128), \
        "light cyan": (224, 255, 255), \
        "light goldenrod": (238, 221, 130), \
        "light goldenrod yellow": (250, 250, 210), \
        "light gray": (211, 211, 211), \
        "light green": (144, 238, 144), \
        "light grey": (211, 211, 211), \
        "light pink": (255, 182, 193), \
        "light salmon": (255, 160, 122), \
        "light sea green": (32, 178, 170), \
        "light sky blue": (135, 206, 250), \
        "light slate blue": (132, 112, 255), \
        "light slate gray": (119, 136, 153), \
        "light slate grey": (119, 136, 153), \
        "light steel blue": (176, 196, 222), \
        "light yellow": (255, 255, 224), \
        "LightBlue": (173, 216, 230), \
        "LightBlue1": (191, 239, 255), \
        "LightBlue2": (178, 223, 238), \
        "LightBlue3": (154, 192, 205), \
        "LightBlue4": (104, 131, 139), \
        "LightCoral": (240, 128, 128), \
        "LightCyan": (224, 255, 255), \
        "LightCyan1": (224, 255, 255), \
        "LightCyan2": (209, 238, 238), \
        "LightCyan3": (180, 205, 205), \
        "LightCyan4": (122, 139, 139), \
        "LightGoldenrod": (238, 221, 130), \
        "LightGoldenrod1": (255, 236, 139), \
        "LightGoldenrod2": (238, 220, 130), \
        "LightGoldenrod3": (205, 190, 112), \
        "LightGoldenrod4": (139, 129, 76), \
        "LightGoldenrodYellow": (250, 250, 210), \
        "LightGray": (211, 211, 211), \
        "LightGreen": (144, 238, 144), \
        "LightGrey": (211, 211, 211), \
        "LightPink": (255, 182, 193), \
        "LightPink1": (255, 174, 185), \
        "LightPink2": (238, 162, 173), \
        "LightPink3": (205, 140, 149), \
        "LightPink4": (139, 95, 101), \
        "LightSalmon": (255, 160, 122), \
        "LightSalmon1": (255, 160, 122), \
        "LightSalmon2": (238, 149, 114), \
        "LightSalmon3": (205, 129, 98), \
        "LightSalmon4": (139, 87, 66), \
        "LightSeaGreen": (32, 178, 170), \
        "LightSkyBlue": (135, 206, 250), \
        "LightSkyBlue1": (176, 226, 255), \
        "LightSkyBlue2": (164, 211, 238), \
        "LightSkyBlue3": (141, 182, 205), \
        "LightSkyBlue4": (96, 123, 139), \
        "LightSlateBlue": (132, 112, 255), \
        "LightSlateGray": (119, 136, 153), \
        "LightSlateGrey": (119, 136, 153), \
        "LightSteelBlue": (176, 196, 222), \
        "LightSteelBlue1": (202, 225, 255), \
        "LightSteelBlue2": (188, 210, 238), \
        "LightSteelBlue3": (162, 181, 205), \
        "LightSteelBlue4": (110, 123, 139), \
        "LightYellow": (255, 255, 224), \
        "LightYellow1": (255, 255, 224), \
        "LightYellow2": (238, 238, 209), \
        "LightYellow3": (205, 205, 180), \
        "LightYellow4": (139, 139, 122), \
        "lime green": (50, 205, 50), \
        "LimeGreen": (50, 205, 50), \
        "linen": (250, 240, 230), \
        "magenta": (255, 0, 255), \
        "magenta1": (255, 0, 255), \
        "magenta2": (238, 0, 238), \
        "magenta3": (205, 0, 205), \
        "magenta4": (139, 0, 139), \
        "maroon": (176, 48, 96), \
        "maroon1": (255, 52, 179), \
        "maroon2": (238, 48, 167), \
        "maroon3": (205, 41, 144), \
        "maroon4": (139, 28, 98), \
        "medium aquamarine": (102, 205, 170), \
        "medium blue": (0, 0, 205), \
        "medium orchid": (186, 85, 211), \
        "medium purple": (147, 112, 219), \
        "medium sea green": (60, 179, 113), \
        "medium slate blue": (123, 104, 238), \
        "medium spring green": (0, 250, 154), \
        "medium turquoise": (72, 209, 204), \
        "medium violet red": (199, 21, 133), \
        "MediumAquamarine": (102, 205, 170), \
        "MediumBlue": (0, 0, 205), \
        "MediumOrchid": (186, 85, 211), \
        "MediumOrchid1": (224, 102, 255), \
        "MediumOrchid2": (209, 95, 238), \
        "MediumOrchid3": (180, 82, 205), \
        "MediumOrchid4": (122, 55, 139), \
        "MediumPurple": (147, 112, 219), \
        "MediumPurple1": (171, 130, 255), \
        "MediumPurple2": (159, 121, 238), \
        "MediumPurple3": (137, 104, 205), \
        "MediumPurple4": (93, 71, 139), \
        "MediumSeaGreen": (60, 179, 113), \
        "MediumSlateBlue": (123, 104, 238), \
        "MediumSpringGreen": (0, 250, 154), \
        "MediumTurquoise": (72, 209, 204), \
        "MediumVioletRed": (199, 21, 133), \
        "midnight blue": (25, 25, 112), \
        "MidnightBlue": (25, 25, 112), \
        "mint cream": (245, 255, 250), \
        "MintCream": (245, 255, 250), \
        "misty rose": (255, 228, 225), \
        "MistyRose": (255, 228, 225), \
        "MistyRose1": (255, 228, 225), \
        "MistyRose2": (238, 213, 210), \
        "MistyRose3": (205, 183, 181), \
        "MistyRose4": (139, 125, 123), \
        "moccasin": (255, 228, 181), \
        "navajo white": (255, 222, 173), \
        "NavajoWhite": (255, 222, 173), \
        "NavajoWhite1": (255, 222, 173), \
        "NavajoWhite2": (238, 207, 161), \
        "NavajoWhite3": (205, 179, 139), \
        "NavajoWhite4": (139, 121, 94), \
        "navy": (0, 0, 128), \
        "navy blue": (0, 0, 128), \
        "NavyBlue": (0, 0, 128), \
        "old lace": (253, 245, 230), \
        "OldLace": (253, 245, 230), \
        "olive drab": (107, 142, 35), \
        "OliveDrab": (107, 142, 35), \
        "OliveDrab1": (192, 255, 62), \
        "OliveDrab2": (179, 238, 58), \
        "OliveDrab3": (154, 205, 50), \
        "OliveDrab4": (105, 139, 34), \
        "orange": (255, 165, 0), \
        "orange red": (255, 69, 0), \
        "orange1": (255, 165, 0), \
        "orange2": (238, 154, 0), \
        "orange3": (205, 133, 0), \
        "orange4": (139, 90, 0), \
        "OrangeRed": (255, 69, 0), \
        "OrangeRed1": (255, 69, 0), \
        "OrangeRed2": (238, 64, 0), \
        "OrangeRed3": (205, 55, 0), \
        "OrangeRed4": (139, 37, 0), \
        "orchid": (218, 112, 214), \
        "orchid1": (255, 131, 250), \
        "orchid2": (238, 122, 233), \
        "orchid3": (205, 105, 201), \
        "orchid4": (139, 71, 137), \
        "pale goldenrod": (238, 232, 170), \
        "pale green": (152, 251, 152), \
        "pale turquoise": (175, 238, 238), \
        "pale violet red": (219, 112, 147), \
        "PaleGoldenrod": (238, 232, 170), \
        "PaleGreen": (152, 251, 152), \
        "PaleGreen1": (154, 255, 154), \
        "PaleGreen2": (144, 238, 144), \
        "PaleGreen3": (124, 205, 124), \
        "PaleGreen4": (84, 139, 84), \
        "PaleTurquoise": (175, 238, 238), \
        "PaleTurquoise1": (187, 255, 255), \
        "PaleTurquoise2": (174, 238, 238), \
        "PaleTurquoise3": (150, 205, 205), \
        "PaleTurquoise4": (102, 139, 139), \
        "PaleVioletRed": (219, 112, 147), \
        "PaleVioletRed1": (255, 130, 171), \
        "PaleVioletRed2": (238, 121, 159), \
        "PaleVioletRed3": (205, 104, 127), \
        "PaleVioletRed4": (139, 71, 93), \
        "papaya whip": (255, 239, 213), \
        "PapayaWhip": (255, 239, 213), \
        "peach puff": (255, 218, 185), \
        "PeachPuff": (255, 218, 185), \
        "PeachPuff1": (255, 218, 185), \
        "PeachPuff2": (238, 203, 173), \
        "PeachPuff3": (205, 175, 149), \
        "PeachPuff4": (139, 119, 101), \
        "peru": (205, 133, 63), \
        "pink": (255, 192, 203), \
        "pink1": (255, 181, 197), \
        "pink2": (238, 169, 184), \
        "pink3": (205, 145, 158), \
        "pink4": (139, 99, 108), \
        "plum": (221, 160, 221), \
        "plum1": (255, 187, 255), \
        "plum2": (238, 174, 238), \
        "plum3": (205, 150, 205), \
        "plum4": (139, 102, 139), \
        "powder blue": (176, 224, 230), \
        "PowderBlue": (176, 224, 230), \
        "purple": (160, 32, 240), \
        "purple1": (155, 48, 255), \
        "purple2": (145, 44, 238), \
        "purple3": (125, 38, 205), \
        "purple4": (85, 26, 139), \
        "red": (255, 0, 0), \
        "red1": (255, 0, 0), \
        "red2": (238, 0, 0), \
        "red3": (205, 0, 0), \
        "red4": (139, 0, 0), \
        "rosy brown": (188, 143, 143), \
        "RosyBrown": (188, 143, 143), \
        "RosyBrown1": (255, 193, 193), \
        "RosyBrown2": (238, 180, 180), \
        "RosyBrown3": (205, 155, 155), \
        "RosyBrown4": (139, 105, 105), \
        "royal blue": (65, 105, 225), \
        "RoyalBlue": (65, 105, 225), \
        "RoyalBlue1": (72, 118, 255), \
        "RoyalBlue2": (67, 110, 238), \
        "RoyalBlue3": (58, 95, 205), \
        "RoyalBlue4": (39, 64, 139), \
        "saddle brown": (139, 69, 19), \
        "SaddleBrown": (139, 69, 19), \
        "salmon": (250, 128, 114), \
        "salmon1": (255, 140, 105), \
        "salmon2": (238, 130, 98), \
        "salmon3": (205, 112, 84), \
        "salmon4": (139, 76, 57), \
        "sandy brown": (244, 164, 96), \
        "SandyBrown": (244, 164, 96), \
        "sea green": (46, 139, 87), \
        "SeaGreen": (46, 139, 87), \
        "SeaGreen1": (84, 255, 159), \
        "SeaGreen2": (78, 238, 148), \
        "SeaGreen3": (67, 205, 128), \
        "SeaGreen4": (46, 139, 87), \
        "seashell": (255, 245, 238), \
        "seashell1": (255, 245, 238), \
        "seashell2": (238, 229, 222), \
        "seashell3": (205, 197, 191), \
        "seashell4": (139, 134, 130), \
        "sienna": (160, 82, 45), \
        "sienna1": (255, 130, 71), \
        "sienna2": (238, 121, 66), \
        "sienna3": (205, 104, 57), \
        "sienna4": (139, 71, 38), \
        "sky blue": (135, 206, 235), \
        "SkyBlue": (135, 206, 235), \
        "SkyBlue1": (135, 206, 255), \
        "SkyBlue2": (126, 192, 238), \
        "SkyBlue3": (108, 166, 205), \
        "SkyBlue4": (74, 112, 139), \
        "slate blue": (106, 90, 205), \
        "slate gray": (112, 128, 144), \
        "slate grey": (112, 128, 144), \
        "SlateBlue": (106, 90, 205), \
        "SlateBlue1": (131, 111, 255), \
        "SlateBlue2": (122, 103, 238), \
        "SlateBlue3": (105, 89, 205), \
        "SlateBlue4": (71, 60, 139), \
        "SlateGray": (112, 128, 144), \
        "SlateGray1": (198, 226, 255), \
        "SlateGray2": (185, 211, 238), \
        "SlateGray3": (159, 182, 205), \
        "SlateGray4": (108, 123, 139), \
        "SlateGrey": (112, 128, 144), \
        "snow": (255, 250, 250), \
        "snow1": (255, 250, 250), \
        "snow2": (238, 233, 233), \
        "snow3": (205, 201, 201), \
        "snow4": (139, 137, 137), \
        "spring green": (0, 255, 127), \
        "SpringGreen": (0, 255, 127), \
        "SpringGreen1": (0, 255, 127), \
        "SpringGreen2": (0, 238, 118), \
        "SpringGreen3": (0, 205, 102), \
        "SpringGreen4": (0, 139, 69), \
        "steel blue": (70, 130, 180), \
        "SteelBlue": (70, 130, 180), \
        "SteelBlue1": (99, 184, 255), \
        "SteelBlue2": (92, 172, 238), \
        "SteelBlue3": (79, 148, 205), \
        "SteelBlue4": (54, 100, 139), \
        "tan": (210, 180, 140), \
        "tan1": (255, 165, 79), \
        "tan2": (238, 154, 73), \
        "tan3": (205, 133, 63), \
        "tan4": (139, 90, 43), \
        "thistle": (216, 191, 216), \
        "thistle1": (255, 225, 255), \
        "thistle2": (238, 210, 238), \
        "thistle3": (205, 181, 205), \
        "thistle4": (139, 123, 139), \
        "tomato": (255, 99, 71 ), \
        "tomato1": (255, 99, 71), \
        "tomato2": (238, 92, 66), \
        "tomato3": (205, 79, 57), \
        "tomato4": (139, 54, 38), \
        "turquoise": (64, 224, 208), \
        "turquoise1": (0, 245, 255), \
        "turquoise2": (0, 229, 238), \
        "turquoise3": (0, 197, 205), \
        "turquoise4": (0, 134, 139), \
        "violet": (238, 130, 238), \
        "violet red": (208, 32, 144), \
        "VioletRed": (208, 32, 144), \
        "VioletRed1": (255, 62, 150), \
        "VioletRed2": (238, 58, 140), \
        "VioletRed3": (205, 50, 120), \
        "VioletRed4": (139, 34, 82), \
        "wheat": (245, 222, 179), \
        "wheat1": (255, 231, 186), \
        "wheat2": (238, 216, 174), \
        "wheat3": (205, 186, 150), \
        "wheat4": (139, 126, 102), \
        "white": (255, 255, 255), \
        "white smoke": (245, 245, 245), \
        "WhiteSmoke": (245, 245, 245), \
        "yellow": (255, 255, 0), \
        "yellow green": (154, 205, 50), \
        "yellow1": (255, 255, 0), \
        "yellow2": (238, 238, 0), \
        "yellow3": (205, 205, 0), \
        "yellow4": (139, 139, 0), \
        "YellowGreen": (154, 205, 50)}
    NAME_TO_RGB_LC = {k.lower() : v for (k, v) in NAME_TO_RGB.items()}
    NAME_TO_HEX_LC = {k.lower() : ("#%02x%02x%02x" % v) for (k, v) in NAME_TO_RGB.items()}
    RGB_TO_NAME = {v : k for (k, v) in NAME_TO_RGB.items()}
    HEX_TO_NAME = {("#%02x%02x%02x" % k) : v for (k, v) in RGB_TO_NAME.items()}