# Licensed under a 3-clause BSD style license - see LICENSE.rst

from regions.shapes import (CircleAnnulusPixelRegion, CirclePixelRegion,
                            EllipseAnnulusPixelRegion, EllipsePixelRegion,
                            PointPixelRegion, PolygonPixelRegion,
                            RectanglePixelRegion)

__all__ = []


# mappings from FITS shape to region class and column names/indices
shape_map = {'point': (PointPixelRegion, ('X0', 'Y0')),
             'circle': (CirclePixelRegion, ('X0', 'Y0', 'R0')),
             'ellipse': (EllipsePixelRegion,
                         ('X0', 'Y0', 'R0', 'R1', 'ROTANG0')),
             'annulus': (CircleAnnulusPixelRegion, ('X0', 'Y0', 'R0', 'R1')),
             'elliptannulus': (EllipseAnnulusPixelRegion,
                               ('X0', 'Y0', 'R0', 'R1', 'R2', 'R3',
                                'ROTANG0')),
             'box': (RectanglePixelRegion, ('X0', 'Y0', 'R0', 'R1')),
             'rotbox': (RectanglePixelRegion,
                        ('X0', 'Y0', 'R0', 'R1', 'ROTANG0')),
             'rectangle': (RectanglePixelRegion, ('X0', 'X1', 'Y0', 'Y1')),
             'rotrectangle': (RectanglePixelRegion,
                              ('X0', 'X1', 'Y0', 'Y1', 'ROTANG0')),
             'polygon': (PolygonPixelRegion, ('X', 'Y'))}


class FITSParserError(Exception):
    """
    A custom exception for FITS parsing errors.
    """
