"""
MIVOT vocabulary and regular expressions.
"""
import re
from astropy import units as u
from pyvo.utils import prototype_feature


class Constant:
    """
    Class used to set constant to identify XML attributes added to the MIVOT ATTRIBUTES
    """
    FIRST_TABLE = "first_table"
    FIELD_UNIT = "field_unit"
    COL_INDEX = "col_index"
    ROOT_COLLECTION = "root_collection"
    ROOT_OBJECT = "root_object"
    NOT_SET = "NotSet"
    ANONYMOUS_TABLE = "AnonymousTable"


# Regexp pattern to check no valid mapping is present
NoMapping = re.compile(r".REPORT\s+status=['\"]KO")

unit_mapping = {
    "deg": u.degree,
    "rad": u.radian,
    "hourangle": u.hourangle,
    "arcsec": u.arcsec,
    "mas": u.mas,
    "pc": u.pc,
    "km": u.km,
    "m": u.m,
    "mas/yr": u.mas / u.yr,
    "mas/y": u.mas / u.yr,
    "km/s": u.km / u.s,
}


@prototype_feature('MIVOT')
class Ele:
    """
    Constant used to identify MIVOT Element
    """
    namespace = ""
    VODML = namespace + "VODML"
    MODEL = namespace + "MODEL"
    GLOBALS = namespace + "GLOBALS"
    TEMPLATES = namespace + "TEMPLATES"
    INSTANCE = namespace + "INSTANCE"
    ATTRIBUTE = namespace + "ATTRIBUTE"
    COLLECTION = namespace + "COLLECTION"
    JOIN = namespace + "JOIN"
    REFERENCE = namespace + "REFERENCE"
    WHERE = namespace + "WHERE"
    NOROLE = "NOROLE"


@prototype_feature('MIVOT')
class Att:
    """
    Constant used to identify attributes in MIVOT Element
    """
    dmrole = "dmrole"
    dmtype = "dmtype"
    dmid = "dmid"
    name = "name"
    value = "value"
    dmref = "dmref"
    tableref = "tableref"
    sourceref = "sourceref"
    ref = "ref"
    primarykey = "PRIMARY_KEY"
    foreignkey = "foreignkey"
