import json
import numpy


class MivotJsonEncoder(json.JSONEncoder):
    """
    Custom JSON encoder for NumPy data types.
    This class extends the default JSONEncoder to handle NumPy integers,
    floating-point numbers, and arrays during JSON encoding.
    """
    def default(self, obj):
        """
        Serialize NumPy data types to their Python equivalents for JSON encoding.
        Parameters
        ----------
        obj : Any
            The object to be encoded.
        Returns
        -------
        Any
            The serialized version of the object.
        """
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MivotJsonEncoder, self).default(obj)
