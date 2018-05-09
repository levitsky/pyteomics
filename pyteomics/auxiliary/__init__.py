try:
    basestring = basestring
except NameError:
    basestring = (str, bytes)

from . import patch as __patch

from .structures import (
    PyteomicsError, Charge, ChargeList,
    _parse_charge, BasicComposition,
    unitfloat, unitint, unitstr, cvstr,
    cvquery)

from .constants import _nist_mass

from .file_helpers import (
    _file_obj, _keepstate, _keepstate_method, IteratorContextManager,
    FileReader, _file_reader, _file_writer, _make_chain)

from .math import (
    linear_regression, linear_regression_perpendicular,
    linear_regression_vertical)

from .target_decoy import (
    _calculate_qvalues, _qvalues_df, _decoy_or_pep_label,
    _construct_dtype, _make_qvalues, _make_filter,
    _itercontext, _iter, qvalues, filter, log_factorial,
    _expectation, _confidence_value, _log_pi_r,
    _log_pi, _make_fdr, fdr)

from .utils import (
    print_tree, memoize, BinaryDataArrayTransformer,
    _decode_base64_data_array)
