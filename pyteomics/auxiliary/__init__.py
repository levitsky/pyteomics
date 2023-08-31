try:
    basestring = basestring
except NameError:
    basestring = (str, bytes)

from .structures import (
    PyteomicsError, Charge, ChargeList,
    _parse_charge, _parse_ion, BasicComposition,
    unitfloat, unitint, unitstr, cvstr,
    cvquery)

from .constants import _nist_mass

from .file_helpers import (
    _file_obj, _keepstate, _keepstate_method, IteratorContextManager,
    FileReader, IndexedTextReader, IndexedReaderMixin, TimeOrderedIndexedReaderMixin,
    IndexSavingMixin, OffsetIndex, HierarchicalOffsetIndex, IndexSavingTextReader,
    _file_reader, _file_writer,
    _make_chain, _check_use_index, FileReadingProcess, TaskMappingMixin,
    serializer, ChainBase, TableJoiner)

from .math import (
    linear_regression, linear_regression_perpendicular,
    linear_regression_vertical)

from .target_decoy import (
    _calculate_qvalues, _qvalues_df, _decoy_or_pep_label,
    _construct_dtype, _make_qvalues, _make_filter,
    _itercontext, _iter, qvalues, filter, log_factorial,
    _expectation, _confidence_value, _log_pi_r,
    _log_pi, _make_fdr, fdr, sigma_T, sigma_fdr)

from .utils import (
    print_tree, memoize, BinaryDataArrayTransformer, ArrayConversionMixin, BinaryArrayConversionMixin,
    MaskedArrayConversionMixin, _decode_base64_data_array)
