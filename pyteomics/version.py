"""
version - Pyteomics version information
=======================================

This module is provided for convenience and captures information about the current version number of Pyteomics.

Constants
---------

  :py:const:`version` - a string with the current version.

  :py:const:`version_info` - a tuple with structured information about the current version.

"""

__version__ = '4.5.6'

from collections import namedtuple
import re


class _VersionInfo(namedtuple('_VersionInfo', ('major', 'minor', 'micro', 'releaselevel', 'serial'))):
    """Tuple mimicking :py:const:`sys.version_info`"""
    def __new__(cls, version_str):
        if isinstance(version_str, str):
            groups = re.match(r'(\d+)\.(\d+)(?:\.)?(\d+)?([a-zA-Z]+)?(\d+)?', version_str).groups()
            inst = super(_VersionInfo, cls).__new__(cls, *groups)
        else:
            inst = super(_VersionInfo, cls).__new__(cls, *(str(x) if x is not None else x for x in version_str))
        inst._version_str = version_str
        inst._version_ints = tuple(int(x) if isinstance(x, str) and x.isdigit() else x for x in inst)
        return inst

    def __str__(self):
        return 'Pyteomics version {}'.format(self._version_str)

    def __lt__(self, other):
        if not isinstance(other, _VersionInfo):
            other = _VersionInfo(other)
        return self._version_ints < other._version_ints

    def __gt__(self, other):
        if not isinstance(other, _VersionInfo):
            other = _VersionInfo(other)
        return self._version_ints > other._version_ints

    def __le__(self, other):
        return self == other or self < other

    def __ge__(self, other):
        return self == other or self > other

    def __eq__(self, other):
        if not isinstance(other, _VersionInfo):
            other = _VersionInfo(other)
        return super(_VersionInfo, self).__eq__(other)


version_info = _VersionInfo(__version__)
version = __version__
