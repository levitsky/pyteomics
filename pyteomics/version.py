__version__ = '4.4.0dev1'

from collections import namedtuple
import re


class _VersionInfo(namedtuple('_VersionInfo', ('major', 'minor', 'micro', 'releaselevel', 'serial'))):
    def __new__(cls, version_str):
        groups = re.match(r'(\d+)\.(\d+)(?:\.)?(\d+)?([a-zA-Z]+)?(\d+)?', version_str).groups()
        inst = super(_VersionInfo, cls).__new__(cls, *groups)
        inst._version_str = version_str
        return inst

    def __str__(self):
        return 'Pyteomics version {}'.format(self._version_str)


version_info = _VersionInfo(__version__)
version = __version__
