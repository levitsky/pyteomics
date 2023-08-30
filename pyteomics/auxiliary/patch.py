try:
    from packaging.version import Version
except ImportError:
    from distutils.version import LooseVersion as Version
