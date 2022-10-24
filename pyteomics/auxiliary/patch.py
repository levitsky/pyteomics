try:
    from packaging.version import Version
except ImportError:
    from distutils.version import LooseVersion as Version

try:
    import pandas as pd
except ImportError:
    pd = None
else:
    if hasattr(pd, '_version'):
        pv = pd._version.get_versions()['version']
    else:
        pv = pd.version.version
    if Version(pv) < Version('0.17'):
        pd.DataFrame.sort_values = pd.DataFrame.sort
