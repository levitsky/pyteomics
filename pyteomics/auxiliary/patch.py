from distutils.version import LooseVersion

try:
    import pandas as pd
except ImportError:
    pd = None
else:
    if hasattr(pd, '_version'):
        pv = pd._version.get_versions()['version']
    else:
        pv = pd.version.version
    if LooseVersion(pv) < LooseVersion('0.17'):
        pd.DataFrame.sort_values = pd.DataFrame.sort
