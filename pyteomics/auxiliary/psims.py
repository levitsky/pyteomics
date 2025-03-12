from functools import partial


try:
    from psims.controlled_vocabulary.controlled_vocabulary import (load_psimod, load_xlmod, load_gno, obo_cache, load_unimod)
    _has_psims = True
except ImportError:
    def _needs_psims(name):
        raise ImportError("Loading %s requires the `psims` library. To access it, please install `psims`" % name)

    load_psimod = partial(_needs_psims, 'PSIMOD')
    load_xlmod = partial(_needs_psims, 'XLMOD')
    load_gno = partial(_needs_psims, 'GNO')
    load_unimod = partial(_needs_psims, 'UNIMOD')
    obo_cache = None
    _has_psims = False
