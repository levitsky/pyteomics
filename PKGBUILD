# Maintainer: Lev Levitsky <levlev at mail dot ru>
pkgname=python-pyteomics
pkgver=3.3
pkgrel=1
pkgdesc="A framework for proteomics data analysis."
arch=('any')
url="http://pythonhosted.org/pyteomics"
license=('Apache')
depends=('python' 'python-setuptools')
optdepends=('python-matplotlib: for pylab_aux module'
            'python-sqlalchemy: for mass.unimod module'
            'python-pandas: for convenient filtering of CSV tables from search engines'
            'python-lxml: for XML parsing modules'
            'python-numpy: for most of features, highly recommended')
options=(!emptydirs)
source=("https://pypi.python.org/packages/95/70/a3cf1ce3b5a3b9ea403c3c3f288aa29108aed5feb6d9339df1fa4b8ff322/pyteomics-3.3.tar.gz")
md5sums=('32b9a0fdced9d2f4434d733b3ec4920d')
changelog="CHANGELOG"
package() {
  cd "${srcdir}/pyteomics-${pkgver}"
  python setup.py install --root="$pkgdir/" --optimize=1
}

# vim:set ts=2 sw=2 et:
