# Maintainer: Lev Levitsky <levlev at mail dot ru>
pkgname=python-pyteomics
pkgver=3.4.2
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
source=("https://pypi.io/packages/source/p/pyteomics/pyteomics-${pkgver}.tar.gz")
md5sums=('c44663bdf812114da10c87b3718220ed')
changelog="CHANGELOG"
package() {
  cd "${srcdir}/pyteomics-${pkgver}"
  python setup.py install --root="$pkgdir/" --optimize=1
}

# vim:set ts=2 sw=2 et:
