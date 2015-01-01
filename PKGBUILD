# Maintainer: Lev Levitsky <levlev at mail dot ru>
pkgname=python-pyteomics
pkgver=2.5.4
pkgrel=1
pkgdesc="A framework for proteomics data analysis."
arch=('any')
url="http://pythonhosted.org/pyteomics"
license=('Apache')
depends=('python' 'python-lxml' 'python-numpy' )
optdepends=('python-matplotlib: for pylab_aux module')
options=(!emptydirs)
source=("https://pypi.python.org/packages/source/p/pyteomics/pyteomics-${pkgver}.tar.gz")
md5sums=('a6d1164d6bf1364f218899bc88232c87')
changelog="CHANGELOG"
package() {
  cd "${srcdir}/pyteomics-${pkgver}"
  python setup.py install --root="$pkgdir/" --optimize=1
}

# vim:set ts=2 sw=2 et:
