#!/bin/sh

ver=$(cat VERSION)
echo "Generating source distribution ..."
python setup.py sdist
ln -s dist/pyteomics-${ver}.tar.gz pyteomics-${ver}.tar.gz
echo "Calculating MD5 sum ..."
md5=$(md5sum pyteomics-${ver}.tar.gz | cut -d' ' -f1)
sed -i "s/^pkgver=.*/pkgver=${ver}/" PKGBUILD
sed -i "s/^md5sums=.*/md5sums=('${md5}')/" PKGBUILD
echo "Generating Python 3 AUR ball ..."
mkaurball -f
echo "Patching PKGBUILD ..."
sed -i.old -E "/^source/!s/python(\b)/python2\1/g" PKGBUILD
echo "Generating Python 2 AUR ball ..."
mkaurball -f
echo "Restoring PKGBUILD ..."
mv PKGBUILD.old PKGBUILD
rm pyteomics-${ver}.tar.gz
