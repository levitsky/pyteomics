#!/bin/sh

ver=$(cat VERSION)
sed -i "s/^pkgver=.*/pkgver=${ver}/" PKGBUILD
mkaurball -f
sed -i.old -E "/^source/!s/python(\b)/python2\1/g" PKGBUILD
mkaurball -f
mv PKGBUILD.old PKGBUILD
rm pyteomics-${ver}.tar.gz
