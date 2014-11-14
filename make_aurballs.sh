ver=$(cat VERSION)
sed -i "s/^pkgver=.*/pkgver=${ver}/" PKGBUILD
mkaurball
sed -i.old -E "/^source/!s/python(\b)/python2\1/g" PKGBUILD
mkaurball
mv PKGBUILD.old PKGBUILD
