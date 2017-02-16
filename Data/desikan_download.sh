#/bin/bash

wget -r --no-parent http://openconnecto.me/mrdata/share/dti/ndmg_v0011/$1/desikan/
mv openconnecto.me/mrdata/share/dti/ndmg_v0011/$1/desikan desikan/$1
rm -rf openconnecto.me
