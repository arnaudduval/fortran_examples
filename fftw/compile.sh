# !/bin/bash

# Pass -legacy as option to compare with fft routine from Numerical Recipes


if [[ $1 = "-legacy" ]]
then
    gfortran love_coef.F90 fft.for -o prog -lfftw3f -L/usr/lib/x86_64-linux-gnu -I /usr/include -DLEGACY
else
    gfortran love_coef.F90 -o prog -lfftw3f -L/usr/lib/x86_64-linux-gnu -I /usr/include
fi