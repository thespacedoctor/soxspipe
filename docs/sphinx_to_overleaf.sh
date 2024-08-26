#! /bin/bash



source "/Users/Dave/anaconda/etc/profile.d/conda.sh" > /dev/null 2>&1
conda activate soxspipe
cd ~/git_repos/_packages_/python/soxspipe/docs/
rm -rf build
make latex
make html
cd build
fd \\.sty   --exec perl -0777 -pi -e 's/\\input\{(sphinx[^\/].*?.sty)}/\\input\{sphinx\/$1\}/gs'
fd \\.sty   --exec perl -0777 -pi -e 's/\\RequirePackage\{(sphinx[^\/].*?)}/\\RequirePackage\{sphinx\/$1\}/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\sphinxtableofcontents//gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\sphinxmaketitle//gs'
# rm -rf ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx
# mkdir ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx
cp -r latex/* ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx/
