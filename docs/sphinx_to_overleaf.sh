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
fd \\.tex   --exec perl -0777 -pi -e 's/\\sphinxtableofcontents//gs'    
fd \\.tex   --exec perl -0777 -pi -e 's/\\sphinxincludegraphics(.*?)\{\{([^\}]*)\}(\.png)\}/\\sphinxincludegraphics$1\{\{sphinx\/$2\}$3\}/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\n*\\subsubsubsection\{Utility API\}/\\\\ \n\n\\textbf\{Utility API\} \\\\/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\n*\\subsubsection\{Utility API\}/\\\\ \n\n\\textbf\{Utility API\} \\\\/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\n*\\subsection\{Utility API\}/\\\\ \n\n\\textbf\{Utility API\} \\\\/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\pagestyle\{empty\}//gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\pagestyle\{plain\}//gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\pagestyle\{normal\}//gs'
# fd \\.tex   --exec perl -0777 -pi -e 's/\[width\=\d*\\sphinxpxdimen\]//gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\[width\=601\\sphinxpxdimen\]/\[width\=210\\sphinxpxdimen\]/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\sphinxcaption/\\caption/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\begin{savenotes}/\\begin{figure}\[H\]\\begin{savenotes}/gs'
fd \\.tex   --exec perl -0777 -pi -e 's/\\end{savenotes}/\\end{savenotes}\\end{figure}/gs'


# rm -rf ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx
# mkdir ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx
cp -r latex/* ~/Dropbox/Apps/Overleaf/SOXS-SPE-0022-Pipeline-Design/sphinx/
cp -r latex/* ~/Dropbox/Apps/Overleaf/SOXS-MAN-0006-Pipeline-User-Manual/sphinx/

