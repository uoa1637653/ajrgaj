#!/bin/bash
echo
echo Should create pdf showing differences between current 
echo and previous version of the  plhdb.tex  document.
echo The perl diff algorithm may take many minutes.
echo
cd `dirname $0`
git latexdiff HEAD~1 --main plhdb.tex --output plhdb~1.pdf
touch -am diff1.command
echo ------- finished and sleeping ten
sleep 10
