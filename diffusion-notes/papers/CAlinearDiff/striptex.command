#!/bin/bash
# Very simple strip of LaTeX from LaTeX+Reduce code.
# Requires begin/ends to be on a line by themselves.
# Tony Roberts, 17 Jul 2013
sed '{1,/\\begin{reduce}/d
      /\\end{reduce}/,/\\begin{reduce}/d
     }' $1.tex > $1.red
