# bibtex ids with colons in them (eg R package names) don't work. Remove the pattern :_
sed -i "" "s/:_/_/g" My\ Library.bib

pandoc -o draft.tex --template=mytemplate.latex --bibliography=My\ Library.bib --csl=ecology.csl draft.Rmd

pdflatex draft.tex

pandoc -o Appendices.pdf --bibliography=My\ Library.bib --csl=ecology.csl Appendices.Rmd

open draft.pdf
open Appendices.pdf
