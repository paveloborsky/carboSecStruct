PROGRAM=carbosecstruc
NAME1=(GlcA GlcB)
NAME2=(GlcA GlcB ManB GlcB AllB GlcB GalB GlcB LIdoA GlcB)
CODE1=(down up)
CODE2=(down up up down down up up down down up)
LINK=(1 1 2 2 3 3 4 4 6 6)
CODET1=(GA4P GB4P)
CODET2=(CGAP CGBP MB2P GB2P NB3P GB3P LB4P GB4P gA6P GB6P)

for f in 0 1 # f as first
do
  for s in {0..9} # s as second
  do
    name=${NAME1[$f]}${LINK[$s]}${NAME2[$s]}


    #========input==========#
    cat -<< _eof > ./${name}.ch
    add 1 ${LINK[$s]} 1
    ena 1 O1 ${CODE1[$f]}
    ena 2 O${LINK[$s]} ${CODE2[$s]}
_eof
    #=======end input========#

    if ((${LINK[$s]} == 6))
    then
      sed -i -e 's/O6/C6/g' ${name}.ch
    fi

    #=======run=========#
    #echo "$PROGRAM -m GLCB.ang phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang"
    $PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang

    #change name in cnf
    sed -i -e "s/1 GB4P/1 ${CODET1[$f]}/g" ${name}.cnf
    sed -i -e "s/2 GB4P/2 ${CODET2[$s]}/g" ${name}.cnf
    sed -i -e "s/# CONVERTED  CNF FILE/$name/g" ${name}.cnf

    #=======clean=========#
    rm -f ./${name}.ch
#    rm -f ./${name}.cnf

  done
done
