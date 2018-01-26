PROGRAM=carbosecstruc
NAME1=(GlcB ManB)
NAME2=(GlcB AllB LIdoA LTalA)
CODE2=(down up)
CODE3=(up down up   down)
CODE6=(up up   down down)
CODET1=(GB4P MB4P)
CODET2=(GB4P NB4P iA4P tA4P)
LINK=(4 4 4 4)

for f in 0 1 # f as first
do
  for s in {0..3} # s as second
  do
    name=${NAME1[$f]}${LINK[$s]}${NAME2[$s]}

    #skip GlcB1GlcB it exists
    if [ $f -eq 0 ] && [ $s -eq 0 ]; then
      continue
    fi

    #========input==========#
    cat -<< _eof > ./${name}.ch
    add 1 ${LINK[$s]} 1
    ena 1 O2 ${CODE2[$f]}
    ena 2 O3 ${CODE3[$s]}
    ena 2 C6 ${CODE6[$s]}
_eof
    #=======end input========#


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
