PROGRAM=carbosecstruc
#################################################
##########  CHANGING SECOND RESIDUE  ############
NAME1=(GlcB)
NAME2=(ManB GlcA)
CODE1_3=(up   )
CODE1_4=(down )
CODE1_6=(up   )
CODE2_1=(up down)
CODE2_2=(up down)
LINK=(4 4)
CODET1=(GB4P)
CODET2=(MB4P GA4P)
for ((f=0; f < ${#NAME1[@]}; f++)); do # f as the first
  for ((s=0; s < ${#NAME2[@]}; s++)); do # s as second

    name=${NAME1[$f]}${LINK[$s]}${NAME2[$s]}
    echo "----------------------------------"
    echo "---------${name}----------"

    #========input==========#
    cat -<< _eof > ./${name}.ch
    add 1 ${LINK[$s]} 1
    ena 1 O3 ${CODE1_3[$f]}
    ena 1 O4 ${CODE1_4[$f]}
    ena 1 C6 ${CODE1_6[$f]}
    ena 2 O1 ${CODE2_1[$s]}
    ena 2 O2 ${CODE2_2[$s]}
_eof
    #=======end input========#


    #=======run=========#
    echo "$PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang"
    $PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang

    #change name in cnf
    sed -i -e "s/1 GB4P/1 ${CODET1[$f]}/g" ${name}.cnf
    sed -i -e "s/2 GB4P/2 ${CODET2[$s]}/g" ${name}.cnf
    sed -i -e "s/# CONVERTED  CNF FILE/$name/g" ${name}.cnf

    #=======clean=========#
    rm -f ./${name}.ch
    rm -f ./${name}.cnf

  done
done

#################################################
##########  CHANGING FIRST RESIDUE  ############
NAME1=(AllB GalB LIdoA)
NAME2=(GlcB)
CODE1_3=(down up up  )
CODE1_4=(down up down)
CODE1_6=(up   up down)
CODE2_1=(up)
CODE2_2=(down)
LINK=(4 4 4 4)
CODET1=(NB4P LA4P iA4P)
CODET2=(GB4P)
for ((f=0; f < ${#NAME1[@]}; f++)); do # f as the first
  for ((s=0; s < ${#NAME2[@]}; s++)); do # s as second

    name=${NAME1[$f]}${LINK[$s]}${NAME2[$s]}
    echo "----------------------------------"
    echo "---------${name}----------"

    #========input==========#
    cat -<< _eof > ./${name}.ch
    add 1 ${LINK[$s]} 1
    ena 1 O3 ${CODE1_3[$f]}
    ena 1 O4 ${CODE1_4[$f]}
    ena 1 C6 ${CODE1_6[$f]}
    ena 2 O1 ${CODE2_1[$s]}
    ena 2 O2 ${CODE2_2[$s]}
_eof
    #=======end input========#


    #=======run=========#
    echo "$PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang"
    $PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang

    #change name in cnf
    sed -i -e "s/1 GB4P/1 ${CODET1[$f]}/g" ${name}.cnf
    sed -i -e "s/2 GB4P/2 ${CODET2[$s]}/g" ${name}.cnf
    sed -i -e "s/# CONVERTED  CNF FILE/$name/g" ${name}.cnf

    #=======clean=========#
    rm -f ./${name}.ch
    rm -f ./${name}.cnf

  done
done


