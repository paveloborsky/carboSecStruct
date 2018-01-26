PROGRAM=carbosecstruc
NAME1=(GlcA GlcB ManA ManB)
CODE1_1=(down up down up)
CODE1_2=(down down up up)

#         |link 1                 |link 2                                         |link 3                                         
aNAME2=(   GlcA  GlcB  ManA  ManB  ManB  GlcB  ManA  GlcA  AltA  AllA  AltB  AllB  AllB  GlcB  GulB  GalB  AltB  ManB  IdoB  TalB  )
aCODE2_0=( no    no    no    no    up    up    down  down  down  down  up    up    down  down  down  down  up    up    up    up    )
aCODE2_1=( down  up    down  up    up    down  up    down  up    down  up    down  down  up    down  up    down  up    down  up    )
aCODE2_2=( down  down  up    up    up    up    up    up    down  down  down  down  down  down  up    up    down  down  up    up    )
aLINK=(    1     1     1     1     2     2     2     2     2     2     2     2     3     3     3     3     3     3     3     3     )

#         |link 4                                             |link 6 
bNAME2=(   GalB  GlcB  LAltA  LIdoA  GulB  AllB  LManA  LTalA  LIdoA  GlcB  LAltA  GalB )
bCODE2_0=( up    up    up     up     down  down  down   down   down   down  up     up   )
bCODE2_1=( up    down  up     down   up    down  up     down   down   up    down   up   )
bCODE2_2=( up    up    down   down   up    up    down   down   no     no    no     no   )
bLINK=(    4     4     4      4      4     4     4      4      6      6     6      6    )

NAME2=(   ${aNAME2[@]}    ${bNAME2[@]}   )
CODE2_0=( ${aCODE2_0[@]}  ${bCODE2_0[@]} )
CODE2_1=( ${aCODE2_1[@]}  ${bCODE2_1[@]} )
CODE2_2=( ${aCODE2_2[@]}  ${bCODE2_2[@]} )
LINK=(    ${aLINK[@]}     ${bLINK[@]}    )

#CODET1=(GA4P GB4P MA2P MB2P)
#CODET2=(CGAP CGBP MB2P GB2P NB3P GB3P LB4P GB4P gA6P GB6P)

for ((f=0; f < ${#NAME1[@]}; f++)); do # f as the first
  for ((s=0; s < ${#NAME2[@]}; s++)); do # s as the second
    name=${NAME1[$f]}${LINK[$s]}${NAME2[$s]}

    echo "==============================="
    echo "==============================="
    echo $name
    echo "==============================="
    echo "==============================="
    l0=$((${LINK[$s]}-1))
    l1=${LINK[$s]}
    l2=$((${LINK[$s]}+1))
    if [ $l0 == 5 ]; then # to fix l0 of 1-6 linkage
      l0=4
    fi

    #========input==========#
    cat -<< _eof > ./${name}.ch
    add 1 ${LINK[$s]} 1
    ena 1 O1    ${CODE1_1[$f]}
    ena 1 O2    ${CODE1_2[$f]}
    ena 2 O$l1  ${CODE2_1[$s]}
_eof
    if [ "${CODE2_0[$s]}" != "no" ]; then
      echo "ena 2 O$l0  ${CODE2_0[$s]}" >> ./${name}.ch
    fi

    if [ "${CODE2_2[$s]}" != "no" ]; then
      echo "ena 2 O$l2  ${CODE2_2[$s]}" >> ./${name}.ch
    fi
    

    sed -i -e 's/O5/C6/g' ${name}.ch # to fix l2 of 1-4 linkage
    sed -i -e 's/O6/C6/g' ${name}.ch # to fix l1 of 1-6 linkage
    #=======end input========#

    
    #=======run=========#
    #echo "$PROGRAM -m GLCB.ang phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang"
    $PROGRAM -m ../lib/GLCB.ang ../phex.tpl  ${name}.cnf  ${name}.ch ${name}.ang

    #change name in cnf
#    sed -i -e "s/1 GB4P/1 ${CODET1[$f]}/g" ${name}.cnf
#    sed -i -e "s/2 GB4P/2 ${CODET2[$s]}/g" ${name}.cnf
    sed -i -e "s/# CONVERTED  CNF FILE/$name/g" ${name}.cnf

    #=======clean=========#
    rm -f ./${name}.ch
#    rm -f ./${name}.cnf

  done
done
