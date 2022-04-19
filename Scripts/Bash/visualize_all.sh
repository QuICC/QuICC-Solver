
for n in $(ls | grep -oP 'state\K([0-9]{4})')
do
    echo $n
    filein=state$n.hdf5
    fileout=visState$n.hdf5
    echo $filein
    echo $fileout
    ln $filein state4Visu.hdf5 -f
    ./BoussinesqShellNuttatingCouetteExplicitVisu
    ln visState0000.hdf5 $fileout -f
    python ../../Scripts/Python/pole_patcher.py $fileout
    python ../../Scripts/Python/createXDMF.py -i $fileout --with-components

done
