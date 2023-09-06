#!/bin/bash

scheme="WLFl"
stage=0

data_dir="_data/Framework/LoadSplitter/${scheme}/Tubular"
ref_dir="_refdata/Transform/Worland/Integrator/"
exe="../../Framework/TestSuite/FrameworkLoadSplitterTests"

# Create path if it doesn't exist
mkdir -p "${data_dir}"
mkdir -p "${ref_dir}"

# Clean data
rm "./${data_dir}"/*

# DB 104
db=104
for np in 4;
do
  for r in 0 1 $((${np} -1));
  do
    ${exe} \[${scheme}\] --algorithm tubular --np ${np} --db ${db} --stage ${stage} --id ${r} --dumpData
  done
done
for np in 8;
do
  for r in 0 3 $((${np} -1));
  do
    ${exe} \[${scheme}\] --algorithm tubular --np ${np} --db ${db} --stage ${stage} --id ${r} --dumpData
  done
done

# DB 108, 109, 110, 111
for db in 108 109 110 111;
do
  for n in {0..10};
  do
    for np in $(( 2**n )) $(( 2**n*12 )) $(( 2**n*18 ));
    do
      if [ "$np" -gt "7" ];
      then
        for r in 0 3 7 $((${np}-1));
        do
          ${exe} \[${scheme}\] --algorithm tubular --np ${np} --db ${db} --stage ${stage} --id ${r} --dumpData
        done
      fi
    done
  done
done

#cp -r "./${data_dir}"/* "./${ref_dir}"/
