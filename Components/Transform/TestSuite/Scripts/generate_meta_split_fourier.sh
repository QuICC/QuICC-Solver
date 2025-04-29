#!/bin/bash
#
# RUN FROM build/Components/Transform/TestSuite

scheme="WLFl"
stage=2

cleardata=true
copy2Ref=true

smallDbs=( 104 )
bigDbs=( 108 )

data_dir="_data/Framework/LoadSplitter/${scheme}/Tubular"
ref_dir="_refdata/Transform/Fourier/Mixed/Integrator/"
exe="../../Framework/TestSuite/FrameworkLoadSplitterTests"

logfile="generate_meta_split_fourier.log"

# Create path if it doesn't exist
mkdir -p "${data_dir}"
mkdir -p "${ref_dir}"

# Clear data
if [ "$clearData" = true ];
then
  echo "Clearing old data"
  rm "./${data_dir}"/*
fi

# DB 104
for db in ${smallDbs[@]};
do
  for np in 4 8;
  do
    for r in 0 $(( (${np} - 1)/2 )) $(( ${np} - 1 ));
    do
      echo "Generating metadata for ${db}: np = ${db}, r = ${r}"
      ${exe} \[${scheme}\] --algorithm tubular --np ${np} --db ${db} --stage ${stage} --id ${r} --dumpData > ${logfile}
    done
  done
done

# DB 108, 109, 110, 111
for db in ${bigDbs[@]};
do
  for n in {0..10};
  do
    for np in $(( 2**n )) $(( 2**n*12 )) $(( 2**n*18 ));
    do
      if [ "$np" -gt "3" ];
      then
        for r in 0 3 7 $(( ${np} - 1 ));
        do
          if [ "$np" -gt "$r" ];
          then
            echo "Generating metadata for ${db}: np = ${db}, r = ${r}"
            ${exe} \[${scheme}\] --algorithm tubular --np ${np} --db ${db} --stage ${stage} --id ${r} --dumpData > ${logfile}
          fi
        done
      fi
    done
  done
done

if [ "$copy2Ref" = true ];
then
  echo "Copying new data to ref"
  cp -r "./${data_dir}"/* "./${ref_dir}"/
fi
