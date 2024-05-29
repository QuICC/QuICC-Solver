#!/bin/bash
#
# RUN FROM build/Components/Framework/TestSuite

genAllRanks=true
clearRef=true
copy2Ref=true
genUniform=true
genTriangular=false

# uniform truncation
if [ "$genUniform" = true ];
then
  for sch in "WLFl"; # "SLFl" "WLFm" "SLFm" "TFF";
  do
    data_dir="_data/Framework/LoadSplitter/${sch}/"
    ref_dir="_refdata/Framework/LoadSplitter/${sch}/"

  # Create path if it doesn't exist
  mkdir -p "${data_dir}"
  mkdir -p "${ref_dir}"

  # Clean data
  rm "./${data_dir}"/{Serial,Single1D,Single2D,Tubular}/*
  if [ "$clearRef" = true ];
  then
    rm -r "./{${ref_dir}"/{Serial,Single1D,Single2D,Tubular}
  fi

  # Serial data
  mkdir -p "${data_dir}"/Serial
  mkdir -p "${ref_dir}"/Serial
  for db in 103 104 106 108;
  do
    for st in 0 1 2 3;
    do
      ./FrameworkLoadSplitterTests \[${sch}\] --algorithm serial --np 1 --db ${db} --stage ${st} --dumpData
    done
  done

  # Single1D data
  mkdir -p "${data_dir}/Single1D"
  mkdir -p "${ref_dir}/Single1D"
  # for db in 103;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 104;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 15 16;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 106;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 16 42 64;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 108;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 32 144 256;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # # Generate all ranks
  # if [ "$genAllRanks" = true ];
  # then
  #   for db in 108;
  #   do
  #     for st in 0 1 2 3;
  #     do
  #       for np in 56;
  #       do
  #         ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --id ${np} --np ${np} --db ${db} --stage ${st} --dumpData
  #       done
  #     done
  #   done
  # fi

  # Single2D data
  mkdir -p "${data_dir}/Single2D"
  mkdir -p "${ref_dir}/Single2D"
  # for db in 103;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 104;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 6 8;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 106;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 15 16;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # for db in 108;
  # do
  #   for st in 0 1 2 3;
  #   do
  #     for np in 4 32 144 256;
  #     do
  #       ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --dumpData
  #     done
  #   done
  # done
  # # Generate all ranks
  # if [ "$genAllRanks" = true ];
  # then
  #   for db in 108;
  #   do
  #     for st in 0 1 2 3;
  #     do
  #       for np in 64;
  #       do
  #         ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --id ${np} --np ${np} --db ${db} --stage ${st} --dumpData
  #       done
  #     done
  #   done
  # fi

  # Tubular data
  mkdir -p "${data_dir}/Tubular"
  mkdir -p "${ref_dir}/Tubular"
  for db in 103;
  do
    for st in 0 1 2 3;
    do
      for np in 4;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --dumpData
      done
    done
  done
#   for db in 104;
#   do
#     for st in 0 1 2 3;
#     do
#       for np in 4 24 56;
#       do
#         ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --dumpData
#       done
#     done
#   done
#   for db in 106;
#   do
#     for st in 0 1 2 3;
#     do
#       for np in 6 128 288;
#       do
#         ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --dumpData
#       done
#     done
#   done
#   for db in 108;
#   do
#     for st in 0 1 2 3;
#     do
#       for np in 4 288 512;
#       do
#         ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --dumpData
#       done
#     done
#   done
#   # Generate all ranks
#   if [ "$genAllRanks" = true ];
#   then
#     for db in 108;
#     do
#       for st in 0 1 2 3;
#       do
#         for np in 144;
#         do
#           ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --id ${np} --np ${np} --db ${db} --stage ${st} --dumpData
#         done
#       done
#     done
#   fi
done
fi

# triangular truncation
if [ "$genTriangular" = true ];
then
  for sch in "WLFl" "WLFm";
  do
    data_dir="_data/Framework/LoadSplitter/${sch}/"
    ref_dir="_refdata/Framework/LoadSplitter/${sch}/"

  # Create path if it doesn't exist
  mkdir -p "${data_dir}"
  mkdir -p "${ref_dir}"

  # Serial data
  mkdir -p "${data_dir}"/Serial
  mkdir -p "${ref_dir}"/Serial
  for db in 104 106 108;
  do
    for st in 0 3;
    do
      ./FrameworkLoadSplitterTests \[${sch}\] --algorithm serial --np 1 --db ${db} --stage ${st} --truncation triangular --dumpData
    done
  done

  # Single1D data
  mkdir -p "${data_dir}/Single1D"
  mkdir -p "${ref_dir}/Single1D"
  for db in 104;
  do
    for st in 0 3;
    do
      for np in 4 15 16;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 106;
  do
    for st in 0 3;
    do
      for np in 4 16 42 64;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 108;
  do
    for st in 0 3;
    do
      for np in 4 32 144 256;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  # Generate all ranks
  if [ "$genAllRanks" = true ];
  then
    for db in 108;
    do
      for st in 0  3;
      do
        for np in 56;
        do
          ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single1d --id ${np} --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
        done
      done
    done
  fi

  # Single2D data
  mkdir -p "${data_dir}/Single2D"
  mkdir -p "${ref_dir}/Single2D"
  for db in 104;
  do
    for st in 0 3;
    do
      for np in 4 6 8;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 106;
  do
    for st in 0 3;
    do
      for np in 4 15 16;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 108;
  do
    for st in 0 3;
    do
      for np in 4 32 144 256;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  # Generate all ranks
  if [ "$genAllRanks" = true ];
  then
    for db in 108;
    do
      for st in 0 3;
      do
        for np in 64;
        do
          ./FrameworkLoadSplitterTests \[${sch}\] --algorithm single2d --id ${np} --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
        done
      done
    done
  fi

  # Tubular data
  mkdir -p "${data_dir}/Tubular"
  mkdir -p "${ref_dir}/Tubular"
  for db in 104;
  do
    for st in 0 3;
    do
      for np in 4 24 56;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 106;
  do
    for st in 0 3;
    do
      for np in 6 128 288;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  for db in 108;
  do
    for st in 0 3;
    do
      for np in 4 288 512;
      do
        ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
      done
    done
  done
  # Generate all ranks
  if [ "$genAllRanks" = true ];
  then
    for db in 108;
    do
      for st in 0 3;
      do
        for np in 144;
        do
          ./FrameworkLoadSplitterTests \[${sch}\] --algorithm tubular --id ${np} --np ${np} --db ${db} --stage ${st} --truncation triangular --dumpData
        done
      done
    done
  fi
done
fi

if [ "$copy2Ref" = true ];
then
  cp -r "./${data_dir}"/* "./${ref_dir}"/
fi
