set(tags
  D1
  D1Laplh
  D1Laplhm1
  D1LaplhZD1Overr1D1R1
  D1R1
  D1Overr1
  D1ZP
  D2
  D3
  DfLaplh
  DsLaplh
  Laplh
  LaplhD0
  LaplhZOverr1D1R1
  P
  P0
  Overr1
  Overr1D0
  Overr1D1
  Overr1D1R1
  Overr2D1R1
  OverrSqD1R1
  Overr1D2R1
  Overr1LaplhPm
  Overr1Pm
  Overr2
  OverrSq
  Slapl
  SlaplD0
  Slaplr
  RGradRTh
  RGradRThD0
  Oversin
  OversinD1Sin
  D1OversinDphi
  OversinDphi
  OversinLaplh
  OversinLaplhDphi
  OversinLaplhm1Dphi
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Backward"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "bwd_"
  )
