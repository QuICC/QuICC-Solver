set(tags 
  D1
  D1Laplh
  D1LaplhZD1Overr1D1R1
  D1R1
  D1ZP
  D2
  D3
  DfLaplh
  DsLaplh
  Laplh
  LaplhZOverr1D1R1
  P
  P0
  Overr1
  Overr1D1R1
  Overr1LaplhPm
  Overr1Pm
  Overr2
  Slapl
  Slaplr
  Oversin
  OversinD1Sin
  OversinDphi
  OversinLaplh
  OversinLaplhDphi
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Backward"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "bwd_"
  )
