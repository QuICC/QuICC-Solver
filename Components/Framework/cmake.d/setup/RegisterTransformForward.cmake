set(tags 
  D1
  D1ZP0
  D2
  DfOverlaplh
  I2D1
  I2P
  I2Q
  I2S
  I2T
  I2rQ
  I2rS
  I2ZI2D1
  I4D1
  I4D1ZI2
  I4P
  I4Q
  I4Overr1D1R1
  I4Overr1D1R1ZI2
  I4Overr1Pm
  I4S
  I6Laplh
  I6LaplhZI4D1R1
  I6Overr1D1R1
  I6Overr1D1R1ZI4
  I6Overr1Pm
  Laplh
  Laplh2
  LaplhD1
  LaplhOversin
  LaplhOversinDphi
  Overlaplh
  OverlaplhD1
  OverlaplhOversin
  OverlaplhOversinDphi
  P
  P0
  Pm
  Pol
  Q
  R1
  S
  Oversin
  OversinDphi
  T
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Forward"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "fwd_"
  )
