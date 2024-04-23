set(tags
  NegCurlCurlNl
  CurlCurlNl
  CurlNl
  Empty
  I2CurlCurlNl
  I2CurlNl
  I2ScalarNl
  NegI2CurlCurlNl
  NegI4CurlCurlNl
  Scalar
  ScalarNl
  TorPol
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Path"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "path_"
  )

set(value_tags
  Scalar
  ScalarNl
  TorPol
  Tor
  Pol
  CurlNl
  CurlCurlNl
  NegCurlCurlNl
  )

quicc_register_tags(
  NAMESPACE "Transform/Path/Value"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${value_tags}
  PREFIX "path_value_"
  )

set(insulating_tags
  Scalar
  ScalarNl
  TorPol
  Tor
  Pol
  CurlNl
  CurlCurlNl
  NegCurlCurlNl
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Path/Insulating"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${insulating_tags}
  PREFIX "path_insulatin_"
  )
