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
