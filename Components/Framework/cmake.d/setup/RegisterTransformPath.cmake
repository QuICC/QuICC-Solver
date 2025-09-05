set(tags
  NegCurlCurlNl
  CurlCurlNl
  CurlNl
  Empty
  I2CurlCurlNl
  I2CurlNl
  I2LaplhCurlNl
  I2ScalarNl
  I2LaplhCurlCurlNl
  NegI2CurlCurlNl
  NegI2rCurlCurlNl
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
