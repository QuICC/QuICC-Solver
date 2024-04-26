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
  ValueScalar
  ValueScalarNl
  ValueTorPol
  ValueCurlNl
  ValueCurlCurlNl
  ValueNegCurlCurlNl
  ValueBc1NegCurlCurlNl
  ValueLaplhCurlNl
  ValueLaplhCurlCurlNl
  InsulatingScalar
  InsulatingScalarNl
  InsulatingTorPol
  InsulatingCurlNl
  InsulatingCurlCurlNl
  InsulatingNegCurlCurlNl
  InsulatingLaplhCurlNl
  InsulatingLaplhCurlCurlNl
  StressFreeTorPol
  StressFreeCurlNl
  NoSlipTorPol
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Path"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "path_"
  )
