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
  ValueWithL0CurlNl
  InsulatingBc2NegCurlCurlNl
  InsulatingScalar
  InsulatingScalarNl
  InsulatingTorPol
  InsulatingCurlNl
  InsulatingCurlCurlNl
  InsulatingNegCurlCurlNl
  InsulatingLaplhCurlNl
  InsulatingLaplhCurlCurlNl
  InsulatingWithL0CurlNl
  StressFreeTorPol
  StressFreeCurlNl
  NoSlipTorPol
  NoSlipBc1NegCurlCurlNl
  NoPenetrationTorPol
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Path"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "path_"
  )
