set(tags
  GrowthRate
  MaxIteration
  Nev
  Tolerance
  StabilityMode
  Sort
  WriteMtx
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "NonDimensional"
  BASECLASS "INumber"
  COMMON_DIR ../../Components/Common
  TAGS ${tags}
  VALUE "value"
  REGISTRATOR registerStability
  )
