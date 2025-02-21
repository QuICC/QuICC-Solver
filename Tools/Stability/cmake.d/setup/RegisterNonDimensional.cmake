set(tags
  Nev
  StabilityMode
  Sort
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
