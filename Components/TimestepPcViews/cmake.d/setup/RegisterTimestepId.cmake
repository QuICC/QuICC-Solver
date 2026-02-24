set(tags
  ImexPc2
  ImexEuler
  ImexPc2b
  ImexEulerb
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Timestep/Id"
  BASECLASS "ITimeScheme"
  COMMON_DIR "../Common"
  TAGS ${tags}
  REGISTRATOR registerAllPc
  )
