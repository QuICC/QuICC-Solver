set(tags
  ImexPc2
  ImexPc2b
  ImexEuler
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Timestep/Id"
  BASECLASS "ITimeScheme"
  COMMON_DIR "../Common"
  TAGS ${tags}
  REGISTRATOR registerAllPc
  )
