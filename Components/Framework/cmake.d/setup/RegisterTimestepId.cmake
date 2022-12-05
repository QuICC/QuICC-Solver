set(tags 
  ImexRkCb2
  ImexRkCb3b
  ImexRkCb3c
  ImexRkCb3d
  ImexRkCb3e
  ImexRkCb3f
  ImexRkCb4
  ImexPc2
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Timestep/Id"
  BASECLASS "ITimeScheme"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
