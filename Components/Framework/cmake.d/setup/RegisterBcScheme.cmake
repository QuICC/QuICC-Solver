set(tags 
  Galerkin
  Tau
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Bc/Scheme"
  BASECLASS "IScheme"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
