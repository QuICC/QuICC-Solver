set(tags 
  On
  Off
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Tag/Generic"
  BASECLASS "IGeneric"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
