set(tags 
  GoOn
  Stop
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "RuntimeStatus"
  BASECLASS "IRuntimeStatus"
  TAGS ${tags}
  )
