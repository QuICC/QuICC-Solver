set(tags 
  NoSlip 
  StressFree
  FixedTemperature
  FixedFlux
  Insulating
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Bc/Name"
  BASECLASS "IName"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
