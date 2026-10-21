set(tags 
  NoPenetration
  NoSlip 
  StressFree
  FixedTemperature
  FixedEntropy
  FixedFlux
  Insulating
  InsulatingTopPerfcondBottom
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Bc/Name"
  BASECLASS "IName"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
