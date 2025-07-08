set(tags 
  Energy
  EnergyD1
  EnergyD1R1
  EnergyR2
  EnergySlaplR2
  Power
  PowerD1
  PowerD1R1
  PowerR2
  PowerSlaplR2
  RadialPower
  RadialPowerOverr1
  RadialPowerOverr1D1R1
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Reductor"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "red_"
  )
