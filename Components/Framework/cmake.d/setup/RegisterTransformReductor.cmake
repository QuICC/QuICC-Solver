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
  ValueEnergy
  ValueEnergyD1
  ValueEnergyD1R1
  ValueEnergyR2
  ValueEnergySlaplR2
  ValuePower
  ValuePowerD1
  ValuePowerD1R1
  ValuePowerR2
  ValuePowerSlaplR2
  InsulatingEnergy
  InsulatingEnergyD1
  InsulatingEnergyD1R1
  InsulatingEnergyR2
  InsulatingEnergySlaplR2
  InsulatingPower
  InsulatingPowerD1
  InsulatingPowerD1R1
  InsulatingPowerR2
  InsulatingPowerSlaplR2
  NoSlipEnergy
  NoSlipEnergyD1
  NoSlipEnergyD1R1
  NoSlipEnergyR2
  NoSlipEnergySlaplR2
  NoSlipPower
  NoSlipPowerD1
  NoSlipPowerD1R1
  NoSlipPowerR2
  NoSlipPowerSlaplR2
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Reductor"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  TAGS ${tags}
  PREFIX "red_"
  )
