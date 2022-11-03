set(tags
  Codensity
  Density
  DxMeanTemperature
  DzMeanTemperature
  Entropy
  FluctMagnetic
  FluctMagneticX
  FluctMagneticY
  FluctMagneticZ
  FluctTemperature
  FluctVelocity
  FluctVelocityX
  FluctVelocityY
  FluctVelocityZ
  Magnetic
  MagneticX
  MagneticY
  MagneticZ
  MeanMagnetic
  MeanMagneticX
  MeanMagneticY
  MeanMagneticZ
  MeanTemperature
  MeanVelocity
  MeanVelocityX
  MeanVelocityY
  MeanVelocityZ
  Pressure
  Streamfunction
  Temperature
  Undefined
  Velocity
  VelocityX
  VelocityY
  VelocityZ
  Vorticity
  VorticityX
  VorticityY
  VorticityZ
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "PhysicalNames"
  BASECLASS "IPhysicalName"
  TAGS ${tags}
  )
