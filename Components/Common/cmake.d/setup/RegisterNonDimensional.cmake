set(excluded 
  "Typedefs.hpp"
  )

set(tags
  Alpha
  Beta
  CflAlfvenDamping
  CflAlfvenScale
  CflInertial
  CflTorsional
  Chandrasekhar
  Chi
  Delta
  Eady
  Ekman
  Elevator
  Elsasser
  Epsilon
  Eta
  FastMean
  Gamma
  Heating
  Iota
  Kappa
  Lambda
  Lower1d
  Lower2d
  Lower3d
  MagneticEkman
  MagneticPrandtl
  MagneticReynolds
  ModifiedElsasser
  Mu
  Nu
  Omega
  Omicron
  Phi
  Pi
  Poincare
  Prandtl
  Psi
  RRatio
  Rayleigh
  Rescaled
  Rho
  Roberts
  Rossby
  Sigma
  Tau
  Taylor
  Theta
  Upper1d
  Upper2d
  Upper3d
  Upsilon
  Xi
  Zeta
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "NonDimensional"
  BASECLASS "INumber"
  TAGS ${tags}
  EXCLUDED ${excluded}
  VALUE "value"
  )
