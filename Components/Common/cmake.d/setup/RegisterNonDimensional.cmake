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
  GrowthRate
  Heating
  Iota
  Kappa
  Lambda
  Lehnert
  Lower1d
  Lower2d
  Lower3d
  Lundquist
  MagneticEkman
  MagneticPrandtl
  MagneticReynolds
  MaxIteration
  ModifiedElsasser
  Mu
  Nev
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
  StabilityMode
  Sort
  Tau
  Taylor
  Theta
  Tolerance
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
