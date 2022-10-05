set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "NonDimensional"
  BASECLASS "INumber"
  EXCLUDED ${excluded}
  VALUE "value"
  )
