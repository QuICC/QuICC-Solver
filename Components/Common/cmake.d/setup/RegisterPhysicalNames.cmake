set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "PhysicalNames"
  BASECLASS "IPhysicalName"
  EXCLUDED ${excluded}
  )
