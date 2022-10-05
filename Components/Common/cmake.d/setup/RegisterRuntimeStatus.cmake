set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "RuntimeStatus"
  BASECLASS "IRuntimeStatus"
  EXCLUDED ${excluded}
  )
