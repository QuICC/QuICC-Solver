set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "SolveTiming"
  BASECLASS "ISolveTiming"
  EXCLUDED ${excluded}
  )
