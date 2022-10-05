set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Arithmetics"
  BASECLASS "IArithmetics"
  EXCLUDED ${excluded}
  )
