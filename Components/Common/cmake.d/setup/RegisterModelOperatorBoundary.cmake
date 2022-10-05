set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "ModelOperatorBoundary"
  BASECLASS "IModelOperatorBoundary"
  EXCLUDED ${excluded}
  )
