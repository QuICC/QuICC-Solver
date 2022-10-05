set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "ModelOperator"
  BASECLASS "IModelOperator"
  EXCLUDED ${excluded}
  )
