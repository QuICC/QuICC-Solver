set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Forward"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  EXCLUDED ${excluded}
  )
