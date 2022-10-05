set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Backward"
  BASECLASS "IOperator"
  COMMON_DIR "../Common"
  EXCLUDED ${excluded}
  )
