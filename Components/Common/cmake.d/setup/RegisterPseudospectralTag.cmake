set(excluded 
  "Typedefs.hpp"
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "PseudospectralTag"
  BASECLASS "IPseudospectralTag"
  EXCLUDED ${excluded}
  )
