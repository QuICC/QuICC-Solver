set(tags
  Error
  Explicit
  Implicit
  Intermediate
  Rhs
  Solution
  Temporary
  Influence
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Register"
  BASECLASS "IDataRegister"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
