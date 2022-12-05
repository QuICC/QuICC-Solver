set(tags
  Error
  Explicit
  Implicit
  Intermediate
  Rhs
  Solution
  Temporary
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Register"
  BASECLASS "IDataRegister"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
