set(tags
  Lhs
  Influence
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Tag/Operator"
  BASECLASS "ISolverMatrix"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
