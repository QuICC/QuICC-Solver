set(tags 
  Add
  None
  Set
  SetNeg
  Sub
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Arithmetics"
  BASECLASS "IArithmetics"
  TAGS ${tags}
  )
