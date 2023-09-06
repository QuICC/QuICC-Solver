set(tags 
  After
  Before
  Prognostic
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "SolveTiming"
  BASECLASS "ISolveTiming"
  TAGS ${tags}
  )
