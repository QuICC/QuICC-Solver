set(tags
  Diagnostic
  Prognostic
  Trivial
  Wrapper
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "PseudospectralTag"
  BASECLASS "IPseudospectralTag"
  TAGS ${tags}
  )
