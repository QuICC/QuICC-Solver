set(tags
  Diagnostic
  Prognostic
  Trivial
  Uninitialized
  Wrapper
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "PseudospectralTag"
  BASECLASS "IPseudospectralTag"
  TAGS ${tags}
  )
