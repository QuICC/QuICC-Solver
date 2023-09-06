set(tags
  Boundary
  ExplicitLinear
  ExplicitNextstep
  ExplicitNonlinear
  ImplicitLinear
  Stencil
  Time
  SplitBoundary
  SplitBoundaryValue
  SplitImplicitLinear
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "ModelOperator"
  BASECLASS "IModelOperator"
  TAGS ${tags}
  )
