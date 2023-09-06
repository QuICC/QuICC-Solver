set(tags
  Default
  GaussianQuadrature
  Fft
  Uniform
  Triangular
  Trapezoidal
  )

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "Transform/Setup"
  BASECLASS "IOption"
  COMMON_DIR "../Common"
  TAGS ${tags}
  )
