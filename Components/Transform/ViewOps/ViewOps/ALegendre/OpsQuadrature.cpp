// External includes
//

// Project includes
//
#include "ViewOps/ALegendre/OpsQuadrature.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {
namespace details {


void ALegendreRule::computeQuadrature(Internal::Array& igrid,
   Internal::Array& iweights, const int gSize)
{
   ::QuICC::Polynomial::Quadrature::LegendreRule quad;
   quad.computeQuadrature(igrid, iweights, gSize);
   // scale for spherical harmonics
   iweights.array() *= 2.0 * ::QuICC::Internal::Math::PI;
}

} // namespace details


} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
