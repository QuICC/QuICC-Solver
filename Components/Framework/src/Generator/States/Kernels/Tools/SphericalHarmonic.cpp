/**
 * @file SphericalHarmonic.cpp
 * @brief Source of the implementation of the tools to generate exact spherical harmonics in physical space
 */

// Configuration includes
//

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Tools/SphericalHarmonic.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Constants.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Tools {

   void SphericalHarmonic::Ylm(Array& rField, const int l, const int m, const MHDComplex c, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = c.real();
      MHDFloat im = c.imag();

      rField = 2.0*re*(static_cast<MHDFloat>(m)*phi).array().cos() - 2.0*im*(static_cast<MHDFloat>(m)*phi).array().sin();
      #if defined QUICC_SH_NORM_SCHMIDT
         // Schmidt quasi-normalized spherical harmonic Y_l^m
         MHDFloat leg = std::sph_legendre(l, m, theta)*std::sqrt(4.0*Math::PI/static_cast<MHDFloat>(2*l + 1));
      #elif defined QUICC_SH_NORM_UNITY
         // Normalized spherical harmonic Y_l^m
         MHDFloat leg = std::sph_legendre(l, m, theta);
      #endif //defined QUICC_SH_NORM_SCHMIDT

      rField *= leg;
   }

   void SphericalHarmonic::Torlm(Array& rField, const int l, const int m, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      if(l == 0 && m == 0)
      {
         SphericalHarmonic::Tor00(rField, id, ct, theta, phi);
      } else if(l == 1 && m == 0)
      {
         SphericalHarmonic::Tor10(rField, id, ct, theta, phi);
      } else if(l == 1 && m == 1)
      {
         SphericalHarmonic::Tor11(rField, id, ct, theta, phi);
      } else if(l == 2 && m == 0)
      {
         SphericalHarmonic::Tor20(rField, id, ct, theta, phi);
      } else if(l == 2 && m == 1)
      {
         SphericalHarmonic::Tor21(rField, id, ct, theta, phi);
      } else if(l == 2 && m == 2)
      {
         SphericalHarmonic::Tor22(rField, id, ct, theta, phi);
      } else if(l == 5 && m == 4)
      {
         SphericalHarmonic::Tor54(rField, id, ct, theta, phi);
      } else
      {
         throw std::logic_error("Requested unimplemented toroidal harmonic");
      }
   }

   void SphericalHarmonic::Pollm(Array& rField, const int l, const int m, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      if(l == 0 && m == 0)
      {
         SphericalHarmonic::Pol00(rField, id, cq, cs, theta, phi);
      } else if(l == 1 && m == 0)
      {
         SphericalHarmonic::Pol10(rField, id, cq, cs, theta, phi);
      } else if(l == 1 && m == 1)
      {
         SphericalHarmonic::Pol11(rField, id, cq, cs, theta, phi);
      } else if(l == 2 && m == 0)
      {
         SphericalHarmonic::Pol20(rField, id, cq, cs, theta, phi);
      } else if(l == 2 && m == 1)
      {
         SphericalHarmonic::Pol21(rField, id, cq, cs, theta, phi);
      } else if(l == 2 && m == 2)
      {
         SphericalHarmonic::Pol22(rField, id, cq, cs, theta, phi);
      } else if(l == 4 && m == 3)
      {
         SphericalHarmonic::Pol43(rField, id, cq, cs, theta, phi);
      } else
      {
         throw std::logic_error("Requested unimplemented poloidal harmonic");
      }
   }

   void SphericalHarmonic::Tor00(Array& rField, FieldComponents::Physical::Id id, const MHDComplex, const MHDFloat, const Array& phi)
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor10(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      //MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = re*std::sqrt(3.0/Math::PI)*std::sin(theta);
         rField = Array::Ones(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor11(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = std::sqrt(3.0/(2.0*Math::PI));
         rField = im*phi.array().cos() + re*phi.array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = std::sqrt(3.0/(2.0*Math::PI))*std::cos(theta);
         rField = re*phi.array().cos() - im*phi.array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor20(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      //MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = re*3.0*std::sqrt(5.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = Array::Ones(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor21(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta);
         rField = im*phi.array().cos() + re*phi.array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = std::sqrt(15.0/(2.0*Math::PI))*std::cos(2.0*theta);
         rField = re*phi.array().cos() - im*phi.array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor22(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -std::sqrt(15.0/(2.0*Math::PI))*std::sin(theta);
         rField = im*(2.0*phi).array().cos() + re*(2.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = -std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta)*std::sin(theta);
         rField = re*(2.0*phi).array().cos() - im*(2.0*phi).array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Tor54(Array& rField, FieldComponents::Physical::Id id, const MHDComplex ct, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re = ct.real();
      MHDFloat im = ct.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -(3.0/2.0)*std::sqrt(385.0/(2.0*Math::PI))*std::cos(theta)*std::pow(std::sin(theta),3);
         rField = im*(4.0*phi).array().cos() + re*(4.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = -(3.0/16.0)*std::sqrt(385.0/(2.0*Math::PI))*(3.0 + 5.0*std::cos(2.0*theta))*std::pow(std::sin(theta),3);
         rField = re*(4.0*phi).array().cos() - im*(4.0*phi).array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol00(Array& rField, FieldComponents::Physical::Id id, const MHDComplex, const MHDComplex, const MHDFloat, const Array& phi)
   {
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol10(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      //MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      //MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = re_q*2.0*std::sqrt(3.0/Math::PI)*std::cos(theta);
         rField = Array::Ones(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -re_s*std::sqrt(3.0/Math::PI)*std::sin(theta);
         rField = Array::Ones(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol11(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = -std::sqrt(6.0/Math::PI)*std::sin(theta);
         rField = re_q*phi.array().cos() - im_q*phi.array().sin();

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -std::sqrt(3.0/(2.0*Math::PI))*std::cos(theta);
         rField = re_s*phi.array().cos() - im_s*phi.array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = std::sqrt(3.0/(2.0*Math::PI));
         rField = im_s*phi.array().cos() + re_s*phi.array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol20(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      //MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      //MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = re_q*(3.0/2.0)*std::sqrt(5.0/Math::PI)*(1.0 + 3.0*std::cos(2.0*theta));
         rField = Array::Ones(phi.size());

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -re_s*3.0*std::sqrt(5.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = Array::Ones(phi.size());

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol21(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = -3.0*std::sqrt(30.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = re_q*phi.array().cos() - im_q*phi.array().sin();

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -std::sqrt(15.0/(2.0*Math::PI))*std::cos(2.0*theta);
         rField = re_s*phi.array().cos() - im_s*phi.array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta);
         rField = im_s*phi.array().cos() + re_s*phi.array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol22(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = 3.0*std::sqrt(15.0/(2.0*Math::PI))*std::pow(std::sin(theta),2);
         rField = re_q*(2.0*phi).array().cos() - im_q*(2.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta)*std::sin(theta);
         rField = re_s*(2.0*phi).array().cos() - im_s*(2.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = -std::sqrt(15.0/(2.0*Math::PI))*std::sin(theta);
         rField = im_s*(2.0*phi).array().cos() + re_s*(2.0*phi).array().sin();
      }
      rField *= factor;
   }

   void SphericalHarmonic::Pol43(Array& rField, FieldComponents::Physical::Id id, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi)
   {
      MHDFloat re_q = cq.real();
      MHDFloat im_q = cq.imag();
      MHDFloat re_s = cs.real();
      MHDFloat im_s = cs.imag();
      MHDFloat factor = 1.0;
      if(id == FieldComponents::Physical::R)
      {
         factor = -15.0*std::sqrt(35.0/Math::PI)*std::cos(theta)*std::pow(std::sin(theta),3);
         rField = re_q*(3.0*phi).array().cos() - im_q*(3.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::THETA)
      {
         factor = -(3.0/4.0)*std::sqrt(35.0/Math::PI)*(1.0 + 2.0*std::cos(2.0*theta))*std::pow(std::sin(theta),2);
         rField = re_s*(3.0*phi).array().cos() - im_s*(3.0*phi).array().sin();

      } else if(id == FieldComponents::Physical::PHI)
      {
         factor = (9.0/4.0)*std::sqrt(35.0/Math::PI)*std::cos(theta)*std::pow(std::sin(theta),2);
         rField = im_s*(3.0*phi).array().cos() + re_s*(3.0*phi).array().sin();
      }
      rField *= factor;
   }

}
}
}
}
