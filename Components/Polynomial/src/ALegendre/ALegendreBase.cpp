/**
 * @file ALegendreBase.cpp
 * @brief Source of the implementation of the associated Legendre polynomial
 * 
 * Implements recurrence relations to calculate Plm and derived functions, such as d Plm / d\theta
 * 
 * notes:
 * 
 * ** Normalisation factors, cs **
 *   they are chosen so that, assuming:
 *          Ylm = clm Plm exp(im\phi), clm = ||Ylm||  sqrt( (2l+1) (l-m)! / 2 / (l+m)! ) / sqrt(2 Pi)
 *   a recurrence relation such as
 *          (l-m) Plm = (2l-1) x Pl-1m - (l+m-1) Pl-2m
 *   is coded like this:
 *          iplm =  (||Ylm|| / clm) ^(-1) [ (||Yl_1m|| / cl_1m) (2l-1) x ipl_1m - (||Yl_2m|| / cl_2m) (l+m-1) ipl_2m ]
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/ALegendre/ALegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   ALegendreBase::ALegendreBase()
   {
   }

   ALegendreBase::~ALegendreBase()
   {
   }

   void ALegendreBase::Plm(Eigen::Ref<Internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& ipl_1m, const Eigen::Ref<const Internal::Matrix>& ipl_2m, const Internal::Array& igrid, NormalizerML norm)
   {
      // Safety assert
      assert(l-m > 0);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::MHDFloat dm = Internal::MHDFloat(m);
      Internal::Array cs = norm(dm, dl);

      // (l-m) Plm = (2l-1) x Pl-1m - (l+m-1) Pl-2m

      iplm.array() = cs(1)*(cs(0)*ipl_2m.array() + (igrid.array()*ipl_1m.array()).array());
   }

   void ALegendreBase::Pmm(Eigen::Ref<Internal::Matrix> ipmm, const int m, const Internal::Array& igrid, NormalizerM norm)
   {
      Internal::MHDFloat dm = Internal::MHDFloat(m);
      Internal::Array cs = norm(dm);

      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         ipmm.setConstant(cs(0));
      } else
      {
         Internal::MHDFloat di = MHD_MP(1.0);

         for(int i = 1; i <= m; i++)
         {
            cs(0) *= -Internal::Math::sqrt(di/(di + MHD_MP(1.0)));
            di += MHD_MP(2.0);
         }
         ipmm.array() = cs(0)*(MHD_MP(1.0) - igrid.array().pow(2)).pow(dm/2);
      }
   }

   void ALegendreBase::Pm1m(Eigen::Ref<Internal::Matrix> ipm1m, const int m, const Eigen::Ref<const Internal::Matrix>& ipmm, const Internal::Array& igrid, NormalizerM norm)
   {
      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         Internal::MHDFloat dm = Internal::MHDFloat(m);
         Internal::Array cs = norm(dm);

         ipm1m.array() = (cs(0)*igrid).array()*ipmm.array();
      }
   }

   void ALegendreBase::dPl0(Eigen::Ref<Internal::Matrix> idpl0, const int l, const Eigen::Ref<const Internal::Matrix>& ipl1, NormalizerL norm)
   {
      // Safety assert
      assert(l > 0);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::Array cs = norm(dl);

      idpl0.array() = cs(0)*ipl1.array();
   }

   void ALegendreBase::d2P10(Eigen::Ref<Internal::Matrix> id2p10, const Eigen::Ref<const Internal::Matrix>& ip10, Normalizer norm)
   {
      Internal::Array cs = norm();

      id2p10.array() = cs(0)*ip10.array();
   }


   void ALegendreBase::d2Pl0(Eigen::Ref<Internal::Matrix> id2pl0, const int l, const Eigen::Ref<const Internal::Matrix>& ipl2, const Eigen::Ref<const Internal::Matrix>& ipl0, NormalizerL norm)
   {
      // Safety assert
      assert(l > 0);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::Array cs = norm(dl);

      id2pl0.array() = -cs(0)*ipl0.array() + cs(1)*ipl2.array();
   }

   void ALegendreBase::d2Pl1(Eigen::Ref<Internal::Matrix> id2pl1, const int l, const Eigen::Ref<const Internal::Matrix>& ipl3, const Eigen::Ref<const Internal::Matrix>& ipl1, NormalizerL norm)
   {
      // Safety assert
      assert(l > 0);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::Array cs = norm(dl);

      id2pl1.array() = -cs(0)*ipl1.array() + cs(1)*ipl3.array();
   }


   void ALegendreBase::dPmm(Eigen::Ref<Internal::Array> idpmm, const int m, const Eigen::Ref<const Internal::Array>& iplm_1, NormalizerM norm)
   {
      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         idpmm.setConstant(MHD_MP(0.0));

      } else
      {
         Internal::MHDFloat dm = Internal::MHDFloat(m);
         Internal::Array cs = norm(dm);

         idpmm = cs(0)*iplm_1;
      }
   }

   void ALegendreBase::dPlm(Eigen::Ref<Internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& iplm_1, const Eigen::Ref<const Internal::Matrix>& iplm1, NormalizerML norm)
   {
      // Safety assert
      assert(l-m > 0);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::MHDFloat dm = Internal::MHDFloat(m);
      Internal::Array cs = norm(dm, dl);

      // dPlm/dtheta = -sqrt(1-x^2) dPlm/dx = -(1/2)[ (l+m)(l-m+1)Plm-1 - Plm+1 ]

      idplm.array() = cs(0)*iplm_1.array() - cs(1)*iplm1.array();
   }

   void ALegendreBase::d2Pmm0m1(Eigen::Ref<Internal::Array> op, const int m, const int l, const Eigen::Ref<const Internal::Array>& iplm_2, const Eigen::Ref<const Internal::Array>& iplm, NormalizerML norm)
   {
      // works for both l=m and l=m+1
      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op.setConstant(MHD_MP(0.0));

      } else
      {
         Internal::MHDFloat dl = Internal::MHDFloat(l);
         Internal::MHDFloat dm = Internal::MHDFloat(m);
         Internal::Array cs = norm(dm, dl);

         op = cs(0)*iplm_2 - cs(1)*iplm;
      }
   }

   void ALegendreBase::d2Plm(Eigen::Ref<Internal::Matrix> id2plm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& iplm_2, const Eigen::Ref<const Internal::Matrix>& iplm, const Eigen::Ref<const Internal::Matrix>& iplm2, NormalizerML norm)
   {
      // Safety assert
      assert(l-m > 1);

      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::MHDFloat dm = Internal::MHDFloat(m);
      Internal::Array cs = norm(dm, dl);

      // d2Plm/dtheta2 = (1/4) A1(l,m)A1(l,m-1) Plm-2
      //                 -(1/4) [A1(l,m) + A1(l,m+1)] Plm
      //                 +(1/4) Plm+2

      id2plm.array() = cs(0)*iplm_2.array() - cs(1)*iplm.array() + cs(2)*iplm2.array();
   }

   void ALegendreBase::sin_1Plm(Eigen::Ref<Internal::Matrix> isin_1plm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& ipl1m1, const Eigen::Ref<const Internal::Matrix>& ipl1m_1, NormalizerML norm)
   {
      Internal::MHDFloat dl = Internal::MHDFloat(l);
      Internal::MHDFloat dm = Internal::MHDFloat(m);
      Internal::Array cs = norm(dm, dl);

      isin_1plm.array() = cs(1)*(ipl1m1.array() + cs(0)*ipl1m_1.array());
   }

   //
   // General polynomial normalizer
   //
   ALegendreBase::NormalizerM ALegendreBase::normPmm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtPmm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitPmm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerM ALegendreBase::normPm1m()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtPm1m;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitPm1m;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerML ALegendreBase::normPlm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtPlm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitPlm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerM ALegendreBase::normdPmm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtdPmm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitdPmm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerL ALegendreBase::normdPl0()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtdPl0;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitdPl0;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::Normalizer ALegendreBase::normd2P10()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtd2P10;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitd2P10;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerL ALegendreBase::normd2Pl0()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtd2Pl0;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitd2Pl0;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerL ALegendreBase::normd2Pl1()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtd2Pl1;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitd2Pl1;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerML ALegendreBase::normdPlm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtdPlm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitdPlm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerML ALegendreBase::normd2Plm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtd2Plm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitd2Plm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   ALegendreBase::NormalizerML ALegendreBase::normsin_1Plm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtsin_1Plm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitsin_1Plm;
      #endif //defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
   }

   Internal::MHDFloat ALegendreBase::A1(const Internal::MHDFloat dl, const Internal::MHDFloat dm)
   {
      // (l,m)-dependent coefficient for the relation:
      // dPlm/dtheta = (1/2) [ A1(l,m) Plm-1 - Plm+1 ]

      return (dl+dm)*(dl-dm+1);
   }

   Internal::Array ALegendreBase::unitPmm(const Internal::MHDFloat dm)
   {
      Internal::Array cs(1);

      if(dm == MHD_MP(0.0))
      {
         cs(0) = Internal::Math::sqrt(MHD_MP(1.0)/(MHD_MP(4.0)*Internal::Math::PI));
      } else
      {
         cs(0) = Internal::Math::sqrt((MHD_MP(2.0)*dm + MHD_MP(1.0))/(MHD_MP(4.0)*Internal::Math::PI));
      }

      return cs;
   }

   Internal::Array ALegendreBase::unitPm1m(const Internal::MHDFloat dm)
   {
      Internal::Array cs(1);

      cs(0) = Internal::Math::sqrt(MHD_MP(2.0)*dm + MHD_MP(3.0));

      return cs;
   }

   Internal::Array ALegendreBase::unitPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = -Internal::Math::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)));
      cs(1) = Internal::Math::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));

      return cs;
   }

   Internal::Array ALegendreBase::unitdPl0(const Internal::MHDFloat dl)
   {
      Internal::Array cs(1);

      cs(0) = Internal::Math::sqrt(dl*(dl + MHD_MP(1.0)));

      return cs;
   }

   Internal::Array ALegendreBase::unitd2P10()
   {
      Internal::Array cs(1);

      cs(0) = -MHD_MP(0.25)* ( ALegendreBase::A1(MHD_MP(1.0), MHD_MP(0.0)) + ALegendreBase::A1(MHD_MP(1.0), MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::unitd2Pl0(const Internal::MHDFloat dl)
   {
      Internal::MHDFloat dm = MHD_MP(0.0);

      Internal::Array cs(2);

      cs(0) = MHD_MP(0.25) * ( ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm + MHD_MP(1.0)) );

      cs(1) = MHD_MP(0.25) * ( MHD_MP(1.0) 
                              + ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm - MHD_MP(1.0)) 
                                 / ( (dl + MHD_MP(2.0)) * (dl + MHD_MP(1.0)) * dl * (dl - MHD_MP(1.0)) ) 
                              ) * Internal::Math::sqrt( (dl + MHD_MP(2.0)) * (dl + MHD_MP(1.0)) * dl * (dl - MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::unitd2Pl1(const Internal::MHDFloat dl)
   {
      Internal::MHDFloat dm = MHD_MP(1.0);

      Internal::Array cs(2);

      cs(0) = MHD_MP(0.25) * (ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm - MHD_MP(1.0)) 
                                 / ( (dl + MHD_MP(1.0)) * dl ) 
                              + ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm + MHD_MP(1.0)) );

      cs(1) = MHD_MP(0.25) *  Internal::Math::sqrt( (dl + MHD_MP(3.0)) *(dl + MHD_MP(2.0)) * (dl - MHD_MP(1.0)) * (dl - MHD_MP(2.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::unitdPmm(const Internal::MHDFloat dm)
   {
      Internal::Array cs(1);

      cs(0) = -Internal::Math::sqrt(dm/MHD_MP(2.0));

      return cs;
   }

   Internal::Array ALegendreBase::unitdPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = -MHD_MP(0.5)*Internal::Math::sqrt((dl - dm + MHD_MP(1.0))*(dl + dm));

      cs(1) = -MHD_MP(0.5)*Internal::Math::sqrt((dl - dm)*(dl + dm + MHD_MP(1.0)));

      return cs;
   }

   Internal::Array ALegendreBase::unitd2Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl) 
   {
      Internal::Array cs(3);

      cs(0) = MHD_MP(0.25) * ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm-MHD_MP(1.0))
                           * Internal::Math::sqrt( MHD_MP(1.0) / ( (dl + dm - MHD_MP(1.0))*(dl + dm) ) )
                           * Internal::Math::sqrt( MHD_MP(1.0) / ( (dl - dm + MHD_MP(2.0))*(dl - dm + MHD_MP(1.0)) ) );

      cs(1) = MHD_MP(0.25) * ( ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm+MHD_MP(1.0)) );

      cs(2) = MHD_MP(0.25) * Internal::Math::sqrt( (dl - dm - MHD_MP(1.0))*(dl - dm) )
                           * Internal::Math::sqrt( (dl + dm + MHD_MP(2.0))*(dl + dm + MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::unitsin_1Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = Internal::Math::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))));

      cs(1) = -Internal::Math::sqrt((MHD_MP(2.0)*dl + MHD_MP(1.0))/(MHD_MP(2.0)*dl + MHD_MP(3.0)))*Internal::Math::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/(2.0*dm);

      return cs;
   }

   Internal::Array ALegendreBase::schmidtPmm(const Internal::MHDFloat)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   Internal::Array ALegendreBase::schmidtPm1m(const Internal::MHDFloat dm)
   {
      Internal::Array cs(1);

      cs(0) = Internal::Math::sqrt(MHD_MP(2.0)*dm + MHD_MP(1.0));

      return cs;
   }

   Internal::Array ALegendreBase::schmidtPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = -Internal::Math::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)));
      cs(1) = Internal::Math::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));

      return cs;
   }

   Internal::Array ALegendreBase::schmidtdPl0(const Internal::MHDFloat dl)
   {
      Internal::Array cs(1);
      // needs checking: is the 0.5 factor correct?
      cs(0) = MHD_MP(0.5)*Internal::Math::sqrt(dl*(dl + MHD_MP(1.0)));

      return cs;
   }

   Internal::Array ALegendreBase::schmidtd2P10()
   {
      Internal::Array cs(1);

      cs(0) = -MHD_MP(0.25)* ( ALegendreBase::A1(MHD_MP(1.0), MHD_MP(0.0)) + ALegendreBase::A1(MHD_MP(1.0), MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::schmidtd2Pl0(const Internal::MHDFloat dl)
   {

      throw std::logic_error("schmidt quasi normalized d2Pl0 norm not implemented");
      // This implementation needs proper checking.

      Internal::Array cs(2);

      Internal::MHDFloat dm = MHD_MP(0.0);

      cs(0) = MHD_MP(0.25) * ( ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm + MHD_MP(1.0)) );

      cs(1) = MHD_MP(0.25) * ( MHD_MP(1.0) 
                              + ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm - MHD_MP(1.0)) 
                                 / ( (dl + MHD_MP(2.0)) * (dl + MHD_MP(1.0)) * dl * (dl - MHD_MP(1.0)) ) 
                              ) * Internal::Math::sqrt( (dl + MHD_MP(2.0)) * (dl + MHD_MP(1.0)) * dl * (dl - MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::schmidtd2Pl1(const Internal::MHDFloat dl)
   {
      Internal::MHDFloat dm = MHD_MP(1.0);

      Internal::Array cs(2);

      cs(0) = MHD_MP(0.25) * (ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm - MHD_MP(1.0)) 
                                 / ( (dl + MHD_MP(1.0)) * dl ) 
                              + ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm + MHD_MP(1.0)) );

      cs(1) = MHD_MP(0.25) *  Internal::Math::sqrt( (dl + MHD_MP(3.0)) *(dl + MHD_MP(2.0)) * (dl - MHD_MP(1.0)) * (dl - MHD_MP(2.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::schmidtdPmm(const Internal::MHDFloat dm)
   {
      Internal::Array cs(1);

      cs(0) = -Internal::Math::sqrt(dm/MHD_MP(2.0));

      return cs;
   }

   Internal::Array ALegendreBase::schmidtdPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = -MHD_MP(0.5)*Internal::Math::sqrt((dl - dm + MHD_MP(1.0))*(dl + dm));

      cs(1) = -MHD_MP(0.5)*Internal::Math::sqrt((dl - dm)*(dl + dm + MHD_MP(1.0)));

      return cs;
   }

   Internal::Array ALegendreBase::schmidtd2Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl) 
   {
      Internal::Array cs(3);

      cs(0) = MHD_MP(0.25) * ALegendreBase::A1(dl, dm)*ALegendreBase::A1(dl, dm-MHD_MP(1.0))
                           * Internal::Math::sqrt( MHD_MP(1.0) / ( (dl + dm - MHD_MP(1.0))*(dl + dm) ) )
                           * Internal::Math::sqrt( MHD_MP(1.0) / ( (dl - dm + MHD_MP(2.0))*(dl - dm + MHD_MP(1.0)) ) );

      cs(1) = MHD_MP(0.25) * ( ALegendreBase::A1(dl, dm) + ALegendreBase::A1(dl, dm+MHD_MP(1.0)) );

      cs(2) = MHD_MP(0.25) * Internal::Math::sqrt( (dl - dm - MHD_MP(1.0))*(dl - dm) )
                           * Internal::Math::sqrt( (dl + dm + MHD_MP(2.0))*(dl + dm + MHD_MP(1.0)) );

      return cs;
   }

   Internal::Array ALegendreBase::schmidtsin_1Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl)
   {
      Internal::Array cs(2);

      cs(0) = Internal::Math::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))));

      cs(1) = -Internal::Math::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/MHD_MP(2.0*dm);

      return cs;
   }

}
}
}
