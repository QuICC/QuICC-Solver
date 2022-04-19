/** 
 * @file ALegendreBase.cpp
 * @brief Source of the implementation of the associated Legendre polynomial
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Polynomial/ALegendre/ALegendreBase.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   ALegendreBase::ALegendreBase()
   {
   }

   ALegendreBase::~ALegendreBase()
   {
   }

   void ALegendreBase::Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid, NormalizerML norm)
   {
      // Safety assert
      assert(l-m > 0);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);
      internal::Array cs = norm(dm, dl);

      iplm.array() = cs(1)*(cs(0)*ipl_2m.array() + (igrid.array()*ipl_1m.array()).array());
   }
   
   void ALegendreBase::Pmm(Eigen::Ref<internal::Matrix> ipmm, const int m, const internal::Array& igrid, NormalizerM norm)
   {
      internal::MHDFloat dm = internal::MHDFloat(m);
      internal::Array cs = norm(dm);

      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         ipmm.setConstant(cs(0));
      } else
      {
         internal::MHDFloat di = MHD_MP(1.0);

         for(int i = 1; i <= m; i++)
         {
            cs(0) *= -precision::sqrt(di/(di + MHD_MP(1.0)));
            di += MHD_MP(2.0);
         }
         ipmm.array() = cs(0)*(MHD_MP(1.0) - igrid.array().pow(2)).pow(dm/2);
      }
   }

   void ALegendreBase::Pm1m(Eigen::Ref<internal::Matrix> ipm1m, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid, NormalizerM norm)
   {
      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         internal::MHDFloat dm = internal::MHDFloat(m);
         internal::Array cs = norm(dm);

         ipm1m.array() = (cs(0)*igrid).array()*ipmm.array();
      }
   }

   void ALegendreBase::dPl0(Eigen::Ref<internal::Matrix> idpl0, const int l, const Eigen::Ref<const internal::Matrix>& ipl1, NormalizerL norm)
   {
      // Safety assert
      assert(l > 0);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::Array cs = norm(dl);

      idpl0.array() = cs(0)*ipl1.array();
   }

   void ALegendreBase::dPmm(Eigen::Ref<internal::Array> idpmm, const int m, const Eigen::Ref<const internal::Array>& iplm_1, NormalizerM norm)
   {
      if(m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         idpmm.setConstant(MHD_MP(0.0));

      } else
      {
         internal::MHDFloat dm = internal::MHDFloat(m);
         internal::Array cs = norm(dm);

         idpmm = cs(0)*iplm_1;
      }
   }

   void ALegendreBase::dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& iplm_1, const Eigen::Ref<const internal::Matrix>& iplm1, NormalizerML norm)
   {
      // Safety assert
      assert(l-m > 0);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);
      internal::Array cs = norm(dm, dl);
 
      idplm.array() = cs(0)*iplm_1.array() - cs(1)*iplm1.array();
   }

   void ALegendreBase::sin_1Plm(Eigen::Ref<internal::Matrix> isin_1plm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl1m1, const Eigen::Ref<const internal::Matrix>& ipl1m_1, NormalizerML norm)
   {
      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);
      internal::Array cs = norm(dm, dl);

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

   ALegendreBase::NormalizerML ALegendreBase::normdPlm()
   {
      #if defined QUICC_ALEGENDRE_NORM_SHSCHMIDT
         return &ALegendreBase::schmidtdPlm;
      #elif defined QUICC_ALEGENDRE_NORM_SHUNITY
         return &ALegendreBase::unitdPlm;
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

   internal::Array ALegendreBase::unitPmm(const internal::MHDFloat dm)
   {
      internal::Array cs(1);

      if(dm == MHD_MP(0.0))
      {
         cs(0) = precision::sqrt(MHD_MP(1.0)/(MHD_MP(4.0)*Precision::PI));
      } else
      {
         cs(0) = precision::sqrt((MHD_MP(2.0)*dm + MHD_MP(1.0))/(MHD_MP(4.0)*Precision::PI));
      }

      return cs;
   }

   internal::Array ALegendreBase::unitPm1m(const internal::MHDFloat dm)
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(2.0)*dm + MHD_MP(3.0));

      return cs;
   }

   internal::Array ALegendreBase::unitPlm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = -precision::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)));
      cs(1) = precision::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));

      return cs;
   }

   internal::Array ALegendreBase::unitdPl0(const internal::MHDFloat dl)
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(dl*(dl + MHD_MP(1.0)));

      return cs;
   }

   internal::Array ALegendreBase::unitdPmm(const internal::MHDFloat dm)
   {
      internal::Array cs(1);

      cs(0) = -precision::sqrt(dm/MHD_MP(2.0));

      return cs;
   }

   internal::Array ALegendreBase::unitdPlm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = -MHD_MP(0.5)*precision::sqrt((dl - dm + MHD_MP(1.0))*(dl + dm));

      cs(1) = -MHD_MP(0.5)*precision::sqrt((dl - dm)*(dl + dm + MHD_MP(1.0)));

      return cs;
   }

   internal::Array ALegendreBase::unitsin_1Plm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = precision::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))));

      cs(1) = -precision::sqrt((MHD_MP(2.0)*dl + MHD_MP(1.0))/(MHD_MP(2.0)*dl + MHD_MP(3.0)))*precision::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/(2.0*dm);

      return cs;
   }

   internal::Array ALegendreBase::schmidtPmm(const internal::MHDFloat)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   internal::Array ALegendreBase::schmidtPm1m(const internal::MHDFloat dm)
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(2.0)*dm + MHD_MP(1.0));

      return cs;
   }

   internal::Array ALegendreBase::schmidtPlm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = -precision::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)));
      cs(1) = precision::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));

      return cs;
   }

   internal::Array ALegendreBase::schmidtdPl0(const internal::MHDFloat dl)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.5)*precision::sqrt(dl*(dl + MHD_MP(1.0)));

      return cs;
   }

   internal::Array ALegendreBase::schmidtdPmm(const internal::MHDFloat dm)
   {
      internal::Array cs(1);

      cs(0) = -precision::sqrt(dm/MHD_MP(2.0));

      return cs;
   }

   internal::Array ALegendreBase::schmidtdPlm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = -MHD_MP(0.5)*precision::sqrt((dl - dm + MHD_MP(1.0))*(dl + dm));

      cs(1) = -MHD_MP(0.5)*precision::sqrt((dl - dm)*(dl + dm + MHD_MP(1.0)));

      return cs;
   }

   internal::Array ALegendreBase::schmidtsin_1Plm(const internal::MHDFloat dm, const internal::MHDFloat dl)
   {
      internal::Array cs(2);

      cs(0) = precision::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))));

      cs(1) = -precision::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/MHD_MP(2.0*dm);

      return cs;
   }

}
}
}
