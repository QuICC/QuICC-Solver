/**
 * @file ALegendreBase.hpp
 * @brief Implementation of the associated Legendre polynomial base
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_ALEGENDREBASE_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_ALEGENDREBASE_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <functional>

// External includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   /**
    * @brief Implementation of the associated Legendre polynomial base
    */
   class ALegendreBase
   {
      public:
         /// Typedef for a m dependent normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat)> NormalizerM;

         /// Typedef for a l dependent normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat)> NormalizerL;

         /// Typedef for a m and l dependent normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat, const Internal::MHDFloat)> NormalizerML;
         /**
          * @brief Constructor
          */
         ALegendreBase();

         /**
          * @brief Destructor
          */
         virtual ~ALegendreBase();

         /**
          * @brief Compute the associated Legendre \f$P_l^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Plm(Eigen::Ref<Internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& ipl_1m, const Eigen::Ref<const Internal::Matrix>& ipl_2m, const Internal::Array& igrid, NormalizerML norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{P_l^m (\cos\theta)}{\sin\theta}\f$
          *
          * Internal computation can be done in multiple precision
          */
         template <typename T> void computeSin_1Plm(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > ipl1m, const int m, const int l, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& ipl1m1, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& ipl1m_1);

         /**
          * @brief Compute the associated Legendre \f$P_m^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pmm(Eigen::Ref<Internal::Matrix> ipmm, const int m, const Internal::Array& igrid, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$P_{m+1}^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pm1m(Eigen::Ref<Internal::Matrix> ipm1m, const int m, const Eigen::Ref<const Internal::Matrix>& ipmm, const Internal::Array& igrid, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_l^m (\cos\theta)\f$ (scheme B)
          *
          * Internal computation can be done in multiple precision. Uses recurrence relation for the wanted expression.
          */
         static void dPlm(Eigen::Ref<Internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& iplm_1, const Eigen::Ref<const Internal::Matrix>& iplm1, NormalizerML norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_l (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision. Uses recurrence relation for the wanted expression.
          */
         static void dPl0(Eigen::Ref<Internal::Matrix> idplm, const int l, const Eigen::Ref<const Internal::Matrix>& iplm1, NormalizerL norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_m^m (\cos\theta)\f$ (scheme B)
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm(Eigen::Ref<Internal::Array> op, const int m, const Eigen::Ref<const Internal::Array>& iplm_1, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$1/\sin\theta P_m^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void sin_1Plm(Eigen::Ref<Internal::Matrix> isin_1plm, const int m, const int l, const Eigen::Ref<const Internal::Matrix>& ipl1m1, const Eigen::Ref<const Internal::Matrix>& ipl1m_1, NormalizerML norm);

         static NormalizerM normPmm();
         static NormalizerM normPm1m();
         static NormalizerML normPlm();
         static NormalizerM normdPmm();
         static NormalizerL normdPl0();
         static NormalizerML normdPlm();
         static NormalizerML normsin_1Plm();

      protected:
         /**
          * @brief rescale column with given array
          */
         template <typename T> void rescale(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > iplm, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& scale);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitPmm(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitPm1m(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitdPl0(const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitdPmm(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitdPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static Internal::Array unitsin_1Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtPmm(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtPm1m(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtdPl0(const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtdPmm(const Internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtdPlm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static Internal::Array schmidtsin_1Plm(const Internal::MHDFloat dm, const Internal::MHDFloat dl);

      private:

   };

   template <typename T> void ALegendreBase::rescale(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > iplm, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& scale)
   {
      if(scale.size() > 0)
      {
         if(scale.cols() == 1)
         {
            iplm = scale.topRows(iplm.rows()).asDiagonal()*iplm;
         } else
         {
            assert(iplm.rows() == scale.rows() && iplm.cols() == scale.cols());
            iplm.array() = scale.array()*iplm.array();
         }
      }
   }
}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_ALEGENDREBASE_HPP
