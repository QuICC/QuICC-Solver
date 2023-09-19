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
#include "Types/Precision.hpp"

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
         typedef std::function<internal::Array(const internal::MHDFloat)> NormalizerM;

         /// Typedef for a l dependent normalizer
         typedef std::function<internal::Array(const internal::MHDFloat)> NormalizerL;

         /// Typedef for a m and l dependent normalizer
         typedef std::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat)> NormalizerML;
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
         static void Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid, NormalizerML norm);

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
         static void Pmm(Eigen::Ref<internal::Matrix> ipmm, const int m, const internal::Array& igrid, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$P_{m+1}^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pm1m(Eigen::Ref<internal::Matrix> ipm1m, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_l^m (\cos\theta)\f$ (scheme B)
          *
          * Internal computation can be done in multiple precision. Uses recurrence relation for the wanted expression.
          */
         static void dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& iplm_1, const Eigen::Ref<const internal::Matrix>& iplm1, NormalizerML norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_l (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision. Uses recurrence relation for the wanted expression.
          */
         static void dPl0(Eigen::Ref<internal::Matrix> idplm, const int l, const Eigen::Ref<const internal::Matrix>& iplm1, NormalizerL norm);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_m^m (\cos\theta)\f$ (scheme B)
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& iplm_1, NormalizerM norm);

         /**
          * @brief Compute the associated Legendre \f$1/\sin\theta P_m^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void sin_1Plm(Eigen::Ref<internal::Matrix> isin_1plm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl1m1, const Eigen::Ref<const internal::Matrix>& ipl1m_1, NormalizerML norm);

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
         static internal::Array unitPmm(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitPm1m(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitPlm(const internal::MHDFloat dm, const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitdPl0(const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitdPmm(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitdPlm(const internal::MHDFloat dm, const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for unit Associated Legendre normalization
          */
         static internal::Array unitsin_1Plm(const internal::MHDFloat dm, const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtPmm(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtPm1m(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtPlm(const internal::MHDFloat dm, const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtdPl0(const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtdPmm(const internal::MHDFloat dm);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtdPlm(const internal::MHDFloat dm, const internal::MHDFloat dl);

         /**
          * @brief Polynomial normalizer for Schmidt Associated Legendre quasi-normalization
          */
         static internal::Array schmidtsin_1Plm(const internal::MHDFloat dm, const internal::MHDFloat dl);

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
