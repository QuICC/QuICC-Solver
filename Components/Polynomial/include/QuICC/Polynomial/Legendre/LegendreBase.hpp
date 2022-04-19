/** 
 * @file LegendreBase.hpp
 * @brief Implementation of the Legendre polynomial base
 */

#ifndef QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP
#define QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP

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
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Polynomial {

namespace Legendre {

   /**
    * @brief Implementation of the Legendre polynomial base
    */ 
   class LegendreBase
   {
      public:
         /// Typedef for a constant normalizer
         typedef std::function<internal::Array()> Normalizer;

         /// Typedef for a n dependent normalizer
         typedef std::function<internal::Array(const internal::MHDFloat)> NormalizerN;

         /**
          * @brief Constructor
          */
         LegendreBase();

         /**
          * @brief Destructor
          */
         virtual ~LegendreBase();

         /**
          * @brief Compute the Legendre \f$P_0 (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<internal::Matrix> ip0, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$P_1 (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<internal::Matrix> ip1, const internal::Array& igrid, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$P_n (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerN norm);

         /**
          * @brief Compute the Legendre \f$\frac{d}{dx} P_n (x)\f$
          *
          * Internal computation can be done in multiple precision.
          */
         static void dP1(Eigen::Ref<internal::Matrix> idp1, const Eigen::Ref<const internal::Matrix>& ipn_1, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$\frac{d}{dx} P_n (x)\f$
          *
          * Internal computation can be done in multiple precision.
          */
         static void dPn(Eigen::Ref<internal::Matrix> idpn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& idpn_1, const internal::Array& igrid, NormalizerN norm);

         static Normalizer normP0();
         static Normalizer normP1();
         static NormalizerN normPn();
         static Normalizer normdP1();
         static NormalizerN normdPn();

      protected:
         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static internal::Array unitP0();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static internal::Array unitP1();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static internal::Array unitPn(const internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static internal::Array unitdP1();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static internal::Array unitdPn(const internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static internal::Array naturalP0();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static internal::Array naturalP1();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static internal::Array naturalPn(const internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static internal::Array naturaldP1();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static internal::Array naturaldPn(const internal::MHDFloat dn);

      private:

   };
}
}
}

#endif // QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP
