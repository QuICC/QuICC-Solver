/**
 * @file LegendreBase.hpp
 * @brief Implementation of the Legendre polynomial base
 */

#ifndef QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP
#define QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP

// System includes
//
#include <functional>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

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
         typedef std::function<Internal::Array()> Normalizer;

         /// Typedef for a n dependent normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat)> NormalizerN;

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
         static void P0(Eigen::Ref<Internal::Matrix> ip0, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$P_1 (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::Array& igrid, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$P_n (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerN norm);

         /**
          * @brief Compute the Legendre \f$\frac{d}{dx} P_n (x)\f$
          *
          * Internal computation can be done in multiple precision.
          */
         static void dP1(Eigen::Ref<Internal::Matrix> idp1, const Eigen::Ref<const Internal::Matrix>& ipn_1, Normalizer norm);

         /**
          * @brief Compute the Legendre \f$\frac{d}{dx} P_n (x)\f$
          *
          * Internal computation can be done in multiple precision.
          */
         static void dPn(Eigen::Ref<Internal::Matrix> idpn, const int n, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& idpn_1, const Internal::Array& igrid, NormalizerN norm);

         static Normalizer normP0();
         static Normalizer normP1();
         static NormalizerN normPn();
         static Normalizer normdP1();
         static NormalizerN normdPn();

      protected:
         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static Internal::Array unitP0();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static Internal::Array unitP1();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static Internal::Array unitPn(const Internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static Internal::Array unitdP1();

         /**
          * @brief Polynomial normalizer for unit Legendre normalization
          */
         static Internal::Array unitdPn(const Internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static Internal::Array naturalP0();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static Internal::Array naturalP1();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static Internal::Array naturalPn(const Internal::MHDFloat dn);

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static Internal::Array naturaldP1();

         /**
          * @brief Polynomial normalizer for natural Legendre normalization
          */
         static Internal::Array naturaldPn(const Internal::MHDFloat dn);

      private:

   };
}
}
}

#endif // QUICC_POLYNOMIAL_LEGENDRE_LEGENDREBASE_HPP
