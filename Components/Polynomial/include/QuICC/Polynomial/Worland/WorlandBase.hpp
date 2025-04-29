/**
 * @file WorlandBase.hpp
 * @brief Implementation of the Worland polynomial base
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP
#define QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the Worland polynomial base
    */
   class WorlandBase
   {
      public:
         /**
          * @brief Constructor
          */
         WorlandBase();

         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         WorlandBase(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta);

         /**
          * @brief Destructor
          */
         virtual ~WorlandBase() = default;

         /**
          * @brief Get alpha parameter of Jacobi polynomial
          *
          * @param l Harmonic degree l
          */
         Internal::MHDFloat alpha(const int l);

         /**
          * @brief Get dBeta (beta = f(l) + dBeta) parameter of Jacobi polynomial
          */
         Internal::MHDFloat dBeta();

      protected:
         /**
          * @brief Get beta parameter of Jacobi polynomial
          *
          * @param l Harmonic degree l
          */
         Internal::MHDFloat beta(const int l);

         /**
          * @brief Compute \f$W_0^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         void computeW0l(Eigen::Ref<Internal::Matrix> iw0ab, const int l, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, ThreeTermRecurrence::NormalizerAB norm);

         static ThreeTermRecurrence::NormalizerNAB normWPnab();
         static ThreeTermRecurrence::NormalizerAB normWP1ab();
         static ThreeTermRecurrence::NormalizerAB normWP0ab();
         static ThreeTermRecurrence::NormalizerNAB normWDPnab();
         static ThreeTermRecurrence::NormalizerAB normWDP1ab();
         static ThreeTermRecurrence::NormalizerAB normWDP0ab();
         static ThreeTermRecurrence::NormalizerNAB normWD2Pnab();
         static ThreeTermRecurrence::NormalizerAB normWD2P1ab();
         static ThreeTermRecurrence::NormalizerAB normWD2P0ab();
         static ThreeTermRecurrence::NormalizerNAB normWD3Pnab();
         static ThreeTermRecurrence::NormalizerAB normWD3P1ab();
         static ThreeTermRecurrence::NormalizerAB normWD3P0ab();

      private:
         /**
          * @brief Polynomial normalizer for unit Worland normalization
          */
         static Internal::Array unitWPnab(const Internal::MHDFloat dn, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Polynomial n=0 normalizer for unit Worland normalization
          */
         static Internal::Array unitWP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Polynomial n=1 normalizer for unit Worland normalization
          */
         static Internal::Array unitWP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative normalizer for unit Worland normalization
          */
         static Internal::Array unitWDPnab(const Internal::MHDFloat dn, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative n=0 normalizer for unit Worland normalization
          */
         static Internal::Array unitWDP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative n=1 normalizer for unit Worland normalization
          */
         static Internal::Array unitWDP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative normalizer for unit Worland normalization
          */
         static Internal::Array unitWD2Pnab(const Internal::MHDFloat dn, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative n=0 normalizer for unit Worland normalization
          */
         static Internal::Array unitWD2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative n=1 normalizer for unit Worland normalization
          */
         static Internal::Array unitWD2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative normalizer for unit Worland normalization
          */
         static Internal::Array unitWD3Pnab(const Internal::MHDFloat dn, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative n=0 normalizer for unit Worland normalization
          */
         static Internal::Array unitWD3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative n=1 normalizer for unit Worland normalization
          */
         static Internal::Array unitWD3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Polynomial normalizer for natural normalization
          */
         static Internal::Array naturalWPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Polynomial n=0 normalizer for natural normalization
          */
         static Internal::Array naturalWP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Polynomial n=1 normalizer for natural normalization
          */
         static Internal::Array naturalWP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative normalizer for natural normalization
          */
         static Internal::Array naturalWDPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative n=0 normalizer for natural normalization
          */
         static Internal::Array naturalWDP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief First derivative n=1 normalizer for natural normalization
          */
         static Internal::Array naturalWDP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative normalizer for natural normalization
          */
         static Internal::Array naturalWD2Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative n=0 normalizer for natural normalization
          */
         static Internal::Array naturalWD2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Second derivative n=1 normalizer for natural normalization
          */
         static Internal::Array naturalWD2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative normalizer for natural normalization
          */
         static Internal::Array naturalWD3Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative n=0 normalizer for natural normalization
          */
         static Internal::Array naturalWD3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Third derivative n=1 normalizer for natural normalization
          */
         static Internal::Array naturalWD3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b);

         /**
          * @brief Alpha parameter
          */
         Internal::MHDFloat mAlpha;

         /**
          * @brief Constant part of beta parameter
          */
         Internal::MHDFloat mDBeta;

   };
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP
