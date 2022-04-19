/** 
 * @file WorlandBase.hpp
 * @brief Implementation of the Worland polynomial base
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP
#define QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
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
          */
         WorlandBase(const internal::MHDFloat alpha, const internal::MHDFloat dBeta);

         /**
          * @brief Destructor
          */
         virtual ~WorlandBase();

         /**
          * @brief Get alpha parameter of Jacobi polynomial
          */
         internal::MHDFloat alpha(const int l);

         /**
          * @brief Get dBeta (beta = f(l) + dBeta) parameter of Jacobi polynomial
          */
         internal::MHDFloat dBeta();

         /// Alpha parameter for Chebyshev type
         static const internal::MHDFloat ALPHA_CHEBYSHEV;
         /// dBeta parameter for Chebyshev type
         static const internal::MHDFloat DBETA_CHEBYSHEV;

         /// Alpha parameter for Legendre type
         static const internal::MHDFloat ALPHA_LEGENDRE;
         /// dBeta parameter for Legendre type
         static const internal::MHDFloat DBETA_LEGENDRE;

         /// Alpha parameter for CylEnergy type
         static const internal::MHDFloat ALPHA_CYLENERGY;
         /// dBeta parameter for CylEnergy type
         static const internal::MHDFloat DBETA_CYLENERGY;

         /// Alpha parameter for SphEnergy type
         static const internal::MHDFloat ALPHA_SPHENERGY;
         /// dBeta parameter for SphEnergy type
         static const internal::MHDFloat DBETA_SPHENERGY;

      protected:
         /**
          * @brief Get beta parameter of Jacobi polynomial
          */
         internal::MHDFloat beta(const int l);

         /**
          * @brief Compute \f$W_0^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         void computeW0l(Eigen::Ref<internal::Matrix> iw0ab, const int l, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, ThreeTermRecurrence::NormalizerAB norm);

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
         static internal::Array unitWPnab(const internal::MHDFloat dn, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWDPnab(const internal::MHDFloat dn, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWDP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWDP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWD2Pnab(const internal::MHDFloat dn, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWD3Pnab(const internal::MHDFloat dn, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial normalizer for natural normalization
          */
         static internal::Array naturalWPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=0 normalizer for natural normalization
          */
         static internal::Array naturalWP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=1 normalizer for natural normalization
          */
         static internal::Array naturalWP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative normalizer for natural normalization
          */
         static internal::Array naturalWDPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalWDP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalWDP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative normalizer for natural normalization
          */
         static internal::Array naturalWD2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalWD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalWD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative normalizer for natural normalization
          */
         static internal::Array naturalWD3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalWD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalWD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Alpha parameter
          */
         internal::MHDFloat mAlpha;

         /**
          * @brief Constant part of beta parameter
          */
         internal::MHDFloat mDBeta;

   };
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_WORLANDBASE_HPP
