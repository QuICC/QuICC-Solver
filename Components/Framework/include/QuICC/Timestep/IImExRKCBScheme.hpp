/**
 * @file IImExRKCBScheme.hpp
 * @brief Interface for a generic implicit/explicit Runge-Kutta scheme (Cavaglieri & Bewley, 2015)
 */

#ifndef QUICC_TIMESTEP_IIMEXRKCBSCHEME_HPP
#define QUICC_TIMESTEP_IIMEXRKCBSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Timestep/IScheme.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Interface of generic of an implicit/explicit Runge-Kutta scheme from CB
    */
   class IImExRKCBScheme: public IScheme
   {
      public:
         /**
          * @brief Constructor
          */
         IImExRKCBScheme() = default;

         /**
          * @brief Destructor
          */
         virtual ~IImExRKCBScheme() = default;

         /**
          * @brief Butcher's tableau a_ij factor for implicit scheme
          */
         virtual MHDFloat aIm(const int i, const int j) const = 0;

         /**
          * @brief Butcher's tableau b_i factor for implicit scheme
          */
         virtual MHDFloat bIm(const int i) const = 0;

         /**
          * @brief Butcher's tableau a_ij factor for explicit scheme
          */
         virtual MHDFloat aEx(const int i, const int j) const = 0;

         /**
          * @brief Butcher's tableau b_i factor for explicit scheme
          */
         virtual MHDFloat bEx(const int i) const = 0;

         /**
          * @brief Butcher's tableau c_i factor for explicit scheme
          */
         virtual MHDFloat cEx(const int i) const = 0;

         /**
          * @brief Butcher's tableau b_i factor for implicit embedded lower order scheme
          */
         virtual MHDFloat bImErr(const int i) const = 0;

         /**
          * @brief Butcher's tableau b_i factor for explicit embedded lower order scheme
          */
         virtual MHDFloat bExErr(const int i) const = 0;

      protected:

      private:

   };

   /// Typedef for a shared pointer IImExRKCBScheme
   typedef std::shared_ptr<IImExRKCBScheme> SharedIImExRKCBScheme;

}
}

#endif // QUICC_TIMESTEP_IIMEXRKCBSCHEME_HPP
