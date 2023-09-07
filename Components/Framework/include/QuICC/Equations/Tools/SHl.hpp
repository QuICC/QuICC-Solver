/**
 * @file SHl.hpp
 * @brief Implementation of some tools for schemes with spherical harmonics expansions with l spectral ordering
 */

#ifndef QUICC_EQUATIONS_TOOLS_SHL_HPP
#define QUICC_EQUATIONS_TOOLS_SHL_HPP

// System includes
//
#include<vector>
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/Tools/ICoupling.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   /**
    * @brief Tools for equations with spherical harmonic expansions with l spectral ordering
    */
   class SHl: public ICoupling
   {
      public:
         /**
          * @brief Constructor
          */
         SHl() = default;

         /**
          * @brief Destructor
          */
         virtual ~SHl() = default;

      private:
         /**
          * @brief Set index values
          */
         virtual std::vector<MHDFloat> identifyIndexes(const Resolution& res, const int matIdx) const;

         /**
          * @brief Compute the number of matrix operators for field coupling
          *
          * @param spRes   Shared resolution
          */
         virtual int computeNMat(const Resolution& res) const;

         /**
          * @brief Interpret Tau resolution provided by python code
          */
         virtual void interpretTauN(ArrayI& rTauNs, const Resolution& res) const;

         /**
          * @brief Interpret Galerkin resolution provided by python code
          */
         virtual void interpretGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const;

         /**
          * @brief Interpret number of RHS provided by python code
          */
         virtual void interpretRhsN(ArrayI& rRhsCols, const Resolution& res) const;

         /**
          * @brief Interpret system size provided by python code
          *
          * @param rSystemNs  System sizes
          * @param res        Resolution object
          * @param nFields    Number of coupled fields in system
          */
         virtual void interpretSystemN(ArrayI& rSystemNs, const Resolution& res, const int nFields) const;

   };

   /// Typedef for a shared SHl
   typedef std::shared_ptr<SHl> SharedSHl;

}
}
}

#endif // QUICC_EQUATIONS_TOOLS_SHL_HPP
