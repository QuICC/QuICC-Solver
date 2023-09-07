/**
 * @file Regular1D.hpp
 * @brief Implementation of some tools for schemes with a regular single dimension index
 */

#ifndef QUICC_EQUATIONS_TOOLS_REGULAR1D_HPP
#define QUICC_EQUATIONS_TOOLS_REGULAR1D_HPP

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
    * @brief Tools for equations with a regular single dimensional index
    */
   class Regular1D: public ICoupling
   {
      public:
         /**
          * @brief Constructor
          */
         Regular1D() = default;

         /**
          * @brief Destructor
          */
         virtual ~Regular1D() = default;

      private:
         /**
          * @brief Set index values
          */
         virtual std::vector<MHDFloat> identifyIndexes(const Resolution& res, const int matIdx) const;

         /**
          * @brief Compute the number of matrix operators for field coupling
          *
          * @param res   Shared resolution
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

   /// Typedef for a shared Regular1D
   typedef std::shared_ptr<Regular1D> SharedRegular1D;

}
}
}

#endif // QUICC_EQUATIONS_TOOLS_REGULAR1D_HPP
