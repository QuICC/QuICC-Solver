/**
 * @file Regular3D.hpp
 * @brief Implementation of some tools for schemes with regular three dimensional index
 */

#ifndef QUICC_EQUATIONS_TOOLS_REGULAR3D_HPP
#define QUICC_EQUATIONS_TOOLS_REGULAR3D_HPP

// System includes
//
#include<vector>
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/Tools/ICoupling.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

/**
 * @brief Tools for equations with regular three dimensional index
 */
   class Regular3D: public ICoupling
   {
      public:
         /**
          * @brief Constructor
          */
         Regular3D() = default;

         /**
          * @brief Destructor
          */
         virtual ~Regular3D() = default;

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

   /// Typedef for a shared Regular3D
   typedef std::shared_ptr<Regular3D> SharedRegular3D;

}
}
}

#endif // QUICC_EQUATIONS_TOOLS_REGULAR3D_HPP
