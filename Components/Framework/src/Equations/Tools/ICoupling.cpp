/**
 * @file ICoupling.cpp
 * @brief Source of the interface to the eigen direction tools
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/ICoupling.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> ICoupling::getIndexes(const Resolution& res, const int matIdx) const
   {
      return this->identifyIndexes(res, matIdx);
   }

   int ICoupling::nMat(const Resolution& res) const
   {
      return this->computeNMat(res);
   }

   void ICoupling::setTauN(ArrayI& rTauNs, const Resolution& res) const
   {
      this->interpretTauN(rTauNs, res);
   }

   void ICoupling::setGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const
   {
      this->interpretGalerkinN(rGalerkinNs, res);
   }

   void ICoupling::setRhsN(ArrayI& rRhsCols, const Resolution& res) const
   {
      this->interpretRhsN(rRhsCols, res);
   }

   void ICoupling::setSystemN(ArrayI& rSystemNs, const Resolution& res) const
   {
      this->interpretSystemN(rSystemNs, res);
   }

}
}
}
