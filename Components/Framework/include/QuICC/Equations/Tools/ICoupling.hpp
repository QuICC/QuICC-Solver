/**
 * @file ICoupling.hpp
 * @brief Interface for the equations coupling tools
 */

#ifndef QUICC_EQUATIONS_TOOLS_ICOUPLING_HPP
#define QUICC_EQUATIONS_TOOLS_ICOUPLING_HPP

// System includes
//
#include <vector>
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   /**
    * @brief Interface for equations coupling tools
    */
   class ICoupling
   {
      public:
         /**
          * @brief Constructor
          */
         ICoupling() = default;

         /**
          * @brief Destructor
          */
         virtual ~ICoupling() = default;

         /**
          * @brief Get index values
          */
         std::vector<MHDFloat> getIndexes(const Resolution& res, const int matIdx) const;

         /**
          * @brief Get the number of matrix operators for field coupling
          *
          * @param res   Shared resolution
          */
         int nMat(const Resolution& res) const;

         /**
          * @brief Set Tau resolution
          */
         void setTauN(ArrayI& rTauNs, const Resolution& res) const;

         /**
          * @brief Set Galerkin resolution provided by python code
          */
         void setGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const;

         /**
          * @brief Set number of RHS provided by python code
          */
         void setRhsN(ArrayI& rRhsCols, const Resolution& res) const;

         /**
          * @brief Set system size provided by python code
          *
          * @param rSystemNs  System sizes
          * @param res        Resolution object
          * @param nFields    Number of coupled fields in system
          */
         void setSystemN(ArrayI& rSystemNs, const Resolution& res, const int nFields) const;

      private:
         /**
          * @brief Set of index values
          */
         virtual std::vector<MHDFloat> identifyIndexes(const Resolution& res, const int matIdx) const = 0;

         /**
          * @brief Compute the number of matrix operators for field coupling
          *
          * @param res   Shared resolution
          */
         virtual int computeNMat(const Resolution& res) const = 0;

         /**
          * @brief Interpret Tau resolution provided by python code
          */
         virtual void interpretTauN(ArrayI& rTauNs, const Resolution& res) const = 0;

         /**
          * @brief Interpret Galerkin resolution provided by python code
          */
         virtual void interpretGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const = 0;

         /**
          * @brief Interpret number of RHS provided by python code
          */
         virtual void interpretRhsN(ArrayI& rRhsCols, const Resolution& res) const = 0;

         /**
          * @brief Interpret system size provided by python code
          */
         virtual void interpretSystemN(ArrayI& rSystemNs, const Resolution& res, const int nFields) const = 0;
   };

   /// Typedef for a shared ICoupling
   typedef std::shared_ptr<ICoupling> SharedICoupling;

}
}
}

#endif // QUICC_EQUATIONS_ICOUPLING_HPP
