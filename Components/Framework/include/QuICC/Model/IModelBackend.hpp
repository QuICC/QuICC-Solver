/** 
 * @file IModelBackend.hpp
 * @brief Interface for a model backend
 */

#ifndef QUICC_MODEL_IMODELBACKEND_HPP
#define QUICC_MODEL_IMODELBACKEND_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for a model backend
    */
   class IModelBackend
   {
      public:
         /// List of field ids
         typedef std::vector<SpectralFieldId> SpectralFieldIds;
         /// Map for boundary conditions
         typedef std::map<std::string,int> BcMap;

         /**
          * @brief Constructor
          */
         IModelBackend() = default;

         /**
          * @brief Destructor
          */
         virtual ~IModelBackend() = default;

         /**
          * @brief Get vector of IDs for the physical fields
          */
         std::vector<std::size_t> fieldIds() const;

         /**
          * @brief Get vector of IDs for the nondimensional parameters
          */
         std::vector<std::size_t> paramIds() const;

         /**
          * @brief Get vector of names for the physical fields
          */
         virtual std::vector<std::string> fieldNames() const = 0;

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         virtual std::vector<std::string> paramNames() const = 0;

         /**
          * @brief Get vector of bools about periodic box
          */
         virtual std::vector<bool> isPeriodicBox() const = 0;

         /**
          * @brief Enable galerkin basis
          */
         virtual void enableGalerkin(const bool flag) = 0;

         /**
          * @brief Get automatically computed parameters
          */
         virtual std::map<std::string,MHDFloat> automaticParameters(const std::map<std::string,MHDFloat>& cfg) const;

         /**
          * @brief Get equation information
          */
         virtual void equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const = 0;

         /**
          * @brief Get operator information
          */
         virtual void operatorInfo(ArrayI& tauN, ArrayI& galN, MatrixI& galShift, ArrayI& rhsCols, ArrayI& sysN, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const = 0;

         /**
          * @brief Build model matrix
          */
         virtual void modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const = 0;

         /**
          * @brief Build galerkin stencil
          */
         virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const = 0;

         /**
          * @brief Build explicit block
          */
         virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const = 0;

      protected:
         /**
          * @brief Validate field names
          */
         void checkFieldNames(const std::vector<std::string>& names) const;

         /**
          * @brief Validate parameters
          */
         void checkParamNames(const std::vector<std::string>& names) const;

      private:
   };

}
}

#endif // QUICC_MODEL_IMODELBACKEND_HPP
