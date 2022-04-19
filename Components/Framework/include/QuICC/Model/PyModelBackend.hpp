/** 
 * @file PyModelBackend.hpp
 * @brief Interface for a Python model backend
 */

#ifndef QUICC_MODEL_PYMODELBACKEND_HPP
#define QUICC_MODEL_PYMODELBACKEND_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/PyQuICC/ModelWrapper.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for Python model backend
    */
   class PyModelBackend: public IModelBackend
   {
      public:
         /**
          * @brief Constructor
          */
         PyModelBackend(const std::string pyModule, const std::string pyClass);

         /**
          * @brief Destructor
          */
         virtual ~PyModelBackend() = default;

         /**
          * @brief Get vector of names for the physical fields
          */
         virtual std::vector<std::string> fieldNames() const override;

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         virtual std::vector<std::string> paramNames() const override;

         /**
          * @brief Get vector of bools about periodic box
          */
         virtual std::vector<bool> isPeriodicBox() const override;

         /**
          * @brief Enable galerkin basis
          */
         virtual void enableGalerkin(const bool flag) override;

         /**
          * @brief Get vector of bools about periodic box
          */
         virtual std::map<std::string,MHDFloat> automaticParameters(const std::map<std::string,MHDFloat>& cfg) const override;

         /**
          * @brief Get equation information
          */
         virtual void equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const override;

         /**
          * @brief Get operator information
          */
         virtual void operatorInfo(ArrayI& tauN, ArrayI& galN, MatrixI& galShift, ArrayI& rhsCols, ArrayI& sysN, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const override;

         /**
          * @brief Build model matrix
          */
         virtual void modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build galerkin stencil
          */
         virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build explicit block
          */
         virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

      protected:
         /**
          * @brief Build base arguments for Python wrapper
          */
         PyObject*  baseArguments(const int tupleSize, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, BcMap bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Shared pointer to Python model wrapper
          */
         std::shared_ptr<PyQuICC::ModelWrapper>   mpWrapper;

      private:
   };

}
}

#endif // QUICC_MODEL_PYMODELBACKEND_HPP
