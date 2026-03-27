/**
 * @file Dispatchers.hpp
 * @brief Implementation of functors used in Interface
 */

#pragma once

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

void  dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams);

void dispatchGalerkinStencil(const std::size_t fieldName, FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const Resolution& res, const bool makeSquare, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams);

void dispatchExplicitBlock(const std::size_t fieldName, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams);

void dispatchCoupling(Equations::CouplingInformation& cinfo, std::size_t fieldName, FieldComponents::Spectral::Id compId, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features, const Resolution& res, const Model::IModelBackend& backend, const std::map<std::size_t, std::size_t>& bcIds);

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
namespace debug {
   /// Create filename to write model operator to MatrixMarket file
   void filenameWriteModelMatrix(const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids);

   /// Write decoupled complex model operator to MatrixMarket file
   void writeModelMatrix(const DecoupledZSparse& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids);
}
#endif

} // namespace Equations
} // namespace QuICC
