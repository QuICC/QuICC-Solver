/**
 * @file Dispatchers.cpp
 * @brief Implementation of functors used in Interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/Dispatchers.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include <unsupported/Eigen/SparseExtra>
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX

namespace QuICC {

namespace Equations {

void dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams)
{
   // Get list of implicit fields
   Equations::CouplingInformation::FieldId_range imRange = cinfo.implicitRange();

   auto&& eigs = cinfo.couplingTools().getIndexes(res, matIdx);
   backend.modelMatrix(rModelMatrix, opId, imRange, matIdx, bcType, res, eigs, bcIds, eqParams);

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   auto opName = ModelOperator::Coordinator::tag(opId);
   auto tags =  std::vector<SpectralFieldId>(imRange.first, imRange.second);

   std::vector<int> fileIdx;
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   if(eigs.size() == 1)
   {
      fileIdx.push_back(tRes.idx<Dimensions::Data::DAT3D>(matIdx));
   }
   else
   {
      ArrayI mode = tRes.mode(matIdx);
      fileIdx.push_back(mode(0));
      fileIdx.push_back(mode(1));
   }
   debug::writeModelMatrix(rModelMatrix, opName, tags, fileIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
}

void dispatchGalerkinStencil(const SpectralFieldId fieldId, SparseMatrix &mat, const int matIdx, const Resolution& res, const bool makeSquare, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams)
{
   auto&& eigs = cinfo.couplingTools().getIndexes(res, matIdx);
   backend.galerkinStencil(mat, fieldId, matIdx, res, eigs, makeSquare, bcIds, eqParams);

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   std::string opName = "galerkin_stencil";
   if(makeSquare)
   {
      opName += "_sq";
   }
   std::vector<SpectralFieldId> tags = {fieldId};

   std::vector<int> fileIdx;
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   if(eigs.size() == 1)
   {
      fileIdx.push_back(tRes.idx<Dimensions::Data::DAT3D>(matIdx));
   }
   else
   {
      ArrayI mode = tRes.mode(matIdx);
      fileIdx.push_back(mode(0));
      fileIdx.push_back(mode(1));
   }
   debug::writeModelMatrix(mat, opName, tags, fileIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
}

void dispatchExplicitBlock(const SpectralFieldId fieldId, DecoupledZSparse& mat, const std::size_t opId,  const SpectralFieldId exId, const int matIdx, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation& cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams)
{
   auto&& eigs = cinfo.couplingTools().getIndexes(res, matIdx);
   backend.explicitBlock(mat, fieldId, opId, exId, matIdx, res, eigs, bcIds, eqParams);

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   auto opName = ModelOperator::Coordinator::tag(opId);
   std::vector<SpectralFieldId> tags = {fieldId, exId};

   std::vector<int> fileIdx;
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   if(eigs.size() == 1)
   {
      fileIdx.push_back(tRes.idx<Dimensions::Data::DAT3D>(matIdx));
   }
   else
   {
      ArrayI mode = tRes.mode(matIdx);
      fileIdx.push_back(mode(0));
      fileIdx.push_back(mode(1));
   }
   debug::writeModelMatrix(mat, opName, tags, fileIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
}

void dispatchCoupling(Equations::CouplingInformation& cinfo, const SpectralFieldId fieldId, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features, const Resolution& res, const Model::IModelBackend& backend, const std::map<std::size_t, std::size_t>& bcIds)
{
   bool hasNL = features.at(CouplingFeature::Nonlinear);
   bool hasSource = features.at(CouplingFeature::Source);
   bool hasBoundaryValue = features.at(CouplingFeature::BoundaryValue);
   bool allowExplicit = features.at(CouplingFeature::AllowExplicit);

   Model::EquationInfo eqInfo;
   backend.equationInfo(eqInfo, fieldId, res);

   // Compute effective starting index for local CPU
   int cpuIZero = iZero;
   if(iZero == 1)
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      if(tRes.idx<Dimensions::Data::DAT3D>(0) == 0 && tRes.idx<Dimensions::Data::DAT2D>(0,0) == 0)
      {
         cpuIZero = 1;
      } else
      {
         cpuIZero = 0;
      }
   } else if(iZero > 1)
   {
      throw std::logic_error("Matrix starting index > 1 is not implemented yet!");
   }

   // General setup: equation type? real/complex solver? start from m = ?
   cinfo.setGeneral(eqType, eqInfo.isComplex, cpuIZero, eqInfo.isSplitEquation);

   // Set source flag: has source term?
   cinfo.setSource(hasSource);

   // Set boundary value flag: has boundary value?
   cinfo.setBoundaryValue(hasBoundaryValue);

   // Set index type: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
   auto idxType = safe_CouplingIndexType_cast(eqInfo.indexMode);
   auto spCoupling = res.sim().ss().createCouplingTools(idxType);
   cinfo.setIndexType(idxType, spCoupling);

   // Create implicit field coupling
   int nFields = std::distance(eqInfo.im.begin(), eqInfo.im.end());
   for(auto fIt = eqInfo.im.cbegin(); fIt != eqInfo.im.cend(); ++fIt)
   {
      cinfo.addImplicitField(fIt->first, fIt->second);
   }

   // Create explicit fields
   bool hasQI = false;
   if(allowExplicit)
   {
      // explicit linear
      for(auto fIt = eqInfo.exL.cbegin(); fIt != eqInfo.exL.cend(); ++fIt)
      {
         cinfo.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitLinear::id());
      }

      // explicit nonlinear
      for(auto fIt = eqInfo.exNL.cbegin(); fIt != eqInfo.exNL.cend(); ++fIt)
      {
         if(!(fIt->first == fieldId.first && fIt->second == fieldId.second))
         {
            cinfo.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNonlinear::id());
         }
      }

      // explicit nextstep
      for(auto fIt = eqInfo.exNS.cbegin(); fIt != eqInfo.exNS.cend(); ++fIt)
      {
         cinfo.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNextstep::id());
      }

      // Extract quasi inverse
      auto fIt = std::find(eqInfo.exNL.begin(), eqInfo.exNL.end(), fieldId);
      if(fIt != eqInfo.exNL.end())
      {
         hasQI = true;
      }
   }

   // Set nonlinear flags: has nonlinear term? has quasi-inverse?
   cinfo.setNonlinear(hasNL, hasNL && hasQI);

   // Sort implicit fields
   cinfo.sortImplicitFields(fieldId.first, fieldId.second);

   // Get number of matrices
   int nMat = cinfo.couplingTools().nMat(res);

   // Set field coupling information
   Model::OperatorInfo opInfo(nMat);
   backend.operatorInfo(opInfo, fieldId, res, cinfo.couplingTools(), bcIds);

   cinfo.couplingTools().setTauN(opInfo.tauN, res);
   cinfo.couplingTools().setGalerkinN(opInfo.galN, res);
   cinfo.couplingTools().setRhsN(opInfo.rhsCols, res);
   cinfo.couplingTools().setSystemN(opInfo.sysN, res, nFields);
   cinfo.setSizes(nMat, opInfo.tauN, opInfo.galN, opInfo.galShift, opInfo.rhsCols, opInfo.sysN);
}

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
namespace debug {
   /// Create filename to write model operator to MatrixMarket file
   void filenameWriteModelMatrix(const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids);

   /// Write decoupled complex model operator to MatrixMarket file
   void writeModelMatrix(const DecoupledZSparse& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids);

   /// Write real model operator to MatrixMarket file
   void writeModelMatrix(const SparseMatrix& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids);

   std::string filenameModelMatrix(const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids)
   {
      std::stringstream ss;
      ss << opName;
      for(const auto& f: tags)
      {
         ss << "_" << PhysicalNames::Coordinator::tag(f.first);
         if(f.second == FieldComponents::Spectral::TOR)
         {
            ss <<  "_tor";
         }
         else if(f.second == FieldComponents::Spectral::POL)
         {
            ss <<  "_pol";
         }
      }
      for(auto idx: ids)
      {
         ss << "_" << idx;
      }

      return ss.str();
   }

   void writeModelMatrix(const DecoupledZSparse& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids)
   {
      auto baseName = filenameModelMatrix(opName, tags, ids);

      std::string matName = baseName + "_re.mtx";
      Eigen::saveMarket(mat.real(), matName);
      matName = baseName + "_im.mtx";
      Eigen::saveMarket(mat.imag(), matName);
   }

   void writeModelMatrix(const SparseMatrix& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const std::vector<int>& ids)
   {
      auto baseName = filenameModelMatrix(opName, tags, ids);
      std::string matName = baseName + ".mtx";
      Eigen::saveMarket(mat, matName);
   }
}
#endif

} // namespace Equations   
} // namespace QuICC
