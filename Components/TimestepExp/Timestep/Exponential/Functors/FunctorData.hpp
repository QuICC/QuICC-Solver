/**
 * @file FunctorData.hpp
 * @brief Data for functors
 */

#pragma once

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "QuICC/Equations/SolutionUpdater.hpp"
#include "Timestep/Exponential/MinEqInfo.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Data for functors
 */
struct FunctorData
{
   /**
    * @brief ctor
    */
   FunctorData() = default;

   /**
    * @brief dtor
    */
   ~FunctorData() = default;

   /**
    * @brief Global resolution
    */
   std::shared_ptr<Resolution> spRes;

   /**
    * @brief Resolution
    */
   const Resolution& res() const {return *spRes;};

   /**
    * @brief Map of boundary condition IDs
    */
   std::map<std::size_t, std::size_t> bcIdMap;

   /**
    * @brief Map of equation parameters
    */
   std::map<std::size_t, NonDimensional::SharedINumber> eqParamsMap;

   /**
    * @brief Model backend
    */
   std::shared_ptr<Model::IModelBackend> spBackend;

   /**
    * @brief Model backend
    */
   const Model::IModelBackend& backend() const {return *spBackend;};

   /**
    * @brief Map of coupling information
    */
   std::map<SpectralFieldId, Equations::CouplingInformation> cInfos; 

   /**
    * @brief Map of scalar field pointers
    */
   std::map<SpectralFieldId, Framework::Selector::ComplexScalarField*> fields; 

   /**
    * @brief Vector of equation descriptors
    */
   std::vector<MinEqInfo> eqInfos;

   /**
    * @brief Map of Real explicit terms
    */
   std::map<std::size_t, std::map<SpectralFieldId, std::map<SpectralFieldId, std::vector<SparseMatrix>>>> exDTerm;

   /**
    * @brief Map of complex explicit terms
    */
   std::map<std::size_t, std::map<SpectralFieldId, std::map<SpectralFieldId, std::vector<SparseMatrixZ>>>> exZTerm;

   /**
    * @brief Galerkin stencils
    */
   std::map<SpectralFieldId, std::vector<SparseMatrix>> stencils;

   /**
    * @brief Spectral constraint kernels
    */
   std::map<SpectralFieldId, std::shared_ptr<Spectral::Kernel::ISpectralKernel>> constraints;

   /**
    * @brief Spectral boundary value kernels
    */
   std::map<SpectralFieldId, std::shared_ptr<Spectral::Kernel::ISpectralKernel>> bcvalues;

   /**
    * @brief Spectral source kernels
    */
   std::map<SpectralFieldId, std::shared_ptr<Spectral::Kernel::ISpectralKernel>> sources;

   /**
    * @brief Solution updaters
    */
   std::map<SpectralFieldId, std::shared_ptr<Equations::SolutionUpdater>> solups;
};


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
