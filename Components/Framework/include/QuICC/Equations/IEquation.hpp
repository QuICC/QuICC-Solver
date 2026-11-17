/**
 * @file IEquation.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_IEQUATION_HPP
#define QUICC_EQUATIONS_IEQUATION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Arithmetics/Utility.hpp"
#include "Types/Typedefs.hpp"
#include "Arithmetics/Utility.hpp"
#include "Arithmetics/LinearAlgebra.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/EquationData.hpp"
#include "QuICC/Equations/SolutionUpdater.hpp"
#include "QuICC/Variables/VariableRequirement.hpp"
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformPath.hpp"
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"

namespace QuICC {

namespace Equations {

   /// Forward declaration of IEquation
   class IEquation;

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IEquation : public EquationData
   {
      public:
         /// Typedef for the the spectral field component ID iterator
         typedef std::vector<FieldComponents::Spectral::Id>::const_iterator   SpectralComponent_iterator;

         /// Typedef for the the spectral field component ID iterator range
         typedef std::pair<SpectralComponent_iterator,SpectralComponent_iterator>  SpectralComponent_range;

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEquation() = default;

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const = 0;

         /**
          * @brief Access the shared resolution
          */
         virtual const Resolution& res() const = 0;

         /**
          * @brief Get the number of spectral components
          */
         virtual int nSpectral() const = 0;

         /**
          * @brief Get vector spectral component range
          */
         virtual SpectralComponent_range spectralRange() const = 0;

         /**
          * @brief Initialise the equation
          */
         virtual void init(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Initialise the nonlinear kernel
          */
         virtual void initNLKernel(const bool force = false);

         /**
          * @brief Get forward transform paths
          */
         virtual std::vector<Transform::TransformPath> forwardPaths();

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() = 0;

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const; // = 0;

         /**
          * @brief Get nonlinear kernel
          */
         Physical::Kernel::SharedIPhysicalKernel spNLKernel() const;

         /**
          * @brief Initialise the spectral equation matrices
          */
         virtual void initSpectralMatrices() = 0;

         /**
          * @brief Set spectral constraint kernel
          */
         virtual void setConstraintKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel);

         /**
          * @brief Get spectral constraint kernel
          */
         Spectral::Kernel::SharedISpectralKernel spConstraintKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Initialize constraint kernels
          *
          * @param spMesh  Physical mesh
          */
         virtual void initConstraintKernel(const std::shared_ptr<std::vector<Array> > spMesh);

         /**
          * @brief Set spectral source kernel
          */
         virtual void setSrcKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel);

         /**
          * @brief Get spectral source kernel
          */
         Spectral::Kernel::SharedISpectralKernel spSrcKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Initialize source kernels
          */
         virtual void initSrcKernel();

         /**
          * @brief Link with other equation
          *
          * @param spEq Equation to link to
          */
         virtual void linkEquation(std::shared_ptr<IEquation> spEq);

         /**
          * @brief Write equation diagnostic
          *
          * @param isAsciiTime  Is ASCII writing time?
          * @param isHdf5Time   Is HDF5 writing time?
          */
         virtual void writeDiagnostics(const bool isAsciiTime, const bool isHdf5Time) const;

         /**
          * @brief Update constraint kernels
          *
          * @param time       Simulation time
          * @param timestep   Simulation timestep
          * @param isFinished Full timestep was computed?
          */
         virtual void updateConstraintKernel(const MHDFloat time, const MHDFloat timestep, const bool isFinished);

      protected:
         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

         /**
          * @brief Initialize source kernels
          */
         virtual void initSolutionUpdater();

         /**
          * @brief Set the default nonlinear components
          */
         virtual void setNLComponents() = 0;

         /**
          * @brief Transform steps object
          */
         virtual std::shared_ptr<Transform::ITransformSteps> transformSteps() const;

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<bool> disabledBackwardPaths() const = 0;

         /**
          * @brief Initialise the spectral equation matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Shared physical interaction kernel
          */
         Physical::Kernel::SharedIPhysicalKernel mspNLKernel;

         /**
          * @brief Solution updaters
          */
          std::map<FieldComponents::Spectral::Id,std::shared_ptr<SolutionUpdater>> mSolUps;

         /**
          * @brief Shared spectral source kernel for each component
          */
          std::map<FieldComponents::Spectral::Id,Spectral::Kernel::SharedISpectralKernel> mSrcKernel;

         /**
          * @brief Shared spectral boundary value kernel for each component
          */
          std::map<FieldComponents::Spectral::Id,Spectral::Kernel::SharedISpectralKernel> mBoundaryKernel;

         /**
          * @brief Shared spectral constraint kernel for each component
          */
          std::map<FieldComponents::Spectral::Id,Spectral::Kernel::SharedISpectralKernel> mConstraintKernel;

      private:

         /**
          * @brief Initialise the Galerkin stencisl
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initGalerkinStencils(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the quasi-inverse spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initQIMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the explicit spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          * @param opId    Type of explicit operator
          */
         void initExplicitMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId, const std::size_t opId);

         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const; // = 0;

         /**
          * @brief Set the explicit matrix operator
          */
         virtual void setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const; // = 0;

   };

   /// Typedef for a smart IEquation
   typedef std::shared_ptr<IEquation> SharedIEquation;

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IEQUATION_HPP
