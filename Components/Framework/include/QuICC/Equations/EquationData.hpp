/**
 * @file EquationData.hpp
 * @brief Lowest building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_EQUATIONDATA_HPP
#define QUICC_EQUATIONS_EQUATIONDATA_HPP

// System includes
//
#include <limits>
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/EquationOptions.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Variables/VariableRequirement.hpp"
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/Model/IModelBackend.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Lowest building block for the implementation of an equation
    */
   class EquationData
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit EquationData(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit EquationData(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~EquationData() = default;

         /**
          * @brief Get name ID of the unknown
          */
         std::size_t name() const;

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         void setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField);

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         void setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Get the galerkin stencil matrix
          *
          * @param compId  Field component ID
          * @param j       Matrix index
          */
         const SparseMatrix& galerkinStencil(const FieldComponents::Spectral::Id compId, const int j) const;

         /**
          * @brief Check if real quasi inverse matrices exist
          *
          * @param compId  Field component ID
          */
         bool hasQID(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Check if complex quasi inverse matrices exist
          *
          * @param compId  Field component ID
          */
         bool hasQIZ(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Check if real explicit matrices exist
          *
          * @param opId    Operator ID
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitDTerm(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Check if complex explicit matrices exist
          *
          * @param opId    Operator ID
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitZTerm(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Get the quasi inverse matrices
          *
          * @param compId  Field component ID
          * @param j       Matrix index
          */
         template <typename TOperator> const TOperator& quasiInverse(const FieldComponents::Spectral::Id compId, const int j) const;

         /**
          * @brief Get the explicit matrices
          *
          * @param opId    Operator ID
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          * @param j       Matrix index
          */
         template <typename TOperator> const TOperator& explicitOperator(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const;

         /**
          * @brief Get the coupling information
          *
          * @param compId  Field component ID
          */
         const CouplingInformation& couplingInfo(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Get map of field storage requirements information
          *
          * \mhdBug Ultimatively this should depend on component
          */
         const VariableRequirement& requirements() const;

         /**
          * @brief Get map of imposed field storage requirements information
          *
          * \mhdBug Ultimatively this should depend on component
          */
         const VariableRequirement& imposedRequirements() const;

         /**
          * @brief Get map of field storage requirements information
          */
         const FieldRequirement& requirements(std::size_t id) const;

         /**
          * @brief Get map of imposed field storage requirements information
          */
         const FieldRequirement& imposedRequirements(std::size_t id) const;

         /**
          * @brief Get the equation parameters
          */
         const EquationParameters& eqParams() const;

         /**
          * @brief Get the equation parameters
          */
         const SharedCEquationParameters spEqParams() const;

         /**
          * @brief Get the list of boundary conditions
          */
         const SimulationBoundary& bcIds() const;

         /**
          * @brief Set the solver index for the coupling information
          *
          * @param compId Spectral component IDa
          * @param idx    Matrix index
          */
         void setSolverIndex(const FieldComponents::Spectral::Id compId, const int idx);

         /**
          * @brief Timing of the solver for the equation
          */
         std::size_t  solveTiming() const;

         /**
          * @brief Current simulation time to allow for timedependent implementations
          */
         MHDFloat  time() const;

         /**
          * @brief Set current simulation time to allow for timedependent implementations
          *
          * @param time       Current simulation time
          * @param finished   Flag for completed multistage timestep
          */
         virtual void  setTime(const MHDFloat time, const bool finished);

         /**
          * @brief Get the nonlinear integration components order
          */
         const std::vector<std::pair<FieldComponents::Spectral::Id,std::size_t> >& nlComponents() const;

         /**
          * @brief Get equation type
          */
         CouplingInformation::EquationTypeId equationType() const;

         /**
          * @brief Cleanup the backend
          */
         void cleanupBackend();

         /**
          * @brief Get equation options
          */
         const EquationOptions& options() const;

      protected:
         enum ForwardPathsId
         {
            FWD_IS_CUSTOM = 0,
            FWD_IS_FIELD = 1,
            FWD_IS_NONLINEAR = 2,
         };

         /**
          * @brief Set the unknown name of equation
          */
         void setForwardPathsType(const ForwardPathsId id);

         /**
          * @brief Set the unknown name of equation
          */
         void setName(std::size_t name);

         /**
          * @brief Set solver timing
          */
         void setSolveTiming(const std::size_t timeId);

         /**
          * @brief Get spatial scheme
          */
         const SpatialScheme::ISpatialScheme& ss() const;

         /**
          * @brief Get model backend
          */
         const Model::IModelBackend& backend() const;

         /**
          * @brief Get equation options pointer
          */
         std::shared_ptr<EquationOptions> spOptions() const;

         /**
          * @brief Add a nonlinear integration component
          */
         void addNLComponent(const FieldComponents::Spectral::Id compId, const std::size_t flag);

         /**
          * @brief Update field requirements information
          */
         FieldRequirement& updateFieldRequirements(std::size_t id);

         /**
          * @brief Set map of component and explicit matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> >& rEDMatrices(const std::size_t opId);

         /**
          * @brief Set map of component and explicit matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> >& rEZMatrices(const std::size_t opId);

         /**
          * @brief Get shared scalar variable
          *
          * @param name Physical name of the field
          */
         Framework::Selector::VariantSharedScalarVariable spScalar(std::size_t name) const;

         /**
          * @brief Get shared vector variable
          *
          * @param name Physical name of the field
          */
         Framework::Selector::VariantSharedVectorVariable spVector(std::size_t name) const;

         /**
          * @brief Type of forward transform
          */
         ForwardPathsId mForwardPathsType;

         /**
          * @brief Storage for the variable requirements
          *
          * \mhdBug Ultimatively this should depend on component
          */
         VariableRequirement mRequirements;

         /**
          * @brief Storage for the imposed variable requirements
          *
          * \mhdBug Ultimatively this should depend on component
          */
         VariableRequirement mImposedRequirements;

         /**
          * @brief Coupling information of the equation
          */
         std::map<FieldComponents::Spectral::Id, CouplingInformation>  mCouplingInfos;

         /**
          * @brief Map of component and quasi inverse matrices (real operators)
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mQIDMatrices;

         /**
          * @brief Map of component and quasi inverse matrices (complex operators)
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrixZ> > mQIZMatrices;

         /**
          * @brief Map of component and galerkin stencil matrices
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mGStencils;

         /**
          * @brief Storage for the shared boundary condition list
          */
         SharedSimulationBoundary mspBcIds;

         /**
          * @brief Storage for the solve timing
          */
         std::size_t   mSolveTiming;

         /**
          * @brief Nonlinear integration component order
          */
         std::vector<std::pair<FieldComponents::Spectral::Id,std::size_t> >   mNLComponents;

      private:
         /**
          * @brief Storage for smart equation parameters
          */
         SharedEquationParameters   mspEqParams;

         /**
          * @brief Storage for smart equation parameters
          */
         SpatialScheme::SharedCISpatialScheme   mspSpatialScheme;

         /**
          * @brief Model generator
          */
         std::shared_ptr<Model::IModelBackend> mspBackend;

         /**
          * @brief Equation options
          */
         std::shared_ptr<EquationOptions> mspOptions;

         /**
          * @brief Name ID of the unknown
          */
         MHDFloat mTime;

         /**
          * @brief Map of name and pointer for the scalar variables
          */
         std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>  mScalars;

         /**
          * @brief Map of name and pointer for the vector variables
          */
         std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>  mVectors;

         /**
          * @brief Name ID of the unknown
          */
         std::size_t mName;

         /**
          * @brief Map of component and explicit linear matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mELDMatrices;

         /**
          * @brief Map of component and explicit linear matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mELZMatrices;

         /**
          * @brief Map of component and explicit nonlinear matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mENLDMatrices;

         /**
          * @brief Map of component and explicit nonlinear matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mENLZMatrices;

         /**
          * @brief Map of component and explicit nextstep matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mENSDMatrices;

         /**
          * @brief Map of component and explicit nextstep matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mENSZMatrices;
   };

   /// Typedef for a smart EquationData
   typedef std::shared_ptr<EquationData> SharedEquationData;

   /**
    * @brief Update time average
    */
   template <typename TData> TData incrementTimeAverage(const TData avg, const TData newData, const MHDFloat time, const MHDFloat timestep);
   MHDFloat incrementTimeAverage(const MHDComplex avg, const MHDFloat newData, const MHDFloat time, const MHDFloat timestep);

   /**
    * @brief Don't update time average
    */
   template <typename TData> TData noupdateTimeAverage(const TData avg, const TData newData);
   MHDFloat noupdateTimeAverage(const MHDComplex avg, const MHDFloat newData);


   template <typename TData> TData incrementTimeAverage(const TData avg, const TData newData, const MHDFloat time, const MHDFloat timestep)
   {
      MHDFloat stepWeight = timestep/time;
      MHDFloat avgWeight = (time-timestep)/time;

      return avgWeight*avg + stepWeight*newData;
   }

   template <typename TData> TData noupdateTimeAverage(const TData avg, const TData)
   {
      return avg;
   }
} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EQUATIONDATA_HPP
