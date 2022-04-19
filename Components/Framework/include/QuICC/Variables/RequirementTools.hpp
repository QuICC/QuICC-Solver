/**
 * @file RequirementTools.hpp
 * @brief Implementation of requirement tools for equations and variables
 */

#ifndef QUICC_REQUIREMENTTOOLS_HPP
#define QUICC_REQUIREMENTTOOLS_HPP

// System includes
//
#include <map>
#include <cassert>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Variables/VariableRequirement.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

namespace QuICC {

   /**
    * @brief Implementation of requirement tools for equations and variables
    */
   class RequirementTools
   {
      public:
         /**
          * @brief Merge variable requirements from equations
          *
          * @param scalarEqs     Scalar equations
          * @param vectorEqs     Vector equations
          */
         static VariableRequirement mergeRequirements(const std::vector<std::shared_ptr<Equations::IScalarEquation> >& scalarEqs, const std::vector<std::shared_ptr<Equations::IVectorEquation> >& vectorEqs);

         /**
          * @brief Merge imposed variable requirements from equations
          *
          * @param scalarEqs     Scalar equations
          * @param vectorEqs     Vector equations
          */
         static VariableRequirement mergeImposedRequirements(const std::vector<std::shared_ptr<Equations::IScalarEquation> >& scalarEqs, const std::vector<std::shared_ptr<Equations::IVectorEquation> >& vectorEqs);

         /**
          * @brief Initialise variables and variable requirements from equations
          *
          * @param rScalarVars   Scalar variables
          * @param rVectorVars   Vector variables
          * @param varInfo       Variable requirements
          * @param spRes      Shared resolution
          */
         static void initVariables(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& rScalarVars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& rVectorVars, const VariableRequirement& varInfo, SharedResolution spRes);

         /**
          * @brief Map variables to the corresponding equation
          *
          * @param rScalarEqs          Scalar equations
          * @param rVectorEqs          Vector equations
          * @param scalarVars          Scalar variables
          * @param vectorVars          Vector variables
          */
         static void mapEquationVariables(std::vector<std::shared_ptr<Equations::IScalarEquation> >& rScalarEqs, std::vector<std::shared_ptr<Equations::IVectorEquation> >& rVectorEqs, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars, const SharedSimulationBoundary spBcs, std::vector<std::size_t>& rUnmapped);

         /**
          * @brief Map imposed variables to the corresponding equation
          *
          * @param rScalarEqs          Scalar equations
          * @param rVectorEqs          Vector equations
          * @param scalarVars          Scalar variables
          * @param vectorVars          Vector variables
          */
         static void mapImposedVariables(std::vector<std::shared_ptr<Equations::IScalarEquation> >& rScalarEqs, std::vector<std::shared_ptr<Equations::IVectorEquation> >& rVectorEqs, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars, std::vector<std::size_t>& rUnmapped);

         /**
          * @brief Build backward transform tree
          *
          * @param backwardTree Transform tree for backward projection
          * @param varInfo      Variable requirements
          */
         static void buildBackwardTree(std::vector<Transform::TransformTree>& backwardTree, const std::vector<std::shared_ptr<Equations::IScalarEquation> >& scalarEqs, const std::vector<std::shared_ptr<Equations::IVectorEquation> >& vectorEqs);


         /**
          * @brief Build forward transform tree
          *
          * @param forwardTree Transform tree for forward integration
          * @param rScalarEqs          Scalar equations
          * @param rVectorEqs          Vector equations
          */
         static void buildForwardTree(std::vector<Transform::TransformTree>& forwardTree, std::vector<std::shared_ptr<Equations::IScalarEquation> >& rScalarEqs, std::vector<std::shared_ptr<Equations::IVectorEquation> >& rVectorEqs);

         /**
          * @brief Build backward transform tree
          *
          * @param backwardTree Transform tree for backward projection
          * @param varInfo      Variable requirements
          */
         static void buildBackwardTree(std::vector<Transform::TransformTree>& backwardTree, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars);

         /**
          * @brief Get backward transform paths
          */
         static std::vector<Transform::TransformPath> backwardPaths(Framework::Selector::VariantSharedScalarVariable spScalar);

         /**
          * @brief Get backward transform paths
          */
         static std::vector<Transform::TransformPath> backwardPaths(Framework::Selector::VariantSharedVectorVariable spVector);

      protected:

      private:
         /**
          * @brief Constructor
          */
         RequirementTools() = default;

         /**
          * @brief Destructor
          */
         ~RequirementTools() = default;
   };

}

#endif // QUICC_REQUIREMENTTOOLS_HPP
