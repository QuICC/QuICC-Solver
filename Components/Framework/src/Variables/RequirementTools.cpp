/**
 * @file RequirementTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 */

// System includes
//

// Project includes
//
#include "QuICC/Variables/RequirementTools.hpp"

namespace QuICC {

   void RequirementTools::mergeRequirements(VariableRequirement& req, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs)
   {
      // Loop over all scalar equations
      for(auto scalEqIt = scalarEqs.cbegin(); scalEqIt < scalarEqs.cend(); scalEqIt++)
      {
         req.merge((*scalEqIt)->requirements());
      }

      // Loop over all vector equations
      for(auto vectEqIt = vectorEqs.cbegin(); vectEqIt < vectorEqs.cend(); vectEqIt++)
      {
         req.merge((*vectEqIt)->requirements());
      }
   }

   void RequirementTools::mergeImposedRequirements(VariableRequirement& req, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs)
   {
      // Loop over all scalar equations
      for(auto scalEqIt = scalarEqs.cbegin(); scalEqIt < scalarEqs.cend(); scalEqIt++)
      {
         req.merge((*scalEqIt)->imposedRequirements());
      }

      // Loop over all vector equations
      for(auto vectEqIt = vectorEqs.cbegin(); vectEqIt < vectorEqs.cend(); vectEqIt++)
      {
         req.merge((*vectEqIt)->imposedRequirements());
      }
   }

   void RequirementTools::initVariables(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& rScalarVars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& rVectorVars, const VariableRequirement& varInfo, SharedResolution spRes)
   {
      //
      // Create the required variables
      //

      // Initialise variables
      for(auto infoIt = varInfo.cbegin(); infoIt != varInfo.cend(); infoIt++)
      {
         // Check if spectral variable is required
         if(infoIt->second.needSpectral())
         {
            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Create the shared scalar variable
               auto spScalar = spRes->sim().ss().createSVar(spRes);
               rScalarVars.insert(std::make_pair(infoIt->first, spScalar));
               std::visit([&](auto&& p){p->initSpectral(infoIt->second.spectralIds());}, rScalarVars.find(infoIt->first)->second);
            } else
            {
               // Create the shared vector variable
               auto spVector = spRes->sim().ss().createVVar(spRes);
               rVectorVars.insert(std::make_pair(infoIt->first, spVector));
               std::visit([&](auto&& p){p->initSpectral(infoIt->second.spectralIds());}, rVectorVars.find(infoIt->first)->second);
            }

            // Initialise the physical values if required
            if(infoIt->second.needPhysical())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  std::visit([&](auto&& p){p->initPhysical(infoIt->second.mapPhysicalComps());}, rScalarVars.at(infoIt->first));
               } else
               {
                  std::visit([&](auto&& p){p->initPhysical(infoIt->second.mapPhysicalComps());}, rVectorVars.at(infoIt->first));
               }
            }

            // Initialise the physical gradient values if required
            if(infoIt->second.needPhysicalGradient())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  std::visit([&](auto&& p){p->initPhysicalGradient(FieldComponents::Spectral::SCALAR, infoIt->second.mapGradientComps(FieldComponents::Spectral::SCALAR));}, rScalarVars.at(infoIt->first));
               } else
               {
                  for(auto it = infoIt->second.spectralIds().cbegin(); it != infoIt->second.spectralIds().cend(); ++it)
                  {
                     std::visit([&](auto&& p){p->initPhysicalGradient(*it, infoIt->second.mapGradientComps(*it));}, rVectorVars.at(infoIt->first));
                  }
               }
            }

            // Initialise the physical 2nd order gradient values if required
            if(infoIt->second.needPhysicalGradient2())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  std::visit([&](auto&& p){p->initPhysicalGradient2(FieldComponents::Spectral::SCALAR, infoIt->second.mapGradient2Comps(FieldComponents::Spectral::SCALAR));}, rScalarVars.at(infoIt->first));
               } else
               {
                  throw std::logic_error("2nd order Vector gradient is not implemented!");
               }
            }

            // Initialise the physical curl values if required
            if(infoIt->second.needPhysicalCurl())
            {
               // Separate scalar and vector fields
               if(infoIt->second.isScalar())
               {
                  throw std::logic_error("Can't initialise curl on scalar field!");
               } else
               {
                  std::visit([&](auto&& p){p->initPhysicalCurl(infoIt->second.mapCurlComps());}, rVectorVars.at(infoIt->first));
               }
            }

            // Separate scalar and vector fields
            if(infoIt->second.isScalar())
            {
               // Initialise to zero
               std::visit([&](auto&& p){p->setZeros();}, rScalarVars.at(infoIt->first));

               std::visit([&](auto&& p){p->profileStorage();}, rScalarVars.at(infoIt->first));

            } else
            {
               // Initialise to zero
               std::visit([&](auto&& p){p->setZeros();}, rVectorVars.at(infoIt->first));

               std::visit([&](auto&& p){p->profileStorage();}, rVectorVars.at(infoIt->first));

            }
         }
      }
   }

   void RequirementTools::mapEquationVariables(std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars, const SharedSimulationBoundary spBcs, std::vector<std::size_t>& rUnmapped)
   {
      // Loop over all scalar variables
      for(auto scalIt = scalarVars.cbegin(); scalIt != scalarVars.cend(); scalIt++)
      {
         auto varName = scalIt->first;

         // Loop over scalar equations
         for(auto scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set scalar variable as unknown scalar field
            if((*scalEqIt)->name() == varName)
            {
               (*scalEqIt)->setUnknown(scalarVars.at(varName));

               // Finish initialisation of equation
               (*scalEqIt)->init(spBcs);
            }

            // Set scalar variable as additional scalar field
            if((*scalEqIt)->requirements(varName).needPhysical() || (*scalEqIt)->requirements(varName).needPhysicalGradient() || (*scalEqIt)->requirements(varName).needPhysicalGradient2())
            {
               (*scalEqIt)->setField(varName, scalarVars.at(varName));
            }
         }

         // Loop over vector equations
         for(auto vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->requirements(varName).needPhysical() || (*vectEqIt)->requirements(varName).needPhysicalGradient() || (*vectEqIt)->requirements(varName).needPhysicalGradient2())
            {
               (*vectEqIt)->setField(varName, scalarVars.at(varName));
            }
         }
      }

      // Loop over all vector variables
      for(auto vectIt = vectorVars.cbegin(); vectIt != vectorVars.cend(); vectIt++)
      {
         auto varName = vectIt->first;

         // Loop over scalar equations
         for(auto scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->requirements(varName).needPhysical() || (*scalEqIt)->requirements(varName).needPhysicalGradient() || (*scalEqIt)->requirements(varName).needPhysicalCurl() || (*scalEqIt)->requirements(varName).needPhysicalGradient2())
            {
               (*scalEqIt)->setField(varName, vectorVars.at(varName));
            }
         }

         // Loop over vector equations
         for(auto vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set vector variable as unknown vector field
            if((*vectEqIt)->name() == varName)
            {
               (*vectEqIt)->setUnknown(vectorVars.at(varName));

               // Finish initialisation of equation
               (*vectEqIt)->init(spBcs);
            }

            // Set vector variable as additional vector field
            if((*vectEqIt)->requirements(varName).needPhysical() || (*vectEqIt)->requirements(varName).needPhysicalGradient() || (*vectEqIt)->requirements(varName).needPhysicalCurl() || (*vectEqIt)->requirements(varName).needPhysicalGradient2())
            {
               (*vectEqIt)->setField(varName, vectorVars.at(varName));
            }
         }
      }
   }

   void RequirementTools::buildBackwardTree(std::vector<Transform::TransformTree>& backwardTree, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs)
   {
      std::map<std::size_t, std::vector<Transform::TransformPath> > branches;

      // Loop over scalar equations
      for(auto scalEqIt = scalarEqs.begin(); scalEqIt < scalarEqs.end(); scalEqIt++)
      {
         auto eqBranches = (*scalEqIt)->backwardPaths();
         if(eqBranches.size() > 0)
         {
            auto eqName = (*scalEqIt)->name();
            branches.insert(std::make_pair(eqName, std::vector<Transform::TransformPath>()));
            branches.find(eqName)->second.insert(branches.find(eqName)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Loop over vector equations
      for(auto vectEqIt = vectorEqs.begin(); vectEqIt < vectorEqs.end(); vectEqIt++)
      {
         auto eqBranches = (*vectEqIt)->backwardPaths();
         if(eqBranches.size() > 0)
         {
            auto eqName = (*vectEqIt)->name();
            branches.insert(std::make_pair(eqName, std::vector<Transform::TransformPath>()));
            branches.find(eqName)->second.insert(branches.find(eqName)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Create the backward transform tree(s)
      Transform::TransformTreeTools::generateTrees(backwardTree, branches, TransformDirection::BACKWARD);
   }

   void RequirementTools::buildForwardTree(std::vector<Transform::TransformTree>& forwardTree, std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs)
   {
      std::map<std::size_t, std::vector<Transform::TransformPath> > branches;

      // Loop over scalar equations
      for(auto scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
      {
         auto eqBranches = (*scalEqIt)->forwardPaths();
         if(eqBranches.size() > 0)
         {
            auto eqName = (*scalEqIt)->name();
            branches.insert(std::make_pair(eqName, std::vector<Transform::TransformPath>()));
            branches.find(eqName)->second.insert(branches.find(eqName)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Loop over vector equations
      for(auto vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
      {
         auto eqBranches = (*vectEqIt)->forwardPaths();
         if(eqBranches.size() > 0)
         {
            auto eqName = (*vectEqIt)->name();
            branches.insert(std::make_pair(eqName, std::vector<Transform::TransformPath>()));
            branches.find(eqName)->second.insert(branches.find(eqName)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Create the integrator tree(s)
      Transform::TransformTreeTools::generateTrees(forwardTree, branches, TransformDirection::FORWARD);
   }

   void RequirementTools::mapImposedVariables(std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars, std::vector<std::size_t>& rUnmapped)
   {
      // Loop over all scalar variables
      for(auto scalIt = scalarVars.cbegin(); scalIt != scalarVars.cend(); scalIt++)
      {
         // Loop over scalar equations
         for(auto scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*scalEqIt)->imposedRequirements(scalIt->first).needPhysical() || (*scalEqIt)->imposedRequirements(scalIt->first).needPhysicalGradient() || (*scalEqIt)->imposedRequirements(scalIt->first).needPhysicalGradient2())
            {
               (*scalEqIt)->setField(scalIt->first, scalarVars.at(scalIt->first));
            }
         }

         // Loop over vector equations
         for(auto vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set scalar variable as additional scalar field
            if((*vectEqIt)->imposedRequirements(scalIt->first).needPhysical() || (*vectEqIt)->imposedRequirements(scalIt->first).needPhysicalGradient() || (*vectEqIt)->imposedRequirements(scalIt->first).needPhysicalGradient2())
            {
               (*vectEqIt)->setField(scalIt->first, scalarVars.at(scalIt->first));
            }
         }
      }

      // Loop over all vector variables
      for(auto vectIt = vectorVars.cbegin(); vectIt != vectorVars.cend(); vectIt++)
      {
         // Loop over scalar equations
         for(auto scalEqIt = rScalarEqs.begin(); scalEqIt < rScalarEqs.end(); scalEqIt++)
         {
            // Set vector variable as additional vector field
            if((*scalEqIt)->imposedRequirements(vectIt->first).needPhysical() || (*scalEqIt)->imposedRequirements(vectIt->first).needPhysicalGradient() || (*scalEqIt)->imposedRequirements(vectIt->first).needPhysicalCurl() || (*scalEqIt)->imposedRequirements(vectIt->first).needPhysicalGradient2())
            {
               (*scalEqIt)->setField(vectIt->first, vectorVars.at(vectIt->first));
            }
         }

         // Loop over vector equations
         for(auto vectEqIt = rVectorEqs.begin(); vectEqIt < rVectorEqs.end(); vectEqIt++)
         {
            // Set vector variable as additional vector field
            if((*vectEqIt)->imposedRequirements(vectIt->first).needPhysical() || (*vectEqIt)->imposedRequirements(vectIt->first).needPhysicalGradient() || (*vectEqIt)->imposedRequirements(vectIt->first).needPhysicalCurl() || (*vectEqIt)->imposedRequirements(vectIt->first).needPhysicalGradient2())
            {
               (*vectEqIt)->setField(vectIt->first, vectorVars.at(vectIt->first));
            }
         }
      }
   }

   void RequirementTools::buildBackwardTree(std::vector<Transform::TransformTree>& backwardTree, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalarVars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectorVars)
   {
      std::map<std::size_t, std::vector<Transform::TransformPath> > branches;

      // Loop over scalar equations
      for(auto scalIt = scalarVars.begin(); scalIt != scalarVars.end(); scalIt++)
      {
         auto eqBranches = RequirementTools::backwardPaths(scalIt->second);
         if(eqBranches.size() > 0)
         {
            branches.insert(std::make_pair(scalIt->first, std::vector<Transform::TransformPath>()));
            branches.find(scalIt->first)->second.insert(branches.find(scalIt->first)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Loop over vector equations
      for(auto vectIt = vectorVars.begin(); vectIt != vectorVars.end(); vectIt++)
      {
         auto eqBranches = RequirementTools::backwardPaths(vectIt->second);
         if(eqBranches.size() > 0)
         {
            branches.insert(std::make_pair(vectIt->first, std::vector<Transform::TransformPath>()));
            branches.find(vectIt->first)->second.insert(branches.find(vectIt->first)->second.end(),eqBranches.begin(), eqBranches.end());
         }
      }

      // Create the backward transform tree(s)
      Transform::TransformTreeTools::generateTrees(backwardTree, branches, TransformDirection::BACKWARD);
   }

   std::vector<Transform::TransformPath> RequirementTools::backwardPaths(Framework::Selector::VariantSharedScalarVariable spScalar)
   {
      std::vector<Transform::TransformPath> paths;

      std::shared_ptr<Transform::ITransformSteps>  spSteps;
      std::visit([&](auto&& p){spSteps = Transform::createTransformSteps(p->dom(0).res().sim().spSpatialScheme());}, spScalar);

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, spScalar))
      {
         std::map<FieldComponents::Physical::Id,bool> compsMap;
         compsMap.insert(std::make_pair(FieldComponents::Physical::SCALAR, true));
         auto b = spSteps->backwardScalar(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, spScalar))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).grad().enabled());}, spScalar);
         auto b = spSteps->backwardGradient(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, spScalar))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>{return (p->dom(0).grad2().enabled());}, spScalar);
         auto b = spSteps->backwardGradient2(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      return paths;
   }


   std::vector<Transform::TransformPath> RequirementTools::backwardPaths(Framework::Selector::VariantSharedVectorVariable spVector)
   {
      std::vector<Transform::TransformPath> paths;

      std::shared_ptr<Transform::ITransformSteps>  spSteps;
      std::visit([&](auto&& p){spSteps = Transform::createTransformSteps(p->dom(0).res().sim().spSpatialScheme());}, spVector);

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, spVector))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).phys().enabled());}, spVector);
         auto branches = spSteps->backwardVector(compsMap);
         paths.insert(paths.end(), branches.begin(), branches.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, spVector))
      {
         auto specMap = std::visit([&](auto&& p)->std::map<FieldComponents::Spectral::Id,bool>{return (p->dom(0).perturbation().enabled());}, spVector);
         for(auto it = specMap.begin(); it != specMap.end(); ++it)
         {
            if(it->second)
            {
               auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).grad(it->first).enabled());}, spVector);
               auto b = spSteps->backwardVGradient(it->first, compsMap);
               paths.insert(paths.end(), b.begin(), b.end());
            }
         }
      }

// Not yet implemented
//      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, spVector))
//      {
//         auto range = this->spectralRange();
//         for(auto it = range.first; it != range.second; ++it)
//         {
//            auto compsMap = std::visit([&](auto&& p)->std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>{return (p->dom(0).grad2(*it).enabled());}, spVector);
//            auto b = spSteps->backwardVGradient2(*it, compsMap);
//            paths.insert(paths.end(),b.begin(), b.end());
//         }
//      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasCurl());}, spVector))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).curl().enabled());}, spVector);
         auto b = spSteps->backwardCurl(compsMap);
         paths.insert(paths.end(),b.begin(), b.end());
      }

      return paths;
   }

} // QuICC
