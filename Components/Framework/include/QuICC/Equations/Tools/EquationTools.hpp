/**
 * @file EquationTools.hpp
 * @brief Implementation of equation tools
 */

#ifndef QUICC_EQUATIONS_TOOLS_EQUATIONTOOLS_HPP
#define QUICC_EQUATIONS_TOOLS_EQUATIONTOOLS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Equations/Tools/EquationConditions.hpp"
#include "QuICC/Equations/Tools/EquationSorters.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#endif //QUICC_DEBUG

namespace QuICC {

namespace Equations {

/**
 * @brief Implementation of equation tools
 */
namespace Tools {

   /**
    * @brief
    *
    * @tparam T                  Template for scalar or vector equation
    * @param rEqs                Vector of shared equations
    * @param rPrognosticRange    Range of prognostic equations
    * @param rDiagnosticRange    Range of diagnostic equations
    * @param rTrivialRange       Range of trivial equations
    * @param rWrapperRange       Range of wrapper equations
    */
   template <typename T> void sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rWrapperRange);

   /**
    * @brief Identify the solver indexes and set it on the equations
    *
    * @param scalarRange   Range of scalar equations of same type
    * @param vectorRange   Range of vector equations of same type
    *
    * \mhdBug This implementation is very BAD! Needs a cleanup as soon as possible
    */
   template <typename TScalar, typename TVector> void identifySolver(const std::pair<TScalar,TScalar>& scalarRange, const std::pair<TVector,TVector>& vectorRange);

   /**
    * @brief Create a small sparse matrix to select matrix block through Kronecker product
    *
    * @param nBlocks Number of blocks (ie. nBlocks x nBlocks matrix)
    * @param row     Row of the required block
    * @param col     Column of the required block
    */
   SparseMatrix makeBlockMatrix(const int nBlocks, const int row, const int col);

   /**
    * @brief Get nonlinear kernels from equations
    */
   template <typename TIt> void getNonlinearKernels(std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& rKernels, TIt first, TIt last);

   /**
    * @brief Setup physical kernels
    */
   void setupPhysicalKernels(std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& rKernels, const std::vector<Array>& mesh);

//
// Implementation follows
//

   template <typename T> void sortByType(typename std::vector<T>& rEqs, std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator>& rPrognosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rDiagnosticRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rTrivialRange, std::pair<typename std::vector<T>::iterator,typename std::vector<T>::iterator>& rWrapperRange)
   {
      // Sort equations
      std::stable_sort(rEqs.begin(), rEqs.end(), Sorters::EquationType());

      // Initialise empty scalar ranges
      rPrognosticRange = std::make_pair(rEqs.end(), rEqs.end());
      rDiagnosticRange = std::make_pair(rEqs.end(), rEqs.end());
      rTrivialRange = std::make_pair(rEqs.end(), rEqs.end());

      // Determine the ranges for the different types
      typename std::vector<T>::iterator  eqIt;

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsPrognostic());
      if(eqIt != rEqs.end())
      {
         rPrognosticRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsDiagnostic());
      if(eqIt != rEqs.end())
      {
         rDiagnosticRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsTrivial());
      if(eqIt != rEqs.end())
      {
         rTrivialRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }

      eqIt = std::find_if(rEqs.begin(), rEqs.end(), Conditions::IsWrapper());
      if(eqIt != rEqs.end())
      {
         rWrapperRange = std::equal_range(rEqs.begin(), rEqs.end(), *eqIt, Sorters::EquationType());
      }
   }

   template <typename TScalar, typename TVector> void identifySolver(const std::pair<TScalar,TScalar>& scalarRange, const std::pair<TVector,TVector>& vectorRange)
   {
      // Current solver indexes for real and complex solvers
      int dIdx = 0;
      int zIdx = 0;

      // Coupling flag
      bool coupled = false;

      // Field identification
      FieldComponents::Spectral::Id compId = FieldComponents::Spectral::SCALAR;
      SpectralFieldId   fieldId;

      // Loop over the scalar equations
      for(auto scalEqIt = scalarRange.first; scalEqIt != scalarRange.second; ++scalEqIt)
      {
         // Build field identity
         fieldId = std::make_pair((*scalEqIt)->name(), compId);

         // Loop over already identified equations
         for(auto& doneSEqIt: make_range(scalarRange.first,scalEqIt))
         {
            // loop over the implicit range
            for(auto& fIt: make_range(doneSEqIt->couplingInfo(compId).implicitRange()))
            {
               // Check if field is in implicit range
               if(fIt == fieldId)
               {
                  // Set the solver index
                  (*scalEqIt)->setSolverIndex(compId, doneSEqIt->couplingInfo(compId).solverIndex());

                  // Debug statements
                  DebuggerMacro_msg("Identified coupled scalar solver: " + PhysicalNames::Coordinator::formatted((*scalEqIt)->name()), 2);
                  DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());

                  // Set coupling flag and break out
                  coupled = true;
                  break;
               }
            }

            // Break out of loop if already identified
            if(coupled)
            {
               break;
            }
         }

         // All checked and equation is not coupled
         if(!coupled)
         {
            // Set solver index for complex equation
            if((*scalEqIt)->couplingInfo(compId).isComplex())
            {
               // Set complex solver index
               (*scalEqIt)->setSolverIndex(compId, zIdx);

               // Increment complex solver index
               zIdx++;

            // Set solver index for real equation
            } else
            {
               // Set real solver index
               (*scalEqIt)->setSolverIndex(compId, dIdx);

               // Increment real solver index
               dIdx++;
            }

            // Debug statements
            DebuggerMacro_msg("Identified first scalar solver: " + PhysicalNames::Coordinator::formatted((*scalEqIt)->name()), 2);
            DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
            DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());
         }

         // Reset coupling flag
         coupled = false;
      }

      // Loop over the vector equations
      for(auto vectEqIt = vectorRange.first; vectEqIt != vectorRange.second; ++vectEqIt)
      {
         // Get coupled counter for each component
         ArrayI counter((*vectEqIt)->nSpectral());
         counter.setConstant(0);

         // Loop over the (identified) scalar equations
         for(auto& doneSEqIt: make_range(scalarRange))
         {
            // loop over the implicit range
            for(auto& fIt: make_range(doneSEqIt->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange()))
            {
               int i = 0;
               for(auto& compIt: make_range((*vectEqIt)->spectralRange()))
               {
                  // Check if field's first component is in implicit range
                  if(counter(i) == 0 && fIt == std::make_pair((*vectEqIt)->name(), compIt))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compIt, doneSEqIt->couplingInfo(FieldComponents::Spectral::SCALAR).solverIndex());

                     // Debug statements
                     DebuggerMacro_msg("Identified coupled vector solver: " + PhysicalNames::Coordinator::formatted((*vectEqIt)->name()), 2);
                     DebuggerMacro_showValue("---> component: ", 2, compIt);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compIt).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compIt).isComplex());

                     // Set coupling flag and break out
                     counter(i) = 1;
                  }
                  // Increment coupled counter index
                  i++;
               }

               // Break out of loop if already identified
               if(counter.sum() == counter.size())
               {
                  break;
               }
            }

            // Break out of loop if already identified
            if(counter.sum() == counter.size())
            {
               break;
            }
         }

         if(counter.sum() < counter.size())
         {
            // Loop over the (identified) vector equations
            for(auto& doneVEqIt: make_range(vectorRange.first, vectEqIt))
            {
               for(auto& doneIt: make_range(doneVEqIt->spectralRange()))
               {
                  // loop over the implicit range of identified components
                  for(auto& fIt: make_range(doneVEqIt->couplingInfo(doneIt).implicitRange()))
                  {
                     int i = 0;
                     for(auto& compIt: make_range((*vectEqIt)->spectralRange()))
                     {
                        // Check if field's first component is in implicit range
                        if(counter(i) == 0 && fIt == std::make_pair((*vectEqIt)->name(), compIt))
                        {
                           // Set the solver index
                           (*vectEqIt)->setSolverIndex(compIt, doneVEqIt->couplingInfo(doneIt).solverIndex());

                           // Debug statements
                           DebuggerMacro_msg("Identified coupled vector solver: " + PhysicalNames::Coordinator::formatted((*vectEqIt)->name()), 2);
                           DebuggerMacro_showValue("---> component: ", 2, compIt);
                           DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compIt).solverIndex());
                           DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compIt).isComplex());

                           // Set coupling flag and break out
                           counter(i) = 1;
                        }
                        // Increment coupled counter index
                        i++;
                     }


                     // Break out of loop if already identified
                     if(counter.sum() == counter.size())
                     {
                        break;
                     }
                  }

                  // Break out of loop if already identified
                  if(counter.sum() == counter.size())
                  {
                     break;
                  }
               }
            }
         }

         // All checked and equation is not coupled
         if(counter.sum() < counter.size())
         {
            std::vector<FieldComponents::Spectral::Id> selfCoupling;
            FieldComponents::Spectral::Id coupledComp;
            int i = 0;
            for(auto& compIt: make_range((*vectEqIt)->spectralRange()))
            {
               if(counter(i) == 0)
               {
                  bool selfCoupled = false;
                  for(auto& fIt: make_range((*vectEqIt)->couplingInfo(compIt).implicitRange()))
                  {
                     if(fIt.first == (*vectEqIt)->name())
                     {
                        selfCoupled = (std::find(selfCoupling.begin(), selfCoupling.end(), fIt.second) != selfCoupling.end());
                     }

                     if(selfCoupled)
                     {
                        coupledComp = fIt.second;
                        break;
                     }
                  }

                  if(selfCoupled)
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(compIt, (*vectEqIt)->couplingInfo(coupledComp).solverIndex());

                     // Debug statements
                     DebuggerMacro_msg("Identified coupled vector solver: " + PhysicalNames::Coordinator::formatted((*vectEqIt)->name()), 2);
                     DebuggerMacro_showValue("---> component: ", 2, compIt);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compIt).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compIt).isComplex());

                     // Set coupling flag and break out
                     counter(i) = 1;
                  } else
                  {
                     // Set solver index for complex equation
                     if((*vectEqIt)->couplingInfo(compIt).isComplex())
                     {
                        // Set complex solver index
                        (*vectEqIt)->setSolverIndex(compIt, zIdx);

                        // Increment complex solver index
                        zIdx++;

                        // Set solver index for real equation
                     } else
                     {
                        // Set real solver index
                        (*vectEqIt)->setSolverIndex(compIt, dIdx);

                        // Increment real solver index
                        dIdx++;
                     }

                     // Debug statements
                     DebuggerMacro_msg("Identified first vector solver: " + PhysicalNames::Coordinator::formatted((*vectEqIt)->name()), 2);
                     DebuggerMacro_showValue("---> component: ", 2, compIt);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(compIt).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(compIt).isComplex());
                  }
               }
               selfCoupling.push_back(compIt);
               // Increment coupled counter index
               i++;
            }
         }

         // Reset coupling flag
         counter.setConstant(0);
      }
   }

   template <typename TIt> void getNonlinearKernels(std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& rKernels, TIt first, TIt last)
   {
      for(TIt eqIt = first; eqIt != last; ++eqIt)
      {
         (*eqIt)->initNLKernel();
         rKernels.insert(std::make_pair((*eqIt)->name(), (*eqIt)->spNLKernel()));
      }
   }
}
}
}

#endif // QUICC_EQUATIONS_TOOLS_EQUATIONTOOLS_HPP
