/**
 * @file TimestepCoordinator.hpp
 * @brief Implementation of the exponential timestep coordinator
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP

// System includes
//
#include <memory>
#include <string>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "Timestep/Exponential/Tags.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/AugmentedJacobianFunctor.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#include "Types/Typedefs.hpp"
#include "View/ViewDense.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

   template <template <class,class,typename> class TStepper, typename TImpl> class TimestepperCoordinator;

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class,typename> class TStepper> class TimestepperCoordinator<TStepper, base_t>
   {
      public:
         /// Typedef for exponential timestepper
         typedef TStepper<SparseMatrix,Matrix,base_t>  TimestepperType;

         /// Typedef for a shared exponential timestepper
         typedef typename std::shared_ptr<TimestepperType >  SharedTimestepperType;

         /// Typedef for mat of shared exponential timesteppers
         typedef std::map<std::size_t, std::tuple<SharedTimestepperType, std::vector<std::size_t>> > TimestepperMap;

         /**
          * @brief Constructor
          */
         TimestepperCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~TimestepperCoordinator() = default;

         /**
          * @brief Initialise coordinator
          *
          * @param sizeInfo  vector of size information
          */
         template <typename TScheme>
         void init(const MHDFloat dt, const std::vector<TimestepperInfo>& info, std::shared_ptr<TScheme> spScheme, std::shared_ptr<AugmentedJacobianFunctor> pJac);

         /**
          * @brief Reset solvers for next step
          */
         void resetSolvers();

         /**
          * @brief Get current time
          */
         MHDFloat stepFraction();

         /**
          * @brief Step forward
          */
         void stepForward();

         /**
          * @brief Update timestep
          */
         void updateTimestep(const MHDFloat dt);

         /**
          * @brief Get timestepper and start row
          */
         std::pair<TimestepperType*,std::size_t> getStepper(const TimestepperInfo& info);

      protected:
         /**
          * @brief add a timestepper
          */
         void addTimestepper(const TimestepperInfo& info, TimestepperMap& steppers);

         /**
          * @brief Update start row
          *
          * @param sizeInfo  Size information
          */
         void updateStartRow(const TimestepperInfo& info, TimestepperMap& steppers);

         /**
          * @brief init a component of timestepper
          */
         void initComponent(const TimestepperInfo& info, TimestepperMap& steppers);

         /**
          * @brief init a timestepper
          */
         void initTimestepper(const TimestepperInfo& info, TimestepperMap& steppers, std::shared_ptr<AugmentedJacobianFunctor> pJac);

         /**
          * @brief Vector of (coupled) exponential timestepper
          */
         TimestepperMap mSteppers;

      private:
         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

         /**
          * @brief Shared Jacobian functor
          */
         std::shared_ptr<AugmentedJacobianFunctor> mpJac;

   };

   template <template <class,class,typename> class TStepper> TimestepperCoordinator<TStepper, base_t>::TimestepperCoordinator()
   {
   }

   template <template <class, class, typename> class TStepper>
   void TimestepperCoordinator<TStepper, base_t>::updateTimestep(const MHDFloat dt)
   {
      this->mDt = dt;

      this->mpJac->updateTimestep(dt);
   }

   template <template <class,class,typename> class TStepper> template <typename TScheme> void TimestepperCoordinator<TStepper,base_t>::init(const MHDFloat dt, const std::vector<TimestepperInfo>& infos, std::shared_ptr<TScheme> spScheme, std::shared_ptr<AugmentedJacobianFunctor> pJac)
   {
      this->mDt = dt;

      this->mpJac = pJac;
      this->mpJac->updateTimestep(dt);

      //
      // Create timesteppers
      //

      DebuggerMacro_msg("Creating " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
            this->addTimestepper(info, this->mSteppers);
      }
      DebuggerMacro_msg("... done", 2);

      //
      // Update the start rows
      //

      DebuggerMacro_msg("Updating start row for " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
         this->updateStartRow(info, this->mSteppers);
      }
      DebuggerMacro_msg("... done", 2);

      //
      // Set timestepping scheme
      for(auto& tsData: this->mSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.setScheme(spScheme);
      }

      //
      // Init timestepper components
      //

      DebuggerMacro_msg("Initializing " + std::to_string(infos.size()) + " timestepper components", 2);
      for(auto&& info: infos)
      {
         this->initComponent(info, this->mSteppers);
      }
      DebuggerMacro_msg("... done", 2);

      //
      // Init timesteppers
      //

      DebuggerMacro_msg("Initializing " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
         this->initTimestepper(info, this->mSteppers, this->mpJac);
      }
      DebuggerMacro_msg("... done", 2);
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::addTimestepper(const TimestepperInfo& info, TimestepperMap& steppers)
   {
      if(steppers.count(info.solverIndex) == 0)
      {
         auto spSolver = std::make_shared<TimestepperType>();

         std::vector<std::size_t> startRow(info.fieldIndex + 1, 0);
         auto stepData = std::make_tuple(spSolver, startRow);
         steppers.emplace(info.solverIndex, stepData);
      }
      else
      {
         // Add start indexes
         auto& stepData = steppers.at(info.solverIndex);
         auto& startRow = std::get<1>(stepData);
         for(std::size_t i = startRow.size(); i < info.fieldIndex + 1; i++)
         {
            startRow.push_back(0);
         }
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateStartRow(const TimestepperInfo& info, TimestepperMap& steppers)
   {
      assert(steppers.count(info.solverIndex) > 0);

      // Update start indexes
      auto& stepData = steppers.at(info.solverIndex);
      auto& startRow = std::get<1>(stepData);
      for(std::size_t i = info.fieldIndex + 1; i <  startRow.size(); i++)
      {
         startRow.at(i) += info.blockN;
      }
   }

   template <template <class,class,typename> class TStepper> std::pair<typename TimestepperCoordinator<TStepper, base_t>::TimestepperType*,std::size_t> TimestepperCoordinator<TStepper, base_t>::getStepper(const TimestepperInfo& info)
   {
      assert(this->mSteppers.count(info.solverIndex) > 0);
      auto& stepData = this->mSteppers.at(info.solverIndex);
      auto spStepper = std::get<0>(stepData);
      auto& startArr = std::get<1>(stepData);
      assert(startArr.size() > info.fieldIndex);
      std::size_t start = startArr.at(info.fieldIndex) + info.matStart;

      auto ts = std::make_pair(spStepper.get(), start);

      return ts;
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::initComponent(const TimestepperInfo& info, TimestepperMap& steppers)
   {
      DebuggerMacro_msg("Processing field ID =  " + std::to_string(info.fieldIndex) + " for solver index " + std::to_string(info.solverIndex), 3);

      auto& stepData = steppers.at(info.solverIndex);
      auto spStepper = std::get<0>(stepData);
      auto&& startArr = std::get<1>(stepData);
      spStepper->addStorage(info.rows, info.cols);
      spStepper->initMatrices(info.matIds, startArr.at(info.fieldIndex));
      spStepper->buildOperators(info.ops, this->mDt, startArr.at(info.fieldIndex));
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::initTimestepper(const TimestepperInfo& info, TimestepperMap& steppers, std::shared_ptr<AugmentedJacobianFunctor> pJac)
   {
      auto& stepData = steppers.at(info.solverIndex);
      auto spStepper = std::get<0>(stepData);
      if(!spStepper->isInitialized())
      {
         DebuggerMacro_msg("Processing field ID =  " + std::to_string(info.fieldIndex) + " for solver index " + std::to_string(info.solverIndex), 3);
         spStepper->initPhi(pJac);
         spStepper->setInitialized();
         spStepper->initSolver();
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::resetSolvers()
   {
      for(auto& tsData: this->mSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.resetSolver();
      }
   }

   template <template <class,class,typename> class TStepper> MHDFloat TimestepperCoordinator<TStepper,  base_t>::stepFraction()
   {
      return 1.0;
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper,base_t>::stepForward()
   {
      // Run timesteppers
      DebuggerMacro_msg("Running timesteppers", 6);
      for(auto& [key,stepData]: this->mSteppers)
      {
         auto spStepper = std::get<0>(stepData);
         auto&& startArr = std::get<1>(stepData);
         auto tsWrapper = std::make_shared<AugmentedJacobianFunctor::TsFunctor>(startArr);
         tsWrapper->setStepper(*spStepper);
         this->mpJac->setStepper(tsWrapper);
         spStepper->stepForward();
      }
      DebuggerMacro_msg("... done", 6);
   }

} // Exponential
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP
