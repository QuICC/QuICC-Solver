/**
 * @file SparseCoordinatorData.hpp
 * @brief Implementation of the base for a general sparse solver coordinator
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_SPARSECOORDINATORDATA_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_SPARSECOORDINATORDATA_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#include "Types/Typedefs.hpp"
#include "Timestep/PredictorCorrector/Views/Tags.hpp"
#include "View/ViewDense.hpp"

#include <iostream>
namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

   struct TimestepperInfo
   {
      bool isComplex;
      std::size_t solverIndex;
      std::size_t fieldIndex;
      std::size_t rows;
      std::size_t cols;
      std::size_t blockN;
      std::map<std::size_t, DecoupledZSparse> ops;
   };

   template <template <class,class,typename> class TStepper, typename TImpl> class TimestepperCoordinator;

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class,typename> class TStepper> class TimestepperCoordinator<TStepper, base_t>
   {
      public:
         /// Typedef for a real operator timestepper
         typedef TStepper<SparseMatrix,DecoupledZMatrix,base_t>  RealTimestepperType;

         /// Typedef for a complex operator timestepper
         typedef TStepper<SparseMatrixZ,MatrixZ,base_t>  ComplexTimestepperType;

         /// Typedef for a shared real operator timestepper
         typedef typename std::shared_ptr<RealTimestepperType >  SharedRealTimestepperType;

         /// Typedef for a shared complex operator timestepper
         typedef typename std::shared_ptr<ComplexTimestepperType >  SharedComplexTimestepperType;

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
         void init(const MHDFloat dt, const std::vector<TimestepperInfo>& info, std::shared_ptr<TScheme> spScheme);

         /**
          * @brief Clear the RHS data of all solvers
          */
         void clearSolvers();

         /**
          * @brief Get current time
          */
         MHDFloat stepFraction();

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

         /**
          * @brief Get error measure
          */
         MHDFloat error() const;

         /**
          * @brief Check if step is done
          */
         bool finishedStep() const;

         /**
          * @brief Update timestep
          */
         void updateTimestep(const MHDFloat dt);

         /**
          * @brief Update time dependencies
          */
         void updateMatrices();

         /**
          * @brief Update RHS
          */
         void updateRhs(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs);

         /**
          * @brief Update RHS
          */
         void updateSolution(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol);

         /**
          * @brief Get solution
          */
         void getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const TimestepperInfo& info);

      protected:
         /**
          * @brief add a timestepper
          */
         template <typename T> void addTimestepper(const TimestepperInfo& info, T& steppers);

         /**
          * @brief Update start row
          *
          * @param sizeInfo  Size information
          */
         template <typename T> void updateStartRow(const TimestepperInfo& info, T& steppers);

         /**
          * @brief init a timestepper
          */
         template <typename T> void initTimestepper(const TimestepperInfo& info, T& steppers);

         /**
          * @brief Run timesteppers
          */
         template <typename T> std::pair<bool,MHDFloat> runSteppers(T& steppers);

         /**
          * @brief Vector of (coupled) real operator
          */
         std::map<std::size_t, std::tuple<SharedRealTimestepperType, std::vector<std::size_t>> > mRealSteppers;

         /**
          * @brief Vector of (coupled) complex operator
          */
         std::map<std::size_t, std::tuple<SharedComplexTimestepperType, std::vector<std::size_t>> > mComplexSteppers;

         /**
          * @brief Flag to signal end of computation
          */
         bool mFinished;

         /**
          * @brief Computation error
          */
         MHDFloat mError;

      private:
         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

   };

   template <template <class,class,typename> class TStepper> TimestepperCoordinator<TStepper, base_t>::TimestepperCoordinator()
      : mFinished(false), mError(-1.0)
   {
   }

   template <template <class,class,typename> class TStepper> MHDFloat TimestepperCoordinator<TStepper, base_t>::error() const
   {
      return this->mError;
   }

   template <template <class,class,typename> class TStepper> bool TimestepperCoordinator<TStepper, base_t>::finishedStep() const
   {
      return this->mFinished;
   }

   template <template <class, class, typename> class TStepper>
   void TimestepperCoordinator<TStepper, base_t>::updateTimestep(const MHDFloat dt)
   {
      this->mDt = dt;
   }

   template <template <class,class,typename> class TStepper> template <typename TScheme> void TimestepperCoordinator<TStepper,base_t>::init(const MHDFloat dt, const std::vector<TimestepperInfo>& infos, std::shared_ptr<TScheme> spScheme)
   {
      this->mDt = dt;

      //
      // Create timesteppers
      //

      DebuggerMacro_msg("Creating " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
         if(info.isComplex)
         {
            this->addTimestepper(info, this->mComplexSteppers);
         }
         else
         {
            this->addTimestepper(info, this->mRealSteppers);
         }
      }
      DebuggerMacro_msg("... done", 2);

      //
      // Update the start rows
      //

      DebuggerMacro_msg("Updating start row for " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
         if(info.isComplex)
         {
            this->updateStartRow(info, this->mComplexSteppers);
         }
         else
         {
            this->updateStartRow(info, this->mRealSteppers);
         }
      }
      DebuggerMacro_msg("... done", 2);

      //
      // Set timestepping scheme
      for(auto& tsData: this->mComplexSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.setScheme(spScheme);
      }
      for(auto& tsData: this->mRealSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.setScheme(spScheme);
      }

      //
      // Init timesteppers
      //

      DebuggerMacro_msg("Initializing " + std::to_string(infos.size()) + " timesteppers", 2);
      for(auto&& info: infos)
      {
         if(info.isComplex)
         {
            this->initTimestepper(info, this->mComplexSteppers);
         }
         else
         {
            this->initTimestepper(info, this->mRealSteppers);
         }
      }
      DebuggerMacro_msg("... done", 2);
   }

   template <template <class,class,typename> class TStepper> template <typename T> void TimestepperCoordinator<TStepper, base_t>::addTimestepper(const TimestepperInfo& info, T& steppers)
   {
      if(steppers.count(info.solverIndex) == 0)
      {
         auto spSolver = std::make_shared<typename std::tuple_element_t<0,typename T::mapped_type>::element_type>();

         std::vector<std::size_t> startRow(info.fieldIndex + 1, 0);
         auto stepData = std::make_tuple(spSolver, startRow);
         steppers.emplace(info.solverIndex, stepData);
      }
      else
      {
         // Add start indexes
         auto& stepData = steppers.at(info.solverIndex);
         auto& startRow = std::get<1>(stepData);
         for(std::size_t i = 0; i < info.fieldIndex + 1 -startRow.size(); i++)
         {
            startRow.push_back(0);
         }
      }
   }

   template <template <class,class,typename> class TStepper> template <typename T> void TimestepperCoordinator<TStepper, base_t>::updateStartRow(const TimestepperInfo& info, T& steppers)
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

   template <template <class,class,typename> class TStepper> template <typename T> void TimestepperCoordinator<TStepper, base_t>::initTimestepper(const TimestepperInfo& info, T& steppers)
   {
      auto& stepData = steppers.at(info.solverIndex);
      auto spStepper = std::get<0>(stepData);
      if(!spStepper->isInitialized())
      {
         spStepper->addStorage(info.rows, info.cols);
         spStepper->initMatrices();
         spStepper->buildOperators(info.ops, this->mDt, info.rows);
         spStepper->setInitialized();
         spStepper->initSolver();
         spStepper->zeroSolver();
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateSolution(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol)
   {
      if(info.isComplex)
      {
         assert(this->mComplexSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mComplexSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->setSolution(sol, startArr.at(info.fieldIndex));
         spStepper->updateSolutions();
      }
      else
      {
         assert(this->mRealSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mRealSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->setSolution(sol, startArr.at(info.fieldIndex));
         spStepper->updateSolutions();
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateRhs(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs)
   {
      if(info.isComplex)
      {
         assert(this->mComplexSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mComplexSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->addRhs(rhs, startArr.at(info.fieldIndex));
      }
      else
      {
         assert(this->mRealSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mRealSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->addRhs(rhs, startArr.at(info.fieldIndex));
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const TimestepperInfo& info)
   {
      if(info.isComplex)
      {
         assert(this->mComplexSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mComplexSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->getSolution(sol, startArr.at(info.fieldIndex));
      }
      else
      {
         assert(this->mRealSteppers.count(info.solverIndex) > 0);
         auto& stepData = this->mRealSteppers.at(info.solverIndex);
         auto spStepper = std::get<0>(stepData);
         auto& startArr = std::get<1>(stepData);
         assert(startArr.size() > info.fieldIndex);
         spStepper->getSolution(sol, startArr.at(info.fieldIndex));
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::clearSolvers()
   {
      for(auto& tsData: this->mRealSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.zeroSolver();
      }

      for(auto& tsData: this->mComplexSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.zeroSolver();
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateMatrices()
   {
      for(auto& tsData: this->mRealSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.updateTimeMatrix(this->mDt);
         ts.updateSolver();
      }

      for(auto& tsData: this->mComplexSteppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         ts.updateTimeMatrix(this->mDt);
         ts.updateSolver();
      }
   }

   template <template <class,class,typename> class TStepper> MHDFloat TimestepperCoordinator<TStepper,  base_t>::stepFraction()
   {
      if(this->mComplexSteppers.size() > 0)
      {
         auto& s = *std::get<0>(this->mComplexSteppers.begin()->second);
         return s.stepFraction();
      } else
      {
         auto& s = *std::get<0>(this->mRealSteppers.begin()->second);
         return s.stepFraction();
      }
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper,base_t>::solveSystems()
   {
      // Run complex timesteppers
      DebuggerMacro_msg("Running complex timesteppers", 6);
      std::pair<bool,MHDFloat> zStatus = this->runSteppers(this->mComplexSteppers);
      DebuggerMacro_msg("... done", 6);

      // Run real timesteppers
      DebuggerMacro_msg("Running real timesteppers", 6);
      std::pair<bool,MHDFloat> dStatus = this->runSteppers(this->mRealSteppers);
      DebuggerMacro_msg("... done", 6);

      this->mFinished = zStatus.first || dStatus.first;

      if(this->mFinished)
      {
         this->mError = std::max(zStatus.second, dStatus.second);
      }
   }

   template <template <class,class,typename> class TStepper> template <typename T> std::pair<bool,MHDFloat> TimestepperCoordinator<TStepper, base_t>::runSteppers(T& steppers)
   {
      std::pair<bool,MHDFloat>  status = std::make_pair(false, -1.0);
      for(auto& tsData: steppers)
      {
         auto& ts = *std::get<0>(tsData.second);
         bool solving = false;
         do
         {
            // Prepare solve of linear system
            bool needSolve = ts.preSolve();

            if(needSolve)
            {
               // Solve linear system
               ts.solve();

               // Work on fields after solve
               solving = ts.postSolve();

            } else
            {
               solving = false;
            }

         } while (solving);

         status.first = ts.finished();

         if(status.first)
         {
            status.second = std::max(status.second, ts.error());
         }
      }

      return status;
   }

} // Views
} // PredictorCorrector
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_SPARSECOORDINATORDATA_HPP
