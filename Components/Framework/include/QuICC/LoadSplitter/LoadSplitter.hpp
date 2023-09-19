/**
 * @file LoadSplitter.hpp
 * @brief Implementation of the workload splitter over the available CPUs
 */

#ifndef QUICC_PARALLEL_LOADSPLITTER_HPP
#define QUICC_PARALLEL_LOADSPLITTER_HPP

// System includes
//
#include <vector>
#include <set>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the workload splitter over the available CPUs
    */
   class LoadSplitter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU
          * @param nCpu    Number of CPUs used
          */
         LoadSplitter(const int id, const int nCpu);

         /**
          * @brief Destructor
          */
         ~LoadSplitter() = default;

         /**
          * @brief Initialise the algorithms
          *
          * @param spBuilder  Spatial scheme builder
          * @param enabled    List of enabled algorithms
          * @param grp        Communication grouper
          * @param factors    Imposed factorization of CPU
          */
         void init(SpatialScheme::SharedIBuilder spBuilder, const std::set<Splitting::Algorithms::Id>& enabled, const Splitting::Groupers::Id grp, const std::list<int>& factors = {});

         /**
          * @brief Get splitting information of the best splitting
          */
         std::pair<SharedResolution,SplittingDescription> bestSplitting(const bool fillCommStructure = true);

         /**
          * @brief Show description of some splittings
          *
          * @param n Maximum number of splittings to show
          */
         void showSplittings(const int n) const;

      protected:

      private:
         /**
          * @brief ID of CPU
          */
         int mId;

         /**
          * @brief Number of CPUs
          */
         int mNCpu;

         /**
          * @brief Maximum number of stored scores
          */
         int mMaxNScores;

         /**
          * @brief Storage for the actual splitting algorithms used
          */
         std::vector<SharedSplittingAlgorithm>  mAlgorithms;

         /**
          * @brief Storage for the scores of the splitting
          */
         std::multimap<int, std::pair<SharedResolution,SplittingDescription> >  mScores;

         /**
          * @brief Initialise the splitting algorithms
          *
          * @param dim     Dimensions (spectral)
          * @param enabled List of enabled algorithms
          * @param factors Imposed factorization of nCPU
          */
         void initAlgorithms(const ArrayI& dim, const std::set<Splitting::Algorithms::Id>& enabled, const std::list<int>& factors);

         /**
          * @brief Initialise the scores and corresponding splitting
          */
         void initScores(const Splitting::Groupers::Id grp);

         /**
          * @brief Describe the obtained splitting
          *
          * @param descr Splitting description object
          */
         void describeSplitting(const SplittingDescription& descr, const bool isTest = false) const;
   };
}
}

#endif // QUICC_PARALLEL_LOADSPLITTER_HPP
