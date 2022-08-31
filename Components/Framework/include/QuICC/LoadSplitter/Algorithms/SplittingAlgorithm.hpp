/** 
 * @file SplittingAlgorithm.hpp
 * @brief Base of the implementation of the load splitting algorithms
 */

#ifndef QUICC_PARALLEL_SPLITTINGALGORITHM_HPP
#define QUICC_PARALLEL_SPLITTINGALGORITHM_HPP

// Configuration includes
//

// System includes
//
#include <list>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Base of the implementation of the load splitting algorithms
    */
   class SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU
          * @param nCpu    Number of CPUs used
          * @param dim     Dimensions
          * @param algo    ID of the algorithm
          */
         SplittingAlgorithm(const int id, const int nCpu, const ArrayI& dim, const Splitting::Algorithms::Id algo);

         /**
          * @brief Destructor
          */
         virtual ~SplittingAlgorithm();

         /**
          * @brief Compute splitting and the corresponding score
          */
         std::pair<int, std::pair<SharedResolution, SplittingDescription> > scoreSplitting(const Splitting::Groupers::Id grp);

         /**
          * @brief Set the spatial scheme builder
          *
          * @param spBuilder  Spatial scheme builder
          */
         void setScheme(SpatialScheme::SharedIBuilder spBuilder);

         /**
          * @brief Select next set of splitting factors
          */
         bool useNextFactors();

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         virtual bool applicable() const = 0;

         /**
          * @brief Build the exact communication structure
          *
          * @param spRes   Shared resolution object
          */
         static void buildCommunicationStructure(const int localId, SharedResolution spRes, std::vector<std::multimap<int,int> >& commStructure);
         
      protected:
         /**
          * @brief ID of the algorithm
          */
         Splitting::Algorithms::Id mAlgo;

         /**
          * @brief ID of the grouper
          */
         Splitting::Groupers::Id mGrouper;

         /**
          * @brief Shared spatial scheme
          */
         SpatialScheme::SharedIBuilder mspScheme;

         /**
          * @brief Initialise the CPU factors
          *
          * @param nFactors   Number of factors
          */
         void initFactors(const int nFactors);

         /**
          * @brief Split dimension i
          *
          * @param transId Split the ith dimension
          * @param cpuId   ID of the CPU
          * @param status  Status output
          */
         virtual SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status) = 0;

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual Array computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp) = 0;

         /**
          * @brief Compute score related to communication structure
          *
          * @param spRes   Shared resolution object
          * @param details Details of the communication structure
          */
         double communicationScore(SharedResolution spRes, ArrayI& details);

         /**
          * @brief Compute score related to load balancing
          *
          * @param spRes   Shared resolution object
          * @param balance Details of the load balancing (on input they contain weights)
          */
         double balancingScore(SharedResolution spRes, Array& balance);

         /**
          * @brief Get id of core
          */
         int id() const;

         /**
          * @brief Get Number of CPUs
          */
         int nCpu() const;

         /**
          * @brief Dimensionality of the fields
          */
         int dims() const;

         /**
          * @brief Get \f$F_{i}\f$ factor in \f$N_{cpu} = \prod_{i} F_{i}\f$
          *
          * @param i index of the factor
          */
         int factor(const int i) const;

         /**
          * @brief Get \f$F_{i}\f$ factors in \f$N_{cpu} = \prod_{i} F_{i}\f$
          */
         const ArrayI& factors() const;

         /**
          * @brief Get maximum \f$F_{i}\f$ factor in \f$N_{cpu} = \prod_{i} F_{i}\f$
          */
         int maxFactor() const;

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
          * @brief Number of dimensions
          */
         int mDims;

         /**
          * @brief Storage for the \f$N_{cpu}\f$ factorisation factors
          */
         ArrayI mFactors;

         /**
          * @brief List of \f$N_{cpu}\f$ factorisation factors
          */
         std::list<int> mNCpuFactors;

         /**
          * @brief Storage for the simulation resolution
          */
         ArrayI mSimDim;
   };

   /// Typedef for a shared pointer to a SplittingAlgorithm object
   typedef std::shared_ptr<SplittingAlgorithm>   SharedSplittingAlgorithm;

}
}

#endif // QUICC_PARALLEL_SPLITTINGALGORITHM_HPP
