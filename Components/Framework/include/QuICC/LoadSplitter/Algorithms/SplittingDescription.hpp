/**
 * @file SplittingDescription.hpp
 * @brief Implementation basic description of the obtained splitting (for any algorithm)
 */

#ifndef QUICC_PARALLEL_SPLITTINGDESCRIPTION_HPP
#define QUICC_PARALLEL_SPLITTINGDESCRIPTION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Io/Xml/VtpWriter.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation basic description of the obtained splitting (for any algorithm)
    */
   class SplittingDescription
   {
      public:
         /// Flag to enable or disable storing distribution
         static bool store;

         /**
          * @brief Constructor
          */
         SplittingDescription() = default;

         /**
          * @brief Destructor
          */
         ~SplittingDescription() = default;

         /**
          * @brief ID of the algorithm
          */
         Splitting::Algorithms::Id algorithm;

         /**
          * @brief ID of the grouper
          */
         Splitting::Groupers::Id grouper;

         /**
          * @brief Number of dimensions
          */
         int dims;

         /**
          * @brief Storage for the \f$N_{cpu}\f$ factorisation factors
          */
         ArrayI factors;

         /**
          * @brief Score of the splitting
          */
         Array score;

         /**
          * @brief Communication structure
          */
         std::map<Dimensions::Transform::Id,std::multimap<int,int> >   structure;

         /**
          * @brief Add stage to description
          */
         void addStage(const int stageId, const SharedTransformResolution spTRes, const int cpuId);

         /**
          * @brief Save stage descriptions
          */
         void save();

      protected:

      private:
         /**
          * @brief Name stages based on ID
          */
         std::string nameStage(const int stageId) const;

         /**
          * @brief Storage for data distribution visualization files
          */
         std::map<int,Io::Xml::SharedVtpWriter> vtpFiles;
   };

}
}

#endif // QUICC_PARALLEL_SPLITTINGDESCRIPTION_HPP
