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
#include "QuICC/Typedefs.hpp"
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
         /**
          * @brief Constructor
          *
          * @param algorithm  ID of the algorithm
          * @param grouper    ID of the transform grouper
          * @param dims       Number of dimensions
          * @param factors    CPU splitting factors
          * @param score      Score 
          */
         SplittingDescription();

         /**
          * @brief Destructor
          */
         ~SplittingDescription();

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
         std::vector<std::multimap<int,int> >   structure;

         #ifdef QUICC_DEBUG
         /**
          * @brief Storage for data distribution visualization files
          */
         std::vector<Io::Xml::SharedVtpWriter> vtpFiles;
         #endif //QUICC_DEBUG
         
      protected:

      private:
   };

}
}

#endif // QUICC_PARALLEL_SPLITTINGDESCRIPTION_HPP
