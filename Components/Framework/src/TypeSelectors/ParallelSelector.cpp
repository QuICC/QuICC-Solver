/**
 * @file ParallelSelector.cpp
 * @brief Definitions of the ParallelSelector grouper selection
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/TypeSelectors/ParallelSelector.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   void setGrouper(const SplittingDescription& descr, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
   {
      if(descr.grouper == Splitting::Groupers::EQUATION)
      {
         setGrouper<Splitting::Groupers::EQUATION>(descr.algorithm, descr.dims, spFwdGrouper, spBwdGrouper);

   #ifdef QUICC_MPI
      } else if(descr.grouper == Splitting::Groupers::SINGLE1D)
      {
         setGrouper<Splitting::Groupers::SINGLE1D>(descr.algorithm, descr.dims, spFwdGrouper, spBwdGrouper);
      } else if(descr.grouper == Splitting::Groupers::SINGLE2D)
      {
         setGrouper<Splitting::Groupers::SINGLE2D>(descr.algorithm, descr.dims, spFwdGrouper, spBwdGrouper);
      } else if(descr.grouper == Splitting::Groupers::TRANSFORM)
      {
         setGrouper<Splitting::Groupers::TRANSFORM>(descr.algorithm, descr.dims, spFwdGrouper, spBwdGrouper);
   #endif //QUICC_MPI
      } else
      {
         throw std::logic_error("Unknown transform grouper setup");
      }
   }
}
}
