/** 
 * @file TransformCommSelector.hpp
 * @brief Typedefs to setup the correct transforms communicator and coordinator
 */

#ifndef QUICC_PARALLEL_TRANSFORMCOMMSELECTOR_HPP
#define QUICC_PARALLEL_TRANSFORMCOMMSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/TransformCoordinators/TransformCoordinator.hpp"
#include "QuICC/Communicators/Communicator.hpp"

namespace QuICC {

   namespace Transform {
      typedef TransformCoordinator<Parallel::Communicator> TransformCoordinatorType;
   }
}

#endif // QUICC_PARALLEL_TRANSFORMCOMMSELECTOR_HPP
