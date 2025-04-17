/**
 * @file Affinity.hpp
 * @brief Affinity utils for thread pool
 */

#ifndef QUICC_THREADPOOL_DETAILS_AFFINITY_HPP
#define QUICC_THREADPOOL_DETAILS_AFFINITY_HPP

// System includes
//
#include <boost/thread/barrier.hpp>

// Project includes
//

namespace QuICC {

namespace ThreadPool {

namespace details {

   /**
    * @brief Print CPU affinity mask
    */
   int printAffinity();

   /**
    * @brief Count CPU affinity mask
    */
   int countAffinity();

   /**
    * @brief Set CPU affinity mask
    *
    * @param bar Barrier
    * @param idx CPU index
    */
   int setAffinity(boost::barrier& bar, const int idx);

} // namespace details
} // namespace ThreadPool
} // namespace QuICC

#endif // QUICC_THREADPOOL_DETAILS_AFFINITY_HPP
