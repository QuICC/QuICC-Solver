/**
 * @file RuntimeStatistics.hpp
 * @brief Runtime statistics
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_RUNTIMESTATISTICS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_RUNTIMESTATISTICS_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief Struct for holding runtime statistics
    */
   struct RuntimeStatistics
   {
      /**
       * @brief ctor
       */
      RuntimeStatistics();

      /**
       * @brief ctor
       */
      ~RuntimeStatistics() = default;

      /**
       * @brief Update global statistics
       */
      void update();

      /**
       * @brief Reset statistics
       */
      void reset();

      /**
       * @brief Print statistics
       */
      void printInfo();

      int step;
      int krystep;
      int reject;
      int exps;
      int m;
      int mMax;
      int mMin;
      int gram_p;
      double conv;
      double convMax;
      double convMin;
   };

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_RUNTIMESTATISTICS_HPP
