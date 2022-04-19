/**
 * @file SimulationIoTools.hpp
 * @brief Implementation of the tools IO related calculations
 */

#ifndef QUICC_SIMULATIONIOTOOLS_HPP
#define QUICC_SIMULATIONIOTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "QuICC/Io/Stats/IStatisticsAsciiWriter.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the tools for IO related calculations
    */
   class SimulationIoTools
   {
      public:
         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<Io::Variable::SharedIVariableAsciiWriter>::iterator ascii_iterator;

         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<Io::Stats::SharedIStatisticsAsciiWriter>::iterator stats_iterator;

         /**
          * @brief Update heavy calculations for ASCII file output
          */
         static void updateHeavyAscii(ascii_iterator asciiBegin, ascii_iterator asciiEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update pre calculations for statistics file output
          */
         static void updateStatsPre(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update calculations for statistics file output
          */
         static void updateStats(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update post calculations for statistics file output
          */
         static void updateStatsPost(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SimulationIoTools();

         /**
          * @brief Empty destructor
          */
         ~SimulationIoTools();
   };

}

#endif // QUICC_SIMULATIONIOTOOLS_HPP
