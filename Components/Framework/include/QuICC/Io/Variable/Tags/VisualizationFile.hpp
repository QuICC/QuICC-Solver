/** 
 * @file VisualizationFile.hpp
 * @brief Definitions and names use by the visualisation file writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_TAGS_VISUALIZATIONFILE_HPP
#define QUICC_IO_VARIABLE_TAGS_TAGS_VISUALIZATIONFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Variable {

namespace Tags {

   /**
    * @brief Definitions and names use by the visualisation file writers
    */
   class VisualizationFile
   {
      public:
         /**
          * @brief HEADER part for visualisation file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for visualisation file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of visualisation file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of visualisation file
          */
         static const std::string   EXTENSION;

         /**
          * @brief Mesh description for visualisation file
          */
         static const std::string   MESH;

         /**
          * @brief Grid description for visualisation file
          */
         static const std::string   GRID;

      private:
         /**
         * @brief Empty destructor
         */
         VisualizationFile();

         /**
         * @brief Destructor
         */
         ~VisualizationFile();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_TAGS_VISUALIZATIONFILE_HPP
