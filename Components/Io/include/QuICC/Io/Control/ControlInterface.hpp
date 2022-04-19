/**
 * @file ControlInterface.hpp
 * @brief Implementation of a external runtime control file
 */

#ifndef QUICC_IO_CONTROL_CONTROLINTERACE_HPP
#define QUICC_IO_CONTROL_CONTROLINTERACE_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "QuICC/Io/Control/ControlFile.hpp"

namespace QuICC {

namespace Io {

namespace Control {

   /**
    * @brief Implementation of an external runtime control file
    */
   class ControlInterface: public ControlFile
   {
      public:
         /**
          * @brief Constructor
          */
         ControlInterface(const bool allowedIo = true);

         /**
          * @brief Destructor
          */
         ~ControlInterface();

         /**
          * @brief Read input or check existance
          */
         void read();

         /**
          * @brief Current runtime status?
          */
         std::size_t status() const;

      protected:

      private:
         /**
          * @brief Is IO allowed?
          */
         bool mAllowedIo;

         /**
          * @brief Handle to the file
          */
         std::fstream mFile;

         /**
          * @brief Current runtime status
          */
         std::size_t mStatus;

         /**
          * @brief Initialise control file
          */
         void init();

         /**
          * @brief Finalise the file
          */
         void finalize();

         /**
          * @brief Create the file
          */
         void create();

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Delete the file
          */
         void deleteFile();

         /**
          * @brief Update status
          */
         void update();
   };
}
}
}

#endif // QUICC_IO_CONTROL_CONTROLINTERACE_HPP
