/**
 * @file CflWriter.hpp
 * @brief Implementation of a CFL writer
 */

#ifndef QUICC_IO_ASCII_CFLWRITER_HPP
#define QUICC_IO_ASCII_CFLWRITER_HPP

// System includes
//
#include <set>
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Io/Ascii/IAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief Implementation of a CFL writer
    */
   class CflWriter: public Io::Ascii::IAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          */
         CflWriter();

         /**
          * @brief Destructor
          */
         virtual ~CflWriter();

         /**
          * @brief Set the simulation time parameters
          *
          * @param time       Reached simulation time
          * @param cfl        CFL timesteps
          * @param steps      Number of steps with previous timestep
          */
         void setSimTime(const MHDFloat time, const Matrix& cfl, const MHDFloat steps);

      protected:
         /**
          * @brief Write State to file
          */
         virtual void writeContent();

         /**
          * @brief Time
          */
         MHDFloat mTime;

         /**
          * @brief Timestep
          */
         MHDFloat mTimestep;

         /**
          * @brief Number of steps with same timestep
          */
         MHDFloat mSteps;

         /**
          * @brief Error
          */
         MHDFloat mError;

         /**
          * @brief Flag for timestep change
          */
         bool mChanged;

         /**
          * @brief Timestep details
          */
         Matrix mDt;

      private:
         /**
          * @brief Write fancy header
          */
         virtual void fancyHeader();

         /**
          * @brief High precision output
          */
         const int mcIoHigh;

         /**
          * @brief Low precision output
          */
         const int mcIoLow;

         /**
          * @brief Width of integer
          */
         const int mcIoIntW;

         /**
          * @brief Need fancy header?
          */
         bool mNeedFancy;
   };

   /// Typedef for a smart reference counting pointer of a CflWriter writer
   typedef std::shared_ptr<CflWriter>   SharedCflWriter;

}
}
}

#endif // QUICC_IO_ASCII_CFLWRITER_HPP
