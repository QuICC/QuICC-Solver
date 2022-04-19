/**
 * @file VtpWriter.hpp
 * @brief Implementation of the VTK PolyData XML format file writer
 */

#ifndef QUICC_IO_XML_VTPWRITER_HPP
#define QUICC_IO_XML_VTPWRITER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Xml/IXmlWriter.hpp"
#include "QuICC/Io/Xml/IVtpFile.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   /**
    * @brief Implementation of the VTK PolyData XML format file writer
    */
   class VtpWriter: public IVtpFile<Io::Xml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         VtpWriter(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~VtpWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();

         /**
          * @brief Convert resolution to XML for VTP file
          */
         void representResolution(SharedTransformResolution spRes, const int rank);

      protected:

      private:
   };

   /// Typedef for a smart pointer of a VtpWriter
   typedef std::shared_ptr<VtpWriter> SharedVtpWriter;

}
}
}

#endif // QUICC_IO_XML_VTPWRITER_HPP
