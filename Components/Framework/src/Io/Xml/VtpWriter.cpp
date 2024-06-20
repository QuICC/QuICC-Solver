/**
 * @file VtpWriter.cpp
 * @brief Source of the implementation of the VTK PolyData XML format file writer
 */

// Configuration includes
//

// System includes
//
#include <set>

// External includes
//
#include <rapidxml_print.hpp>

// Class include
//
#include "QuICC/Io/Xml/VtpWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/DimensionTools.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   VtpWriter::VtpWriter(const std::string& name)
      : IVtpFile<Io::Xml::IXmlWriter>(name)
   {
   }

   VtpWriter::~VtpWriter()
   {
   }

   void VtpWriter::representResolution(SharedTransformResolution spTRes, const int rank)
   {
      std::stringstream oss;

      // Node variables
      rapidxml::xml_node<>* pPiece;
      rapidxml::xml_node<>* pPoints;
      rapidxml::xml_node<>* pData;

      // Get master PolyData tag
      pPiece = this->mXML.first_node(this->HEADER.c_str());
      rapidxml::xml_node<> *pPoly = pPiece->first_node(this->TYPE.c_str());

      // Piece for 1D distribution
      pPiece = this->mXML.allocate_node(rapidxml::node_element, "Piece");
      pPiece->append_attribute(this->mXML.allocate_attribute("NumberOfVerts", "0"));
      pPiece->append_attribute(this->mXML.allocate_attribute("NumberOfLines", "0"));
      pPiece->append_attribute(this->mXML.allocate_attribute("NumberOfStrips", "0"));
      pPiece->append_attribute(this->mXML.allocate_attribute("NumberOfPolys", "0"));
      pPoly->append_node(pPiece);

      // Points node
      pPoints = this->mXML.allocate_node(rapidxml::node_element, "Points");
      pPiece->append_node(pPoints);

      // DataArray node
      pData = this->mXML.allocate_node(rapidxml::node_element, "DataArray");
      pData->append_attribute(this->mXML.allocate_attribute("type", "Int32"));
      pData->append_attribute(this->mXML.allocate_attribute("NumberOfComponents", "3"));
      pData->append_attribute(this->mXML.allocate_attribute("format", "ascii"));
      int count = 0;
      for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            int j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
            oss << 0 << " " << j_ << " " << k_ << " ";
            count++;
         }
      }
      pData->value(this->mXML.allocate_string(oss.str().c_str()));
      oss.str("");
      pPoints->append_node(pData);

      // Write points count to file
      oss << count;
      pPiece->append_attribute(this->mXML.allocate_attribute("NumberOfPoints", this->mXML.allocate_string(oss.str().c_str())));
      oss.str("");

      // PointData node
      pPoints = this->mXML.allocate_node(rapidxml::node_element, "PointData");
      pPoints->append_attribute(this->mXML.allocate_attribute("Scalars", "rank"));
      pPiece->append_node(pPoints);

      // DataArray node for rank
      pData = this->mXML.allocate_node(rapidxml::node_element, "DataArray");
      pData->append_attribute(this->mXML.allocate_attribute("type", "Int32"));
      pData->append_attribute(this->mXML.allocate_attribute("Name", "rank"));
      pData->append_attribute(this->mXML.allocate_attribute("format", "ascii"));
      for(int n = 0; n < count; n++)
      {
         oss << rank << " ";
      }
      pData->value(this->mXML.allocate_string(oss.str().c_str()));
      oss.str("");
      pPoints->append_node(pData);
   }

   void VtpWriter::write()
   {
      // Do pre writing processing
      this->preWrite();

      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Write xml content to file
         this->mFile << this->mXML;
      }

      // Do post writing processing
      this->postWrite();
   }
}
}
}
