/**
 * @file GxlWriter.cpp
 * @brief Source of the implementation of the GXL format file writer
 */

// Configuration includes
//

// System includes
//
#include <set>
#include <stdexcept>

// External includes
//
#include <rapidxml_print.hpp>

// Class include
//
#include "QuICC/Io/Xml/GxlWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Arithmetics/Set.hpp"
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   GxlWriter::GxlWriter(const std::string& name)
      : IGxlFile<Io::Xml::IXmlWriter>(name)
   {
   }

   GxlWriter::~GxlWriter()
   {
   }

   void GxlWriter::graphCommunication(const std::vector<std::multimap<int,int> >& structure)
   {
      std::stringstream oss;

      // Set subgraph long identifier
      std::vector<std::string> lName;
      lName.push_back("cluster_A");
      lName.push_back("cluster_B");

      // Set subgraph short identifier
      std::vector<std::string> sName;
      sName.push_back("A");
      sName.push_back("B");

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the two transposes
      for(size_t i = 0; i < structure.size(); i++)
      {
         // Create subgraph
         rapidxml::xml_node<> *pNSubgraph = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "N_" << lName.at(i);
         pNSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
         oss.str("");
         rapidxml::xml_node<> *pSubgraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
         pSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(lName.at(i).c_str(),0)));

         if(static_cast<Dimensions::Transform::Id>(i) == Dimensions::Transform::TRA1D)
         {
            this->createAttr(pSubgraph, "label", "Transpose 1D/2D");
         } else
         {
            this->createAttr(pSubgraph, "label", "Transpose 2D/3D");
         }
         this->createAttr(pSubgraph, "color", color.at(i));

         pNSubgraph->append_node(pSubgraph);
         pGraph->append_node(pNSubgraph);

         // Loop over the CPUs
         for(auto itCpu = structure.at(i).cbegin(); itCpu != structure.at(i).cend(); itCpu = structure.at(i).upper_bound(itCpu->first))
         {
            rapidxml::xml_node<> *pCpu = this->mXML.allocate_node(rapidxml::node_element, "node");
            oss << "cpu" << itCpu->first << sName.at(i);
            pCpu->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            oss << "CPU " << itCpu->first;
            this->createAttr(pCpu, "label", oss.str());
            oss.str("");

            pSubgraph->append_node(pCpu);
         }

         // Loop over the CPUs
         for(auto itCpu = structure.at(i).cbegin(); itCpu != structure.at(i).cend(); ++itCpu)
         {
            // Set "from" attribute
            rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
            oss << "cpu" << itCpu->first << sName.at(i);
            pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            // Set "to" attribute
            oss << "cpu" << itCpu->second << sName.at(i);
            pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            // Add edge
            this->createAttr(pEdge, "color", color.at(i));
            pSubgraph->append_node(pEdge);
         }
      }
   }

   void GxlWriter::graphTransformPath(const std::map<std::size_t, std::vector<Transform::TransformPath> >& paths, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");
      color.push_back("red");

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the paths
      int unid = 0;
      for(auto nameIt = paths.cbegin(); nameIt != paths.cend(); ++nameIt)
      {
         for(auto pathIt = nameIt->second.cbegin(); pathIt != nameIt->second.cend(); ++pathIt)
         {
            std::string field;

            rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
            oss << "p" << unid << "_" << nameIt->first << "_" << pathIt->startId();
            field = oss.str();
            oss.str("");
            pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(field.c_str(),0)));

            oss << PhysicalNames::Coordinator::formatted(nameIt->first) << " ";
            if(dir == TransformDirection::FORWARD)
            {
               oss << Tools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(pathIt->startId()));
            } else
            {
               oss << Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(pathIt->startId()));
            }
            this->createAttr(pNode, "label", oss.str());
            oss.str("");

            pGraph->append_node(pNode);

            // Add all path nodes
            std::string curNode = field;
            for(int i = 0; i < pathIt->nEdges(); ++i)
            {
               rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
               oss << field << "_" << i;
               pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));

               // Set "from" attribute
               rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
               pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(curNode.c_str(),0)));

               // Set "to" attribute
               pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(oss.str().c_str(),0)));
               curNode = oss.str();
               oss.str("");

               // Add edge
               this->createAttr(pEdge, "color", color.at(i));
               oss << pathIt->edge(i).opId();
               this->createAttr(pEdge, "label", oss.str());
               oss.str("");
               pGraph->append_node(pEdge);

               if(pathIt->edge(i).arithId() == Arithmetics::Set::id())
               {
                  oss << "==";
               } else if(pathIt->edge(i).arithId() == Arithmetics::SetNeg::id())
               {
                  oss << "=-";
               } else if(pathIt->edge(i).arithId() == Arithmetics::Add::id())
               {
                  oss << "++";
               } else if(pathIt->edge(i).arithId() == Arithmetics::Sub::id())
               {
                  oss << "--";
               }
               if(pathIt->edge(i).outId().at(0) != -1)
               {
                  if(pathIt->fieldId() == FieldType::SCALAR)
                  {
                     oss << std::endl << "Scalar";
                  } else if(pathIt->fieldId() == FieldType::CURL)
                  {
                     oss << std::endl << "Curl";
                  } else if(pathIt->fieldId() == FieldType::GRADIENT)
                  {
                     oss << std::endl << "Gradient";
                  } else if(pathIt->fieldId() == FieldType::GRADIENT2)
                  {
                     oss << std::endl << "Gradient2";
                  } else
                  {
                     oss << std::endl;
                  }

                  for(auto outIt = pathIt->edge(i).outId().cbegin(); outIt != pathIt->edge(i).outId().cend(); ++outIt)
                  {
                     if(dir == TransformDirection::FORWARD)
                     {
                        oss << " " << Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(*outIt));
                     } else
                     {
                        oss << " " << Tools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(*outIt));
                     }
                  }
               }
               this->createAttr(pNode, "label", oss.str());
               oss.str("");

               pGraph->append_node(pNode);
            }

            unid++;
         }
      }
   }

   void GxlWriter::graphTransformTree(const std::vector<Transform::TransformTree>& trees, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");
      color.push_back("red");
      std::vector<std::string>::const_iterator colorIt = color.begin();

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the trees
      for(auto treeIt = trees.cbegin(); treeIt != trees.cend(); ++treeIt)
      {
         rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "t" << "_" << treeIt->name() << "_" << treeIt->comp<int>();
         std::string root = oss.str();
         oss.str("");
         pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(root.c_str(),0)));

         oss << PhysicalNames::Coordinator::formatted(treeIt->name()) << " ";
         if(dir == TransformDirection::FORWARD)
         {
            oss << Tools::IdToHuman::toString(treeIt->comp<FieldComponents::Physical::Id>());
         } else
         {
            oss << Tools::IdToHuman::toString(treeIt->comp<FieldComponents::Spectral::Id>());
         }
         this->createAttr(pNode, "label", oss.str());
         oss.str("");

         pGraph->append_node(pNode);

         graphTransformTreeEdge(treeIt->root(), root, colorIt, pGraph, dir);
      }
   }

   void GxlWriter::graphTransformTreeEdge(const Transform::TransformTreeEdge& edge, const std::string& root, std::vector<std::string>::const_iterator colorIt, rapidxml::xml_node<> * pGraph, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      for(auto edgeIt = edge.edgeRange().first; edgeIt != edge.edgeRange().second; ++edgeIt)
      {
         rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << root << edgeIt->opId();
         std::string nextRoot = oss.str();
         oss.str("");
         pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(nextRoot.c_str(),0)));

         // Set "from" attribute
         rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
         pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(root.c_str(),0)));

         // Set "to" attribute
         pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(nextRoot.c_str(),0)));

         // Add edge
         this->createAttr(pEdge, "color", *colorIt);
         oss << edgeIt->opId();
         if(edgeIt->recoverInput() && edgeIt->holdInput())
         {
            oss << std::endl << "(R,H)";
         } else if(edgeIt->recoverInput())
         {
            oss << std::endl << "(R)";
         } else if(edgeIt->holdInput())
         {
            oss << std::endl << "(H)";
         }
         this->createAttr(pEdge, "label", oss.str());
         oss.str("");
         pGraph->append_node(pEdge);

         if(edgeIt->arithId() == Arithmetics::Set::id())
         {
            oss << "==";
         } else if(edgeIt->arithId() == Arithmetics::SetNeg::id())
         {
            oss << "=-";
         } else if(edgeIt->arithId() == Arithmetics::Add::id())
         {
            oss << "++";
         } else if(edgeIt->arithId() == Arithmetics::Sub::id())
         {
            oss << "--";
         }
         if(edgeIt->combinedOutId() >= 0 || edgeIt->recoverOutId() >= 0)
         {
            oss << "[";
            if(edgeIt->recoverOutId() >= 0)
            {
               oss << edgeIt->recoverOutId() << ", ";
            }
            if(edgeIt->combinedArithId() == Arithmetics::Set::id())
            {
               oss << "==";
            } else if(edgeIt->combinedArithId() == Arithmetics::SetNeg::id())
            {
               oss << "=-";
            } else if(edgeIt->combinedArithId() == Arithmetics::Add::id())
            {
               oss << "++";
            } else if(edgeIt->combinedArithId() == Arithmetics::Sub::id())
            {
               oss << "--";
            }
            if(edgeIt->combinedOutId() >= 0)
            {
               oss << ", ->" << edgeIt->combinedOutId();
            }
            oss << "]";
         }
         if(edgeIt->outId<int>() != -1)
         {
            oss << std::endl;
            if(edgeIt->fieldId() == FieldType::SCALAR)
            {
               oss << "Scalar";
            } else if(edgeIt->fieldId() == FieldType::VECTOR)
            {
               // Just use name
            } else if(edgeIt->fieldId() == FieldType::CURL)
            {
               oss << "Curl";
            } else if(edgeIt->fieldId() == FieldType::GRADIENT)
            {
               oss << "Gradient";
            } else if(edgeIt->fieldId() == FieldType::GRADIENT2)
            {
               oss << "Gradient2";
            } else
            {
               throw std::logic_error("Unknown field type ID in tree");
            }

            for(auto outIt = edgeIt->outIds().cbegin(); outIt != edgeIt->outIds().cend(); ++outIt)
            {
               if(dir == TransformDirection::FORWARD)
               {
                  oss << " " << Tools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(*outIt));
               } else
               {
                  oss << " " << Tools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(*outIt));
               }
            }
         }
         this->createAttr(pNode, "label", oss.str());
         oss.str("");

         pGraph->append_node(pNode);

         colorIt++;
         graphTransformTreeEdge(*edgeIt, nextRoot, colorIt, pGraph, dir);
         colorIt--;
      }
   }

   void GxlWriter::createAttr(rapidxml::xml_node<>* sParent, const std::string& name, const std::string& value)
   {
      rapidxml::xml_node<>* pAttr;
      rapidxml::xml_node<>* pStr;

      pAttr = this->mXML.allocate_node(rapidxml::node_element, "attr");
      pAttr->append_attribute(this->mXML.allocate_attribute("name", this->mXML.allocate_string(name.c_str()),0));
      pStr = this->mXML.allocate_node(rapidxml::node_element, "string");
      pStr->value(this->mXML.allocate_string(value.c_str(),0));

      pAttr->append_node(pStr);
      sParent->append_node(pAttr);
   }

   void GxlWriter::write()
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
