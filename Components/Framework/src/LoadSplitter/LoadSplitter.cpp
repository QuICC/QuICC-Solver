/** 
 * @file LoadSplitter.cpp
 * @brief Source of the workload splitter
 */

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <sstream>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/LoadSplitter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Tools/Formatter.hpp"

// Splitting algorithms
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/LoadSplitter/Algorithms/SerialSplitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SingleSplitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/TubularSplitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/Coupled2DSplitting.hpp"
#include "QuICC/Io/Xml/GxlWriter.hpp"
#include "QuICC/Io/Xml/VtpWriter.hpp"

namespace QuICC {

namespace Parallel {

   LoadSplitter::LoadSplitter(const int id, const int nCpu)
      : mId(id), mNCpu(nCpu), mMaxNScores(1)
   {
   }

   LoadSplitter::~LoadSplitter()
   {
   }

   void LoadSplitter::initAlgorithms(const ArrayI& dim, const std::set<Splitting::Algorithms::Id>& enabled)
   {
      // Check for serial version of code request
      if(this->mNCpu == 1)
      {
         this->mAlgorithms.push_back(std::make_shared<SerialSplitting>(this->mId, this->mNCpu, dim));

      // Setup the parallel version (initialise all algorithms and then choose the best one)
      } else
      {
         // Check that problem is at least 2D
         if(dim.size() > 1)
         {
            // Add the single splitting algorithm for first data exchange
            if(enabled.count(Splitting::Algorithms::SINGLE1D) == 1)
            {
               this->mAlgorithms.push_back(std::make_shared<SingleSplitting>(this->mId, this->mNCpu, dim, Splitting::Locations::FIRST));
            }

            // Add the single splitting algorithm for first data exchange for coupled dimensions
            if(enabled.count(Splitting::Algorithms::COUPLED2D) == 1)
            {
               this->mAlgorithms.push_back(std::make_shared<Coupled2DSplitting>(this->mId, this->mNCpu, dim));
            }


            // Check if problem is 3D
            if(dim.size() == 3)
            {
               // Add the single splitting algorithm for second data exchange
               if(enabled.count(Splitting::Algorithms::SINGLE2D) == 1)
               {
                  this->mAlgorithms.push_back(std::make_shared<SingleSplitting>(this->mId, this->mNCpu, dim, Splitting::Locations::SECOND));
               }

               // Add the tubular splitting algorithm
               if(enabled.count(Splitting::Algorithms::TUBULAR) == 1)
               {
                  this->mAlgorithms.push_back(std::make_shared<TubularSplitting>(this->mId, this->mNCpu, dim));
               }
            }

         // 1D problems can't be parallelised with the current algorithms
         } else
         {
            throw std::logic_error("There is no parallelisation algorithm for 1D problems!");
         }
      }

      // Safety check to make sure at least one algorithm has been initialised
      if(this->mAlgorithms.size() == 0)
      {
         throw std::logic_error("No algorithm has been initialised!");
      }
   }

   void LoadSplitter::init(SpatialScheme::SharedIBuilder spBuilder, const std::set<Splitting::Algorithms::Id>& enabled, const Splitting::Groupers::Id grp)
   {
      // Initialise the splitting algorithms
      this->initAlgorithms(spBuilder->resolution(), enabled);

      // Loop over all initialised algorithms
      for(auto it = this->mAlgorithms.begin(); it != this->mAlgorithms.end(); it++)
      {
         (*it)->setScheme(spBuilder);
      }

      // Compute the core resolutions and corresponding scores
      this->initScores(grp);
   }

   void LoadSplitter::initScores(const Splitting::Groupers::Id grp)
   {
      // Loop over all initialised algorithms
      for(auto it = this->mAlgorithms.begin(); it != this->mAlgorithms.end(); it++)
      {
         // Loop over possible factorisations of nCPU
         while((*it)->useNextFactors())
         {
            // Check if obtained factorisation is applicable to splitting algorithm
            if((*it)->applicable())
            {
               // Get scored resolution object
               std::pair<int, std::pair<SharedResolution, SplittingDescription> > scored = (*it)->scoreSplitting(grp);

               // Only scores bigger than 0 are considered
               if(scored.first > 0)
               {
                  // If not at minimal size yet, add to scores
                  if(this->mScores.size() < static_cast<unsigned int>(this->mMaxNScores))
                  {
                     this->mScores.insert(scored);
                  } else
                  {
                     if(scored.first > this->mScores.begin()->first)
                     {
                        // Remove lowest score
                        this->mScores.erase(this->mScores.begin());
                        // Add new score
                        this->mScores.insert(scored);
                     }
                  }

                  #ifdef QUICC_DEBUG
                  // Show description of tested splittings
                  this->describeSplitting(scored.second.second, true);
                  #endif // QUICC_DEBUG
               }
            }
         }
      }
   }

   std::pair<SharedResolution,SplittingDescription> LoadSplitter::bestSplitting(const bool needCommStructure)
   {
      // Make sure there is at least one successful splitting (score > 0)
      if(this->mScores.size() > 0)
      {
         if(needCommStructure)
         {
            // Build communication structure
            SplittingAlgorithm::buildCommunicationStructure(this->mId, this->mScores.rbegin()->second.first, this->mScores.rbegin()->second.second.structure);

            // Describe the splitting with the highest score
            this->describeSplitting(this->mScores.rbegin()->second.second);

            #ifdef QUICC_DEBUG
            for(auto it = this->mScores.rbegin()->second.second.vtpFiles.cbegin(); it != this->mScores.rbegin()->second.second.vtpFiles.cend(); ++it)
            {
               (*it)->write();
               (*it)->finalize();
            }
            #endif //QUICC_DEBUG
         }

         // Return the splitting with the highest score
         return this->mScores.rbegin()->second;
      } else
      {
         throw std::logic_error("No usable splitting has been found!");
      }
   }

   void LoadSplitter::describeSplitting(const SplittingDescription& descr, const bool isTest) const
   {
      // Output a short description of the selected splitting. Make it look nice ;)
      if(QuICCEnv().allowsIO())
      {
         if(isTest)
         {
            // Print load splitting test header
            Tools::Formatter::printLine(std::cout, '~');
            Tools::Formatter::printCentered(std::cout, "Tested load distribution", '~');
            Tools::Formatter::printLine(std::cout, '~');
         } else
         {
            // Print load splitting header
            Tools::Formatter::printLine(std::cout, '-');
            Tools::Formatter::printCentered(std::cout, "Load distribution", '*');
            Tools::Formatter::printLine(std::cout, '-');
         }

         std::string tmpStr;

         // Print grouper information
         switch(descr.grouper)
         {
            case(Splitting::Groupers::EQUATION):
               tmpStr = "Equation";
               break;
            case(Splitting::Groupers::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splitting::Groupers::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splitting::Groupers::TRANSFORM):
               tmpStr = "Transform";
               break;
         }
         Tools::Formatter::printCentered(std::cout, "Grouper: " + tmpStr);

         // Print Algorithm information
         switch(descr.algorithm)
         {
            case(Splitting::Algorithms::SERIAL):
               tmpStr = "Serial";
               break;
            case(Splitting::Algorithms::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splitting::Algorithms::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splitting::Algorithms::TUBULAR):
               tmpStr = "Tubular";
               break;
            case(Splitting::Algorithms::COUPLED2D):
               tmpStr = "Coupled 2D";
               break;
         }
         Tools::Formatter::printCentered(std::cout, "Algorithm: " + tmpStr);

         // Print factorisation information
         tmpStr = "";
         std::stringstream oss;
         for(int i = 0; i < descr.factors.size();i++)
         {
            oss << descr.factors(i);
            tmpStr += oss.str();
            oss.str("");
            if(i < descr.factors.size() - 1)
            {
               tmpStr += " x ";
            }
         }
         if(descr.algorithm == Splitting::Algorithms::SINGLE1D || descr.algorithm == Splitting::Algorithms::COUPLED2D)
         {
            tmpStr += " x 1";
         } else if(descr.algorithm == Splitting::Algorithms::SINGLE2D)
         {
            tmpStr = "1 x " + tmpStr;
         }
         Tools::Formatter::printCentered(std::cout, "Factorisation: " + tmpStr);

         // Create GXL graph format file
         tmpStr.erase(std::remove(tmpStr.begin(),tmpStr.end(),' '),tmpStr.end());
         Io::Xml::GxlWriter gxl("Communication_graph_" + tmpStr);
         gxl.init();
         gxl.graphCommunication(descr.structure);
         gxl.write();
         gxl.finalize();

         oss.str("");
         oss << static_cast<int>(descr.score.prod());
         Tools::Formatter::printCentered(std::cout, "Score: " + oss.str());
         oss.str("");
         oss << " (" << descr.score(0) << ", " << std::setprecision(2)<< descr.score(1) << ", " << descr.score(2) << ", " << descr.score(3) << ")";
         Tools::Formatter::printCentered(std::cout, oss.str());

         Tools::Formatter::printNewline(std::cout);
      }

      // Synchronize
      QuICCEnv().synchronize();
   }

   void LoadSplitter::showSplittings(const int n) const
   {
      // Get maximum between number of scores and n
      int maxN = std::min(static_cast<int>(this->mScores.size()), n);

      // Create reverse iterator
      std::multimap<int, std::pair<SharedResolution,SplittingDescription> >::const_reverse_iterator rit;

      // Set start iterator
      rit = this->mScores.rbegin();

      // Loop over scores
      for(int i = 0; i < maxN; ++i)
      {
         // Describe the obtained splitting
         this->describeSplitting(rit->second.second);

         // Increment iterator
         rit++;
      }
   }

}
}
