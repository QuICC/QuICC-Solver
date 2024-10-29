/**
 * @file TestModel.hpp
 * @brief Test model setup to validate IO
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_IO_TESTMODEL_HPP
#define QUICC_TESTSUITE_FRAMEWORK_IO_TESTMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Generator/VisualizationGenerator.hpp"
#include "QuICC/Model/IPhysicalPyModel.hpp"
#include "QuICC/Simulation/Simulation.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection sphere
 * model (Toroidal/Poloidal formulation)
 */
class TestModel : public Model::IPhysicalPyModel<Simulation, StateGenerator,
                     VisualizationGenerator>
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::WLFl SchemeType;

   /**
    * @brief Constructor
    */
   TestModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~TestModel() = default;

   /// Python script/module name
   std::string PYMODULE() final;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

   /// Formulation used for vector fields
   virtual VectorFormulation::Id SchemeFormulation() override;

   /**
    * @brief Version string
    */
   std::string version() const final;

   /**
    * @brief Add equation
    *
    * @param spSim   Shared simulation object
    */
   virtual void addEquations(SharedSimulation spSim) override {};

   /**
    * @brief Add the initial state generation equations
    *
    * @param spGen   Shared generator object
    */
   virtual void addStates(SharedStateGenerator spGen) override;

   /**
    * @brief Add visualization
    *
    * @param spGen   Shared generator object
    */
   virtual void addVisualizers(SharedVisualizationGenerator spViz) override {};

   /**
    * @brief Add the required ASCII output files
    *
    * @param spSim   Shared simulation object
    */
   virtual void addAsciiOutputFiles(SharedSimulation spSim) override;

   /**
    * @brief Add the required ASCII output files
    *
    * @param spSim   Shared simulation object
    */
   virtual void addAsciiOutputFiles(SharedStateGenerator spSim);

   /**
    * @brief XML configuration tags
    */
   virtual std::map<std::string, std::map<std::string, int>>
   configTags() const override;

protected:
private:
};

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_IO_TESTMODEL_HPP
