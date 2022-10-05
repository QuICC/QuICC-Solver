#!/usr/bin/env python

def header(descr):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief Reductor transform operator {name}
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_{cName.upper()}_HPP
#define QUICC_TRANSFORM_REDUCTOR_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Reductor/IRegisterId.hpp"

namespace QuICC {{

namespace Transform {{

namespace Reductor {{

   /**
    * @brief Reductor transform operator {name}
    */
   class {cName}: public IRegisterId<{cName}>
   {{
      public:
         /**
          * @brief Constructor
          */
         {cName}();

         friend class IRegisterId<{cName}>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   }};

}} // Reductor
}} // Transform
}} // QuICC

#endif // QUICC_TRANSFORM_REDUCTOR_{cName.upper()}_HPP
'''
    return content

def source(descr):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the reductor transform operator {name}
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace Transform {{

namespace Reductor {{

   std::string {cName}::sTag()
   {{
      return "{tag}";
   }}

   std::string {cName}::sFormatted()
   {{
      return "{name}";
   }}

   {cName}::{cName}()
      : IRegisterId<{cName}>({cName}::sTag(), {cName}::sFormatted())
   {{
   }}

}} // Reductor
}} // Transform
}} // QuICC
'''
    return content

names = []
names.append({'class':"Energy", 'name':"Reductor::Energy", 'tag':"Red::Energy"})
names.append({'class':"EnergyD1", 'name':"Reductor::EnergyD1", 'tag':"Red::EnergyD1"})
names.append({'class':"EnergyD1R1", 'name':"Reductor::EnergyD1R1", 'tag':"Red::EnergyD1R1"})
names.append({'class':"EnergyR2", 'name':"Reductor::EnergyR2", 'tag':"Red::EnergyR2"})
names.append({'class':"EnergySLAPLR2", 'name':"Reductor::EnergySLAPLR2", 'tag':"Red::EnergySLAPLR2"})
names.append({'class':"Power", 'name':"Reductor::Power", 'tag':"Red::Power"})
names.append({'class':"PowerD1", 'name':"Reductor::PowerD1", 'tag':"Red::PowerD1"})
names.append({'class':"PowerD1R1", 'name':"Reductor::PowerD1R1", 'tag':"Red::PowerD1R1"})
names.append({'class':"PowerR2", 'name':"Reductor::PowerR2", 'tag':"Red::PowerR2"})
names.append({'class':"PowerSLAPLR2", 'name':"Reductor::PowerSLAPLR2", 'tag':"Red::PowerSLAPLR2"})
names.append({'class':"RadialPower", 'name':"Reductor::RadialPower", 'tag':"Red::RadialPower"})
names.append({'class':"RadialPowerDivR1", 'name':"Reductor::RadialPowerDivR1", 'tag':"Red::RadialPowerDivR1"})
names.append({'class':"RadialPowerDivR1D1R1", 'name':"Reductor::RadialPowerDivR1D1R1", 'tag':"Red::RadialPowerDivR1D1R1"})

incPath = '../include/QuICC/Transform/Reductor/'
srcPath = '../src/Transform/Reductor/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n))
    f.close()

