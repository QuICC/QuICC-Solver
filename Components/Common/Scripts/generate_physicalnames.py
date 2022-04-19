#!/usr/bin/env python

def header(descr):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief {name} physical name 
 */

#ifndef QUICC_PHYSICALNAMES_{cName.upper()}_HPP
#define QUICC_PHYSICALNAMES_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalNames/IRegisterId.hpp"

namespace QuICC {{

namespace PhysicalNames {{

   /**
    * @brief {name} physical name
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

}}
}}

#endif // QUICC_PHYSICALNAMES_{cName.upper()}_HPP
'''
    return content

def source(descr):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the {name} physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace PhysicalNames {{

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

}}
}}
'''
    return content

names = []
names.append({'class':"Codensity", 'name':"Codensity", 'tag':"codensity"})
names.append({'class':"Density", 'name':"Density", 'tag':"density"})
names.append({'class':"DxMeanTemperature", 'name':"D_x mean temperature", 'tag':"dx_meantemperature"})
names.append({'class':"DzMeanTemperature", 'name':"D_z mean temperature", 'tag':"dz_meantemperature"})
names.append({'class':"Entropy", 'name':"Entropy", 'tag':"entropy"})
names.append({'class':"FluctMagnetic", 'name':"Fluctuating magnetic", 'tag':"fluct_magnetic"})
names.append({'class':"FluctMagneticX", 'name':"Fluctuating magnetic X", 'tag':"fluct_magneticx"})
names.append({'class':"FluctMagneticY", 'name':"Fluctuating magnetic Y", 'tag':"fluct_magneticy"})
names.append({'class':"FluctMagneticZ", 'name':"Fluctuating magnetic Z", 'tag':"fluct_magneticz"})
names.append({'class':"FluctTemperature", 'name':"Fluctuating temperature", 'tag':"fluct_temperature"})
names.append({'class':"FluctVelocity", 'name':"Fluctuating velocity", 'tag':"fluct_velocity"})
names.append({'class':"FluctVelocityX", 'name':"Fluctuating velocity X", 'tag':"fluct_velocityx"})
names.append({'class':"FluctVelocityY", 'name':"Fluctuating velocity Y", 'tag':"fluct_velocityy"})
names.append({'class':"FluctVelocityZ", 'name':"Fluctuating velocity Z", 'tag':"fluct_velocityz"})
names.append({'class':"Magnetic", 'name':"Magnetic", 'tag':"magnetic"})
names.append({'class':"MagneticX", 'name':"Magnetic X", 'tag':"magneticx"})
names.append({'class':"MagneticY", 'name':"Magnetic Y", 'tag':"magneticy"})
names.append({'class':"MagneticZ", 'name':"Magnetic Z", 'tag':"magneticz"})
names.append({'class':"MeanMagnetic", 'name':"Mean magnetic", 'tag':"mean_magnetic"})
names.append({'class':"MeanMagneticX", 'name':"Mean magnetic X", 'tag':"mean_magneticx"})
names.append({'class':"MeanMagneticY", 'name':"Mean magnetic Y", 'tag':"mean_magneticy"})
names.append({'class':"MeanMagneticZ", 'name':"Mean magnetic Z", 'tag':"mean_magneticz"})
names.append({'class':"MeanTemperature", 'name':"Mean temperature", 'tag':"mean_temperature"})
names.append({'class':"MeanVelocity", 'name':"Mean velocity", 'tag':"mean_velocity"})
names.append({'class':"MeanVelocityX", 'name':"Mean velocity X", 'tag':"mean_velocityx"})
names.append({'class':"MeanVelocityY", 'name':"Mean velocity Y", 'tag':"mean_velocityy"})
names.append({'class':"MeanVelocityZ", 'name':"Mean velocity Z", 'tag':"mean_velocityz"})
names.append({'class':"Pressure", 'name':"Pressure", 'tag':"pressure"})
names.append({'class':"Streamfunction", 'name':"Streamfunction", 'tag':"streamfunction"})
names.append({'class':"Temperature", 'name':"Temperature", 'tag':"temperature"})
names.append({'class':"Velocity", 'name':"Velocity", 'tag':"velocity"})
names.append({'class':"VelocityX", 'name':"Velocity X", 'tag':"velocityx"})
names.append({'class':"VelocityY", 'name':"Velocity Y", 'tag':"velocityy"})
names.append({'class':"VelocityZ", 'name':"Velocity Z", 'tag':"velocityz"})
names.append({'class':"Vorticity", 'name':"Vorticity", 'tag':"vorticity"})
names.append({'class':"VorticityX", 'name':"Vorticity X", 'tag':"vorticityx"})
names.append({'class':"VorticityY", 'name':"Vorticity Y", 'tag':"vorticityy"})
names.append({'class':"VorticityZ", 'name':"Vorticity Z", 'tag':"vorticityz"})

incPath = '../include/QuICC/PhysicalNames/'
srcPath = '../src/PhysicalNames/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n))
    f.close()

