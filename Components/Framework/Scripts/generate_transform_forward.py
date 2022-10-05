#!/usr/bin/env python

def header(descr):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief Forward transform operator {name}
 */

#ifndef QUICC_TRANSFORM_FORWARD_{cName.upper()}_HPP
#define QUICC_TRANSFORM_FORWARD_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Forward/IRegisterId.hpp"

namespace QuICC {{

namespace Transform {{

namespace Forward {{

   /**
    * @brief Forward transform operator {name}
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

}} // Forward
}} // Transform
}} // QuICC

#endif // QUICC_TRANSFORM_FORWARD_{cName.upper()}_HPP
'''
    return content

def source(descr):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the forward transform operator {name}
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace Transform {{

namespace Forward {{

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

}} // Forward
}} // Transform
}} // QuICC
'''
    return content

names = []
names.append({'class':"D1", 'name':"Forward::D1", 'tag':"Fwd::D1"})
names.append({'class':"D1ZP0", 'name':"Forward::D1ZP0", 'tag':"Fwd::D1ZP0"})
names.append({'class':"D2", 'name':"Forward::D2", 'tag':"Fwd::D2"})
names.append({'class':"DfLaplh_1", 'name':"Forward::DfLaplh_1", 'tag':"Fwd::DfLaplh_1"})
names.append({'class':"I2D1", 'name':"Forward::I2D1", 'tag':"Fwd::I2D1"})
names.append({'class':"I2P", 'name':"Forward::I2P", 'tag':"Fwd::I2P"})
names.append({'class':"I2Q", 'name':"Forward::I2Q", 'tag':"Fwd::I2Q"})
names.append({'class':"I2S", 'name':"Forward::I2S", 'tag':"Fwd::I2S"})
names.append({'class':"I2T", 'name':"Forward::I2T", 'tag':"Fwd::I2T"})
names.append({'class':"I2ZI2D1", 'name':"Forward::I2ZI2D1", 'tag':"Fwd::I2ZI2D1"})
names.append({'class':"I4D1", 'name':"Forward::I4D1", 'tag':"Fwd::I4D1"})
names.append({'class':"I4D1ZI2", 'name':"Forward::I4D1ZI2", 'tag':"Fwd::I4D1ZI2"})
names.append({'class':"I4P", 'name':"Forward::I4P", 'tag':"Fwd::I4P"})
names.append({'class':"I4Q", 'name':"Forward::I4Q", 'tag':"Fwd::I4Q"})
names.append({'class':"I4R_1D1R1", 'name':"Forward::I4R_1D1R1", 'tag':"Fwd::I4R_1D1R1"})
names.append({'class':"I4R_1D1R1ZI2", 'name':"Forward::I4R_1D1R1ZI2", 'tag':"Fwd::I4R_1D1R1ZI2"})
names.append({'class':"I4R_1Pm", 'name':"Forward::I4R_1Pm", 'tag':"Fwd::I4R_1Pm"})
names.append({'class':"I4S", 'name':"Forward::I4S", 'tag':"Fwd::I4S"})
names.append({'class':"I6Laplh", 'name':"Forward::I6Laplh", 'tag':"Fwd::I6Laplh"})
names.append({'class':"I6LaplhZI4D1R1", 'name':"Forward::I6LaplhZI4D1R1", 'tag':"Fwd::I6LaplhZI4D1R1"})
names.append({'class':"I6R_1D1R1", 'name':"Forward::I6R_1D1R1", 'tag':"Fwd::I6R_1D1R1"})
names.append({'class':"I6R_1D1R1ZI4", 'name':"Forward::I6R_1D1R1ZI4", 'tag':"Fwd::I6R_1D1R1ZI4"})
names.append({'class':"I6R_1Pm", 'name':"Forward::I6R_1Pm", 'tag':"Fwd::I6R_1Pm"})
names.append({'class':"Laplh", 'name':"Forward::Laplh", 'tag':"Fwd::Laplh"})
names.append({'class':"Laplh2", 'name':"Forward::Laplh2", 'tag':"Fwd::Laplh2"})
names.append({'class':"LaplhD1", 'name':"Forward::LaplhD1", 'tag':"Fwd::LaplhD1"})
names.append({'class':"LaplhSin_1", 'name':"Forward::LaplhSin_1", 'tag':"Fwd::LaplhSin_1"})
names.append({'class':"LaplhSin_1Dphi", 'name':"Forward::LaplhSin_1Dphi", 'tag':"Fwd::LaplhSin_1Dphi"})
names.append({'class':"Laplh_1", 'name':"Forward::Laplh_1", 'tag':"Fwd::Laplh_1"})
names.append({'class':"Laplh_1D1", 'name':"Forward::Laplh_1D1", 'tag':"Fwd::Laplh_1D1"})
names.append({'class':"Laplh_1Sin_1", 'name':"Forward::Laplh_1Sin_1", 'tag':"Fwd::Laplh_1Sin_1"})
names.append({'class':"Laplh_1Sin_1Dphi", 'name':"Forward::Laplh_1Sin_1Dphi", 'tag':"Fwd::Laplh_1Sin_1Dphi"})
names.append({'class':"P", 'name':"Forward::P", 'tag':"Fwd::P"})
names.append({'class':"P0", 'name':"Forward::P0", 'tag':"Fwd::P0"})
names.append({'class':"Pm", 'name':"Forward::Pm", 'tag':"Fwd::Pm"})
names.append({'class':"Pol", 'name':"Forward::Pol", 'tag':"Fwd::Pol"})
names.append({'class':"Q", 'name':"Forward::Q", 'tag':"Fwd::Q"})
names.append({'class':"R1", 'name':"Forward::R1", 'tag':"Fwd::R1"})
names.append({'class':"S", 'name':"Forward::S", 'tag':"Fwd::S"})
names.append({'class':"Sin_1", 'name':"Forward::Sin_1", 'tag':"Fwd::Sin_1"})
names.append({'class':"Sin_1Dphi", 'name':"Forward::Sin_1Dphi", 'tag':"Fwd::Sin_1Dphi"})
names.append({'class':"T", 'name':"Forward::T", 'tag':"Fwd::T"})

incPath = '../include/QuICC/Transform/Forward/'
srcPath = '../src/Transform/Forward/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n))
    f.close()

