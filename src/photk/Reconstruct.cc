//__BEGIN_LICENSE__
//  Copyright (c) 2009-2012, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/// \file Reconstruct.cc
///

#include <photk/Reconstruct.h>
using namespace photometry;

// Generic Ostream options for Debugging
std::ostream& photometry::operator<<( std::ostream& os,
                                          GlobalParams const& global ) {

  os << "-- Global Params --\n";
  os << " ReflectanceType   : " << global.reflectanceType   << "\n";
  os << " Shadow threshold  : " << global.shadowThresh      << "\n";
  os << " SpacecraftPosFile : " << global.spacecraftPosFile << "\n";
  os << " SunPosFile        : " << global.sunPosFile        << "\n";

  return os;
}

std::ostream& photometry::operator<<( std::ostream& os,
                                          ModelParams const& model ) {

  os << "-- Model Params --\n";
  os << " Exposure Time       : " << model.exposureTime        << "\n";
  os << " Sun Position        : " << model.sunPosition         << "\n";
  os << " Spacecraft Position : " << model.spacecraftPosition  << "\n";
  os << " Info File           : " << model.infoFilename        << "\n";
  os << " DEMFilename         : " << model.DEMFilename         << "\n";
  os << " meanDEMFilename     : " << model.meanDEMFilename     << "\n";
  os << " sfsDEMFilename      : " << model.sfsDEMFilename      << "\n";
  os << " var2DEMFilename     : " << model.var2DEMFilename     << "\n";
  os << " reliefFilename      : " << model.reliefFilename      << "\n";
  os << " shadowFilename      : " << model.shadowFilename      << "\n";
  os << " errorFilename       : " << model.errorFilename       << "\n";
  os << " errorHeightFilename : " << model.errorHeightFilename << "\n";
  os << " inputFilename       : " << model.inputFilename       << "\n";
  os << " outputFilename      : " << model.outputFilename      << "\n";

  return os;
}
