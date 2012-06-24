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


/// \file Exposure.h

#ifndef __PHOTOMETRY_EXPOSURE_H__
#define __PHOTOMETRY_EXPOSURE_H__

#include <string>
#include <vector>
#include <photk/Reconstruct.h>

namespace photometry {

  //computes the exposure time from image, albedo and DEM
  void ComputeExposure(ModelParams *currModelParams,
                       GlobalParams globalParams);

  //used for mosaicking, with no reflectance model
  void ComputeExposureAlbedo(ModelParams *currModelParams,
                             GlobalParams globalParams);
  /*
  void AppendExposureInfoToFile(std::string exposureFilename,
                                ModelParams currModelParams);
  std::vector<float> ReadExposureInfoFile(std::string exposureFilename,
                                          int numEntries);
  */
  void AppendExposureInfoToFile(ModelParams modelParams);
  void ReadExposureInfoFromFile(ModelParams *modelParams);

} // end photometry

#endif//__PHOTOMETRY_EXPOSURE_H__
