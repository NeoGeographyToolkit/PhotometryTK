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


/// \file Weights.h

#ifndef __PHOTOMETRY_WEIGHTS_H__
#define __PHOTOMETRY_WEIGHTS_H__

#include <vw/Math/Vector.h>
#include <string>

namespace photometry {

  void ComputeImageCenterLines(struct ModelParams & modelParams);
  
  int* ComputeImageHCenterLine(std::string input_img_file,
                              int **r_hMaxDistArray);
  int* ComputeDEMHCenterLine(std::string input_DEM_file, int noDataDEMVal,
                            int **r_hMaxDistArray);
  int* ComputeImageVCenterLine(std::string input_img_file,
                                 int **r_hMaxDistArray);
  int* ComputeDEMVCenterLine(std::string input_DEM_file, int noDataDEMVal,
                               int **r_hMaxDistArray);
  float ComputeLineWeightsH(vw::Vector2  const& pix, std::vector<int> const& hCenterLine,
                           std::vector<int> const& hMaxDistArray);
  float ComputeLineWeightsV(vw::Vector2  const& pix, std::vector<int> const& hCenterLine,
                            std::vector<int> const& hMaxDistArray);
  float ComputeLineWeightsHV    (vw::Vector2  const& pix, struct ModelParams const& modelParams);
  void SaveWeightsParamsToFile  (struct ModelParams const& modelParams);
  void ReadWeightsParamsFromFile(struct ModelParams & modelParams);

} // end photometry

#endif//__PHOTOMETRY_WEIGHTS_H__
