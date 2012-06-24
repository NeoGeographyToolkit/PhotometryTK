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


/// \file ReconstructError.h

#ifndef __PHOTOMETRY_RECONSTRUCTERROR_H__
#define __PHOTOMETRY_RECONSTRUCTERROR_H__

#include <string>
#include <vector>
//#include <photk/Reconstruct.h>

namespace photometry {

  float ComputeError(float intensity, float T,
                          float albedo, float reflectance);
                          //Vector3 /*xyz*/, Vector3 /*xyz_prior*/)

  //reconstruction error functions
  void ComputeReconstructionErrorMap(ModelParams input_img_params,
                                     std::vector<ModelParams> overlap_img_params,
                                     GlobalParams globalParams,
                                     float *avgError, int *totalNumSamples);

} // end photometry

#endif//__PHOTOMETRY_RECONSTRUCTERROR_H__
