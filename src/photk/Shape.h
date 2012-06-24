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


/// \file Shape.h

#ifndef __PHOTOMETRY_SHAPE_H__
#define __PHOTOMETRY_SHAPE_H__

#include <iostream>
#include <fstream>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>

#include <photk/Reconstruct.h>

namespace photometry {

  void InitMeanDEMTile(std::string sampleTile,
                       ImageRecord currTileCorners,
                       std::string meanDEMTileFile,
                       std::vector<ImageRecord> & DEMImages,
                       std::vector<int> & overlap,
                       GlobalParams globalParams);

  // Only old code below
  
  void InitDEM(ModelParams input_img_params,
               std::vector<ModelParams> overlap_img_params,
               GlobalParams globalParams);
  vw::Vector3 computeNormalFrom3DPoints(vw::Vector3 p1, vw::Vector3 p2, vw::Vector3 p3);

} // end photometry

#endif//__PHOTOMETRY_SHAPE_H__
