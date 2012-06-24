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


/// \file ShapeFromShading.h

#ifndef __PHOTOMETRY_SHAPE_FROM_SHADING_H__
#define __PHOTOMETRY_SHAPE_FROM_SHADING_H__

#include <vw/Image/ImageView.h>
#include <vw/Math/Vector.h>
#include <photk/Reconstruct.h>

namespace photometry {

  void UpdateHeightMapTiles(std::string DEMTileFile,
                            std::string albedoTileFile,
                            std::string sfsTileFile,
                            std::vector<ModelParams> & overlapImgParams,
                            GlobalParams globalParams
                            );

  // Does conjugate gradient descent on the DEM, keeping all else fixed.
  void optimize_conjugate_gradient(vw::ImageView<vw::PixelGray<double> > *image_predicted,
                                   vw::ImageView<vw::PixelGray<double> > *image,
                                   vw::ImageView<vw::PixelGray<double> > *dem,
                                   vw::ImageView<vw::PixelGray<double> > *init_dem,
                                   vw::ImageView<vw::PixelGray<double> > *albedo,
                                   vw::Vector3 *light_direction);

  // Old code, using images, not tiles
  void UpdateHeightMapOld(ModelParams input_img_params, std::vector<ModelParams> overlap_img_params, GlobalParams globalParams);
  
} //end photometry

#endif//__PHOTOMETRY_SHAPE_FROM_SHADING_H__
