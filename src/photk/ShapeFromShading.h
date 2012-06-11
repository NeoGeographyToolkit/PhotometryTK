// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
  void UpdateHeightMap(ModelParams input_img_params, std::vector<ModelParams> overlap_img_params, GlobalParams globalParams);
  
} //end photometry

#endif//__PHOTOMETRY_SHAPE_FROM_SHADING_H__
