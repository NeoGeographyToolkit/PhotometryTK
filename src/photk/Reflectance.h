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


/// \file Reflectance.h

#ifndef __PHOTOMETRY_REFLECTANCE_H__
#define __PHOTOMETRY_REFLECTANCE_H__

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <map>
#include <vw/Math/Vector.h>

#include <photk/Reconstruct.h>

namespace photometry {

  vw::Vector3 computeNormalFrom3DPointsGeneral(vw::Vector3 p1, vw::Vector3 p2, vw::Vector3 p3);
  vw::Vector3 computeNormalFrom3DPoints(vw::Vector3 p1, vw::Vector3 p2, vw::Vector3 p3);

  void ReadSunOrSpacecraftPosition(std::string const& filename,             // Input
                                   std::map<std::string, vw::Vector3> & records // Output
                                   );
  
  float computeReflectanceFromNormalOld(vw::Vector3 sunPos, vw::Vector3 xyz,  vw::Vector3 normal);
  float computeLambertianReflectanceFromNormal(vw::Vector3 sunPos,
                                               vw::Vector3 xyz, vw::Vector3 normal);
  float computeLunarLambertianReflectanceFromNormal(vw::Vector3 const& sunPos,
                                                    vw::Vector3 const& viewPos,
                                                    vw::Vector3 const& xyz,
                                                    vw::Vector3 const& normal,
                                                    float phaseCoeffC1, float phaseCoeffC2,
                                                    float & alpha // output, phase angle
                                                    );
  float ComputeReflectance(vw::Vector3 const& normal, vw::Vector3 const& xyz,
                           ModelParams const& input_img_params,
                           GlobalParams const& globalParams,
                           float & phaseAngle // output
                           );

  void computeReflectanceAux(vw::ImageView<vw::PixelGray<float> > const& DEMTile,
                             vw::cartography::GeoReference const& DEMGeo,
                             float noDEMDataValue,
                             ModelParams const& input_img_params,
                             GlobalParams const& globalParams,
                             vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& outputReflectance,
                             bool savePhaseAngle,
                             vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& phaseAngle);

  void computeReflectanceAtPixel(int x, int y,
                                 vw::ImageView<vw::PixelGray<float> > const& DEMTile,
                                 vw::cartography::GeoReference const& DEMGeo,
                                 float noDEMDataValue,
                                 ModelParams const& input_img_params,
                                 GlobalParams const& globalParams,                                 
                                 bool savePhaseAngle,
                                 float &outputReflectance,
                                 float & phaseAngle);
    
  float actOnImage(std::vector<ImageRecord> & DEMTiles,
                   std::vector<ImageRecord> & albedoTiles,
                   std::vector<ImageRecord> & weightsSumTiles,
                   std::vector<int> & overlap,
                   ModelParams & input_img_params,
                   GlobalParams const& globalParams);

  // Only old code below
  
  float computeLunarLambertianReflectanceFromNormalOld(vw::Vector3 sunPos,
                                                       vw::Vector3 viewerPos,
                                                       vw::Vector3 xyz,
                                                       vw::Vector3 normal,
                                                       float B_0, float L);

  float computeImageReflectanceOld(ModelParams const& input_img_params,
                                   GlobalParams const&  globalParams);
  float computeImageReflectanceOld(ModelParams const& input_img_params,
                                   ModelParams const& overlap_img_params,
                                   GlobalParams const& globalParams);

  void computeXYZandSurfaceNormalOld(vw::ImageView<vw::PixelGray<float> > const& DEMTile,
                                     vw::cartography::GeoReference const& DEMGeo,
                                     float noDEMDataValue,
                                     vw::ImageView<vw::Vector3> & dem_xyz,
                                     vw::ImageView<vw::Vector3> & surface_normal
                                     );

  void computeReflectanceAuxOld(vw::ImageView<vw::Vector3> const& dem_xyz,
                                vw::ImageView<vw::Vector3> const& surface_normal,
                                ModelParams const& input_img_params,
                                GlobalParams const& globalParams,
                                vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& outputReflectance,
                                bool savePhaseAngle,
                                vw::ImageView<vw::PixelMask<vw::PixelGray<float> > > & phaseAngle
                                );
  
  
  float computeImageReflectanceNoWriteOld(ModelParams const& input_img_params,
                                          GlobalParams const& globalParams,
                                          vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& output_img);

}

#endif//__PHOTOMETRY_REFLECTANCE_H__
