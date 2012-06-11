// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
  
  float computeReflectanceFromNormal(vw::Vector3 sunPos, vw::Vector3 xyz,  vw::Vector3 normal);
  float computeLambertianReflectanceFromNormal(vw::Vector3 sunPos,
                                               vw::Vector3 xyz, vw::Vector3 normal);
  float computeLunarLambertianReflectanceFromNormalOld(vw::Vector3 sunPos,
                                                    vw::Vector3 viewerPos,
                                                    vw::Vector3 xyz,
                                                    vw::Vector3 normal,
                                                    float B_0, float L);
  float computeLunarLambertianReflectanceFromNormal(vw::Vector3 sunPos,
                                                    vw::Vector3 viewPos,
                                                    vw::Vector3 xyz,
                                                    vw::Vector3 normal,
                                                    float phaseCoeffA1, float phaseCoeffA2,
                                                    float & alpha // output, phase angle
                                                    );
  float computeImageReflectance(ModelParams const& input_img_params,
                                GlobalParams const&  globalParams);
  float ComputeReflectance(vw::Vector3 normal, vw::Vector3 xyz,
                           ModelParams const& input_img_params,
                           GlobalParams const& globalParams,
                           float & phaseAngle // output
                           );
  float computeImageReflectance(ModelParams const& input_img_params,
                                ModelParams const& overlap_img_params,
                                GlobalParams const& globalParams);

  void computeXYZandSurfaceNormal(vw::ImageView<vw::PixelGray<float> > const& DEMTile,
                                  vw::cartography::GeoReference const& DEMGeo,
                                  float noDEMDataValue,
                                  vw::ImageView<vw::Vector3> & dem_xyz,
                                  vw::ImageView<vw::Vector3> & surface_normal
                                  );

  void computeReflectanceAux(vw::ImageView<vw::Vector3> const& dem_xyz,
                             vw::ImageView<vw::Vector3> const& surface_normal,
                             ModelParams const& input_img_params,
                             GlobalParams const& globalParams,
                             vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& outputReflectance,
                             bool savePhaseAngle,
                             vw::ImageView<vw::PixelMask<vw::PixelGray<float> > > & phaseAngle
                             );


  float actOnImage(std::vector<ImageRecord> & DEMTiles,
                   std::vector<ImageRecord> & albedoTiles,
                   std::vector<ImageRecord> & weightsSumTiles,
                   std::vector<int> & overlap,
                   ModelParams & input_img_params,
                   GlobalParams const& globalParams);
  
  float computeImageReflectanceNoWrite(ModelParams const& input_img_params,
                                       GlobalParams const& globalParams,
                                       vw::ImageView<vw::PixelMask<vw::PixelGray<float> > >& output_img);

}

#endif//__PHOTOMETRY_REFLECTANCE_H__
