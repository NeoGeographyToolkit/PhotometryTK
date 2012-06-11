// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
