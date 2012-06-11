// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
