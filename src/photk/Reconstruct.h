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

/// \file Reconstruct.h

#ifndef __PHOTOMETRY_RECONSTRUCT_H__
#define __PHOTOMETRY_RECONSTRUCT_H__

enum {NO_REFL = 0, LAMBERT, LUNAR_LAMBERT};
enum {NO_SHADOW_REMOVAL = 0,
      CONSTANT_THRESHOLD_SHADOW_REMOVAL,
      LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL,
      LUNAR_LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL
};

#include <iostream>
#include <string>
#include <vw/Math/Vector.h>

struct ImageRecord {
  static const double defaultCoord = -999.0;
  bool useImage;
  std::string path;
  double north, west, south, east;
  
  ImageRecord(void) {}
  ImageRecord(const std::string& path_) :
    useImage(true),
    path(path_),
    north(defaultCoord)
  {}
};
 

namespace photometry {

  struct GlobalParams{

    std::string drgDir;
    std::string demDir;
    std::string sunPosFile;
    std::string spacecraftPosFile;

    int initialSetup;
    float tileSize;        // in degrees
    int pixelPadding;      // in pixels
    vw::Vector4 simulationBox; // lonMin, lonMax, latMin, latMax (If not present the entire albedo will be simulated)
    
    int reflectanceType;
    int initDEM;
    int initExposure;
    int initAlbedo;
    int shadowType;

    float shadowThresh;

    //float exposureInitRefValue;//this will be removed
    //int exposureInitRefIndex;//this will be removed
    float TRConst;
    int updateAlbedo, updateExposure, updateHeight;

    // Two parameters used in the formula for the reflectance
    float phaseCoeffC1, phaseCoeffC2;
    // Update the components of the coefficients phaseCoeffC1 and
    // phaseCoeffC2 for each tile.
    int updateTilePhaseCoeffs;
    // Update the phase coefficients by combining the results from all tiles
    int updatePhaseCoeffs;
    
    int useWeights;
    int saveWeights, computeWeightsSum, useNormalizedCostFun;
    int maxNumIter;
    int computeErrors;
    int noDEMDataValue;
    int forceMosaic; // see the description in reconstruct.cc
  };

  struct ModelParams {

    float   exposureTime;

    vw::Vector2 cameraParams; //currently not used
    vw::Vector3 sunPosition; //relative to the center of the Moon
    vw::Vector3 spacecraftPosition;//relative to the center of the planet

    std::vector<int> hCenterLine;
    std::vector<int> hMaxDistArray;
    std::vector<int> vCenterLine;
    std::vector<int> vMaxDistArray;

    int *hCenterLineDEM;
    int *hMaxDistArrayDEM;
    int *vCenterLineDEM;
    int *vMaxDistArrayDEM;
    
    vw::Vector4 corners; // cached bounds to quickly calculate overlap
    
    /*
    vector<int> hCenterLine;
    vector<int> hMaxDistArray;
    vector<int> hCenterLineDEM;
    vector<int> hMaxDistArrayDEM;
    vector<int> vCenterLine;
    vector<int> vMaxDistArray;
    vector<int> vCenterLineDEM;
    vector<int> vMaxDistArrayDEM;
    */
    std::string infoFilename, DEMFilename, meanDEMFilename,
      var2DEMFilename, reliefFilename, shadowFilename,
      errorFilename, inputFilename, outputFilename, 
      sfsDEMFilename, errorHeightFilename, weightFilename, exposureFilename;

    ModelParams(){
      hCenterLineDEM   = NULL;
      hMaxDistArrayDEM = NULL;
      vCenterLineDEM   = NULL;
      vMaxDistArrayDEM = NULL;
    }

    ~ModelParams(){
      // Need to deal with copying of these structures properly before deleting
      // delete[] hCenterLine;       hCenterLine       = NULL;
      // delete[] hMaxDistArray;     hMaxDistArray     = NULL;
      // delete[] hCenterLineDEM;    hCenterLineDEM    = NULL;
      // delete[] hMaxDistArrayDEM;  hMaxDistArrayDEM  = NULL;
      // delete[] vCenterLine;       vCenterLine       = NULL;
      // delete[] vMaxDistArray;     vMaxDistArray     = NULL;
      // delete[] vCenterLineDEM;    vCenterLineDEM    = NULL;
      // delete[] vMaxDistArrayDEM;  vMaxDistArrayDEM  = NULL;
    }
    
  };

  // Generic Ostream options for Debugging
  std::ostream& operator<<( std::ostream& os, GlobalParams const& global );
  std::ostream& operator<<( std::ostream& os, ModelParams const& model );

  // Written by Taemin Kim - START
  #define DYNAMIC_RANGE  256
  // Written by Taemin Kim - END

} // end photometry

#endif//__PHOTOMETRY_RECONSTRUCT_H__
