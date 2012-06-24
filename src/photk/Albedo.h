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

/// \file Albedo.h

#ifndef __PHOTOMETRY_ALBEDO_H__
#define __PHOTOMETRY_ALBEDO_H__

#include <string>
#include <vector>
#include <photk/Reconstruct.h>

namespace photometry {

  struct phaseCoeffsData{
    double phaseCoeffC1_num, phaseCoeffC1_den;
    double phaseCoeffC2_num, phaseCoeffC2_den;
    phaseCoeffsData(){
      phaseCoeffC1_num = phaseCoeffC1_den = phaseCoeffC2_num = phaseCoeffC2_den = 0.0;
    }
    void writeToFile(std::string fileName){
      std::cout << "Writing to: " << fileName << std::endl;
      std::ofstream fh(fileName.c_str());
      fh.precision(16);
      fh << phaseCoeffC1_num << ' ' << phaseCoeffC1_den << ' '
         << phaseCoeffC2_num << ' ' << phaseCoeffC2_den << std::endl;
      fh.close();
    }
    void readFromFile(std::string fileName){
      std::ifstream fh(fileName.c_str());
      if (!fh || !( fh >> phaseCoeffC1_num >> phaseCoeffC1_den >>
                    phaseCoeffC2_num >> phaseCoeffC2_den) ){
        std::cerr << "Could not read phase data from file: " << fileName << std::endl;
        exit(1);
      }
      fh.close();
    }
  };
  
  double actOnTile(bool isLastIter, bool computeErrors,
                   std::string sampleTileFile,
                   ImageRecord albedoTileCorners,
                   std::string DEMTileFile,
                   std::string albedoTileFile,
                   std::string errorTileFile, std::string weightsSumFile, 
                   std::vector<ModelParams> & overlap_img_params,
                   GlobalParams globalParams,
                   phaseCoeffsData & PCD
                   );
  
  void AppendCostFunToFile(double costFunVal, std::string fileName);
  
  // Only obsolete code below

  //image mosaic functions
  void InitImageMosaic(ModelParams input_img_params,
                       std::vector<ModelParams> overlap_img_params,
                       GlobalParams globalParams);

  void InitImageMosaicByBlocks(ModelParams input_img_params,
                               std::vector<ModelParams> overlap_img_params,
                               GlobalParams globalParams);

  void UpdateImageMosaic(ModelParams input_img_params,
                         std::vector<ModelParams> overlap_img_params,
                         GlobalParams globalParams);

  //albedo mosaic functions
  void InitAlbedoMosaic(ModelParams input_img_params,
                        std::vector<ModelParams> overlap_img_params,
                        GlobalParams globalParams);

  //albedo mosaic functions
  void UpdateAlbedoMosaic(ModelParams input_img_params,
                          std::vector<ModelParams> overlap_img_params,
                          GlobalParams globalParams);

} // end photometry

#endif//__PHOTOMETRY_ALBEDO_H__
