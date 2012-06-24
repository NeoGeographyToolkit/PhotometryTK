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


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <math.h>
#include <photk/Shadow.h>
#include <photk/Reconstruct.h>
using namespace photometry;

void photometry::ComputeSaveShadowMap( ModelParams input_img_params,
                                           GlobalParams globalParams) {
  std::string shadowMapFile = input_img_params.shadowFilename;
  std::string origfile = input_img_params.inputFilename;
  DiskImageView<PixelMask<PixelGray<uint8> > > originalImage(origfile);
  ImageView<PixelMask<PixelGray<uint8> > > shadowImage(originalImage.cols(), originalImage.rows());

  GeoReference originalGeo;
  read_georeference(originalGeo, origfile);

  for (int k=0; k < (int)originalImage.rows(); ++k) {
    for (int l=0; l < (int)originalImage.cols(); ++l) {


      Vector2 sample_pix(l,k);

      if ( is_valid(originalImage(l,k)) ) {
        shadowImage(l, k) = 0;
        //printf("shadowThresh = %d\n", globalParams.shadowThresh);
        if (originalImage(l, k) < globalParams.shadowThresh){
          shadowImage(l, k) = 255;
        }
      }

    }
  }

  write_georeferenced_image(shadowMapFile,
                            channel_cast<uint8>(shadowImage),
                            originalGeo, TerminalProgressCallback("photometry","Processing:"));
}


//input_img_file is the original image
//output_img_file is the brightness compensated image file with invalid values for shadow
//this is also the filename of the output image where shadows are added
//
void photometry::AddShadows(std::string input_img_file,
                                std::string output_img_file,
                                std::string shadow_file) {
  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  output_img(output_img_file);
  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

  ImageView<PixelMask<PixelGray<uint8> > > r_img (input_img.cols(), input_img.rows());
  int l,k;
  //initialize  output_img, and numSamples
  for (k = 0 ; k < input_img.rows(); ++k) {
    for (l = 0; l < input_img.cols(); ++l) {
      if ( (is_valid(input_img(l,k))) && (shadowImage(l,k) == 255)  ){
        r_img(l,k) = (uint8)(input_img(l,k));
      }
      else{
        r_img(l,k) = (uint8)(output_img(l,k));
      }
    }
  }

  //write in the previous DEM
  write_georeferenced_image(output_img_file,
                            channel_cast<uint8>(r_img),
                            input_img_geo, TerminalProgressCallback("photometry","Processing:"));

}
