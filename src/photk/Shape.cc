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

#include <iostream>

#include <math.h>
#include <time.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <photk/Shape.h>
#include <photk/Reconstruct.h>
#include <photk/Reflectance.h>
#include <photk/Weights.h>
#include <photk/Misc.h>
using namespace photometry;

//initializes the DEM tile
void photometry::InitMeanDEMTile(std::string sampleTileFile,
                                 ImageRecord currTileCorners,
                                 std::string meanDEMTileFile,
                                 std::vector<ImageRecord> & DEMImages,
                                 std::vector<int> & overlap,
                                 GlobalParams globalParams) {

  // For a given tile, find all the input DEM tiles overlapping with
  // it, combine them into one combined DEM, and then interpolate that
  // combined DEM at all pixels of the desired tile.

  // To make things a bit more efficient, when creating the combined
  // DEM, we take only pixels not too far from the desired tile.

  // Try to read the noDEMDataValue from the DEM images. If that fails, use
  // the value provided in the parameter file.
  float noDEMDataValue;
  if ( overlap.size() == 0 || !readNoDEMDataVal(DEMImages[overlap[0]].path, noDEMDataValue)){
    noDEMDataValue = globalParams.noDEMDataValue;
  }
  std::cout << "using noDEMDataValue: " << noDEMDataValue << std::endl;

  // Create the DEM tile based on the current tile corners, and the sample tile
  // from which we extract and then adjust the georeference.
  std::cout << "Reading " << sampleTileFile << std::endl;
  GeoReference DEMTileGeo;
  read_georeference(DEMTileGeo, sampleTileFile);
  ImageView<PixelGray<float> > meanDEMTile;
  createGeoRefAndTileWithGivenCorners(currTileCorners.west, currTileCorners.east,
                                      currTileCorners.south, currTileCorners.north,
                                      DEMTileGeo, meanDEMTile // outputs
                                      );
  for (int k = 0 ; k < (int)meanDEMTile.rows(); ++k) {
    for (int l = 0; l < (int)meanDEMTile.cols(); ++l) {
      meanDEMTile(l, k) = noDEMDataValue;
    }
  }

//  system("echo top1 is $(top -u $(whoami) -b -n 1|grep reconstruct)");

  // Find the combined DEM image as mentioned earlier.
  ImageView<PixelGray<float> > combined_DEM;
  GeoReference combined_DEM_geo;
  Vector2 begTileLonLat = DEMTileGeo.pixel_to_lonlat(Vector2(0, 0));
  Vector2 endTileLonLat = DEMTileGeo.pixel_to_lonlat(Vector2(meanDEMTile.cols()-1, meanDEMTile.rows()-1));
  std::vector<std::string> overlapDEMVec;
  for (int i = 0; i < (int)overlap.size(); i++){
    overlapDEMVec.push_back(DEMImages[overlap[i]].path);
  }
  readDEMTilesIntersectingBox(noDEMDataValue, begTileLonLat, endTileLonLat, overlapDEMVec, // Inputs
                              combined_DEM, combined_DEM_geo                               // Outputs
                              );
  
//  system("echo top2 is $(top -u $(whoami) -b -n 1|grep reconstruct)");

  InterpolationView<EdgeExtensionView<ImageView<PixelGray<float> >, ConstantEdgeExtension>, BicubicInterpolation>
    interp_combined_DEM = interpolate(combined_DEM, BicubicInterpolation(), ConstantEdgeExtension());
  // Wrong below
  //ImageViewRef<PixelGray<float> >  interp_combined_DEM
  // = interpolate(edge_extend(combined_DEM.impl(),  ConstantEdgeExtension()),
  //BilinearInterpolation());

  // Interpolate
  for (int k = 0 ; k < (int)meanDEMTile.rows(); ++k) {
    for (int l = 0; l < (int)meanDEMTile.cols(); ++l) {
      
      Vector2 input_DEM_pix(l,k);
      
      //check for overlap between the output image and the input DEM image
      Vector2 combined_pix = combined_DEM_geo.lonlat_to_pixel(DEMTileGeo.pixel_to_lonlat(input_DEM_pix));
      float x = combined_pix[0];
      float y = combined_pix[1];
        
      //check for valid DEM coordinates
      if ((x>=0) && (x <= combined_DEM.cols()-1) && (y>=0) && (y<= combined_DEM.rows()-1)){

        // Check that all four grid points used for interpolation are valid
        if ( combined_DEM( floor(x), floor(y) ) != noDEMDataValue &&
             combined_DEM( floor(x), ceil(y)  ) != noDEMDataValue &&
             combined_DEM( ceil(x),  floor(y) ) != noDEMDataValue &&
             combined_DEM( ceil(x),  ceil(y)  ) != noDEMDataValue
             ){
          meanDEMTile(l, k) = (float)interp_combined_DEM(x, y);
        }
      }
    }
  }
  
//  system("echo top3 is $(top -u $(whoami) -b -n 1|grep reconstruct)");

  // Write the DEM together with the noDEMDataValue
  std::cout << "Writing: " << meanDEMTileFile << std::endl;
  DiskImageResourceGDAL::Options gdal_options;
  gdal_options["COMPRESS"] = "LZW";
  DiskImageResourceGDAL rsrc(meanDEMTileFile, meanDEMTile.format(), Vector2i(256, 256), gdal_options);
  rsrc.set_nodata_write(noDEMDataValue);
  write_georeference(rsrc, DEMTileGeo);
  write_image(rsrc, meanDEMTile, TerminalProgressCallback("{Core}", "Processing:"));
  
  return;
}

// Only old code below

float ComputeGradient_DEM(float /*intensity*/, float T,
                          float albedo, Vector3 s,
                          Vector3 p, Vector3 p_left,
                          Vector3 p_top, Vector3 /*xyz_prior*/) {
  float grad;
  float temp  = (p[2]-p_left[2])*(p[2]-p_left[2]) + (p[2]-p_top[2])*(p[2]-p_top[2]) + 1;
  float temp1 = (p[2]-p_left[2])*s[0]+(p[2]-p_top[2])*s[1] +s[2];
  float temp2 = (p[2]-p_left[2])+(p[2]-p_top[2]);
  grad = T*albedo*((s[0]+s[1])*sqrt(temp)-(temp1*temp2/sqrt(temp)))/temp;

  return grad;
}

float ComputeError_DEM(float intensity, float T, float albedo,
                       float reflectance, Vector3 /*xyz*/, Vector3 xyz_prior) {
  float error;
  error = (intensity-T*albedo*reflectance) + (xyz_prior[2]-xyz_prior[2]);
  return error;
}

//initializes the DEM file by getting the average DEM values in the overlapping areas of consecutive DEM files.
void photometry::InitDEM( ModelParams input_img_params,
                              std::vector<ModelParams> overlap_img_params,
                              GlobalParams globalParams) {

  int i;
  unsigned l, k;

  std::string input_DEM_file = input_img_params.DEMFilename;
  std::string mean_DEM_file  = input_img_params.meanDEMFilename;
  std::string var2_DEM_file  = input_img_params.var2DEMFilename;

  DiskImageView<PixelGray<float> >  input_DEM_image(input_DEM_file);
  GeoReference input_DEM_geo;
  read_georeference(input_DEM_geo, input_DEM_file);

  ImageView<PixelGray<float> > mean_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());
  ImageView<PixelMask<PixelGray<float> > >var2_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());

  ImageView<PixelGray<int> > numSamples(input_DEM_image.cols(), input_DEM_image.rows());
  ImageView<PixelGray<float> > norm(input_DEM_image.cols(), input_DEM_image.rows());

  printf("corrections to var2\n");
  //initialize  mean_DEM-image, var2_DEM_image and numSamples and norm
  for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
    for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

      mean_DEM_image(l, k) = globalParams.noDEMDataValue; //-10000;
      var2_DEM_image(l, k) = 0;
      numSamples(l, k) = 0;
      norm(l,k) = 0;

      Vector2 input_DEM_pix(l,k);

      if ( input_DEM_image(l,k) != globalParams.noDEMDataValue ) {

        if (globalParams.useWeights == 0){
          mean_DEM_image(l, k) = (float)input_DEM_image(l,k);
          var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k);
          numSamples(l, k) = 1;
        }
        else{
          float weight = ComputeLineWeightsHV(input_DEM_pix, input_img_params);
          mean_DEM_image(l, k) = (float)input_DEM_image(l,k)*weight;
          //weight added by Ara 08/28
          var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k)*weight;
          //numSamples(l, k) = 1;
          norm(l, k) = weight;
        }
      }
    }
  }

  for (i = 0; i < (int) overlap_img_params.size(); i++){

    printf("DEM = %s\n", overlap_img_params[i].DEMFilename.c_str());
    DiskImageView<PixelGray<float> > overlap_DEM_image(overlap_img_params[i].DEMFilename);

    InterpolationView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_overlap_DEM_image = interpolate(overlap_DEM_image, BilinearInterpolation(), ConstantEdgeExtension());
    // Wrong
    //ImageViewRef<PixelGray<float> >  interp_overlap_DEM_image
    //= interpolate(edge_extend(overlap_DEM_image.impl(), ConstantEdgeExtension()), BilinearInterpolation());
    
    GeoReference overlap_DEM_geo;
    read_georeference(overlap_DEM_geo, overlap_img_params[i].DEMFilename);


    for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
      for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

        Vector2 input_DEM_pix(l,k);

        if ( input_DEM_image(l,k) != globalParams.noDEMDataValue ) {

          //check for overlap between the output image and the input DEM image
          Vector2 overlap_dem_pix = overlap_DEM_geo.lonlat_to_pixel(input_DEM_geo.pixel_to_lonlat(input_DEM_pix));
          float x = overlap_dem_pix[0];
          float y = overlap_dem_pix[1];

          //check for valid DEM coordinates
          if ((x>=0) && (x <= overlap_DEM_image.cols()-1) && (y>=0) && (y<= overlap_DEM_image.rows()-1)){

            // Check that all four grid points used for interpolation are valid
            if ( overlap_DEM_image( floor(x), floor(y) ) != globalParams.noDEMDataValue &&
                 overlap_DEM_image( floor(x), ceil(y)  ) != globalParams.noDEMDataValue &&
                 overlap_DEM_image( ceil(x),  floor(y) ) != globalParams.noDEMDataValue &&
                 overlap_DEM_image( ceil(x),  ceil(y)  ) != globalParams.noDEMDataValue
                 ){
              
              if (globalParams.useWeights == 0){
                mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y);
                var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y);
                numSamples(l, k) = (int)numSamples(l, k) + 1;
              }
              else{
                float weight = ComputeLineWeightsHV(overlap_dem_pix, overlap_img_params[i]);
                mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*weight;
                //weight added by Ara 08/28/
                var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y)*weight;
                norm(l, k) = norm(l,k) + weight;
              }
            }
          }
        }
      }
    }
  }


  for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
    for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

      //compute variance only where the mean DEM is valid
      if ( input_DEM_image(l,k) != globalParams.noDEMDataValue ) {

        if ((globalParams.useWeights == 0) && (numSamples(l,k)!=0)){
          mean_DEM_image(l, k) = mean_DEM_image(l, k)/numSamples(l,k);
          var2_DEM_image(l, k) = var2_DEM_image(l, k)/numSamples(l, k) - mean_DEM_image(l,k)*mean_DEM_image(l,k);
        }

        if ((globalParams.useWeights == 1) && (norm(l,k)!=0)){
          mean_DEM_image(l, k) = mean_DEM_image(l, k)/norm(l,k);
          var2_DEM_image(l, k) = var2_DEM_image(l, k)/norm(l, k) - mean_DEM_image(l,k)*mean_DEM_image(l,k);
          if (var2_DEM_image(l, k) < 0){ //this should never happen
              var2_DEM_image(l, k) = 0;
              var2_DEM_image(l, k).invalidate();
          }
        }
      }

    }
  }

  std::cout << "Writing " << mean_DEM_file << std::endl;
  write_georeferenced_image(mean_DEM_file,
                            mean_DEM_image,
                            input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));

  write_georeferenced_image(var2_DEM_file,
                            var2_DEM_image,
                            input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));

}



//initializes the DEM file by getting the average DEM values in the overlapping areas of consecutive DEM files.
void DetectDEMOutliers( std::string input_DEM_file,
                        std::string /*mean_DEM_file*/,
                        std::string var2_DEM_file,
                        ModelParams input_img_params,
                        std::vector<std::string> overlap_DEM_files,
                        std::vector<ModelParams> /*overlap_img_params*/,
                        GlobalParams /*globalParams*/) {

    DiskImageView<PixelGray<float> >  input_DEM_image(input_DEM_file);
    GeoReference input_DEM_geo;
    read_georeference(input_DEM_geo, input_DEM_file);

    DiskImageView<PixelGray<float> >  mean_DEM_image(input_DEM_file);

    ImageView<PixelMask<PixelGray<float> > >var2_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());
    ImageView<PixelGray<int> > numSamples(input_DEM_image.cols(), input_DEM_image.rows());

    float avgStdDevDEM;

    float *meanDEMOffset = new float[5];
    //read the meanDEMOffset and avgStdDevDEM from file
    FILE *fp = fopen(input_img_params.infoFilename.c_str(),"r");

    fscanf(fp, "%f %f %f %f %f\n", &meanDEMOffset[0], &meanDEMOffset[1], &meanDEMOffset[2], &meanDEMOffset[3], &avgStdDevDEM);
    fclose(fp);


    //initialize  mean_DEM-image, var2_DEM_image and numSamples
    for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
      for (int32 l = 0; l < input_DEM_image.cols(); ++l) {

           numSamples(l, k) = 0;
           Vector2 input_DEM_pix(l,k);

           if ( input_DEM_image(l,k) != -10000 ) {
               var2_DEM_image(l, k) = ((float)input_DEM_image(l,k)-meanDEMOffset[0]) *((float)input_DEM_image(l,k)-meanDEMOffset[0]);
               numSamples(l, k) = 1;
           }

        }
    }

    for (size_t i = 0; i < overlap_DEM_files.size(); i++){

      printf("DEM = %s\n", overlap_DEM_files[i].c_str());

      DiskImageView<PixelGray<float> >  overlap_DEM_image(overlap_DEM_files[i]);
      GeoReference overlap_DEM_geo;
      read_georeference(overlap_DEM_geo, overlap_DEM_files[i]);


      // This is the wrong way of doing interpolation. The the Stereo module for  the right way.
      ImageViewRef<PixelGray<float> >  interp_overlap_DEM_image = interpolate(edge_extend(overlap_DEM_image.impl(),
                                                                              ConstantEdgeExtension()),
                                                                              BilinearInterpolation());

      for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
        for (int32 l = 0; l < input_DEM_image.cols(); ++l) {

          Vector2 input_DEM_pix(l,k);

          if ( input_DEM_image(l,k) != -10000 ) {

              //check for overlap between the output image and the input DEM image
              Vector2 overlap_dem_pix = overlap_DEM_geo.lonlat_to_pixel(input_DEM_geo.pixel_to_lonlat(input_DEM_pix));
              int x = (int)overlap_dem_pix[0];
              int y = (int)overlap_dem_pix[1];

              //check for valid DEM coordinates
              if ((x>=0) && (x < overlap_DEM_image.cols()) && (y>=0) && (y< overlap_DEM_image.rows())){

                if ( overlap_DEM_image(x, y) != -10000 ) {

                    var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + ((float)interp_overlap_DEM_image(x, y)-meanDEMOffset[i+1])*
                                                                         ((float)interp_overlap_DEM_image(x, y)-meanDEMOffset[i+1]);
                    numSamples(l, k) = (int)numSamples(l, k) + 1;
                }
             }
          }
        }
      }
    }


    //compute mean and variance
    int totalNumSamples = 0;
    avgStdDevDEM = 0.0;
    for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
      for (int32 l = 0; l < input_DEM_image.cols(); ++l) {
         if (numSamples(l,k)!=0){

            var2_DEM_image(l, k) = var2_DEM_image(l, k)/numSamples(l,k) - mean_DEM_image(l, k)*mean_DEM_image(l,k);

            //compute the DEM standard deviation
            var2_DEM_image(l, k) = sqrt((float)var2_DEM_image(l,k));
            avgStdDevDEM = avgStdDevDEM + var2_DEM_image(l,k)/numSamples(l,k);
            totalNumSamples++;
         }
       }
    }

    //printf("average DEM error = %f\n", totalVar2/(float)numSamples);

    //write the DEM
    write_georeferenced_image(var2_DEM_file,
                              //channel_cast<float>(var2_DEM_image),
                              channel_cast<uint8>(clamp(var2_DEM_image,0.0,255.0)),
                              //var2_DEM_image,
                              input_DEM_geo, TerminalProgressCallback("photometry","Processing:"));


    //write_georeferenced_image(output_file, tm_image, geo1, TerminalProgressCallback());
}

