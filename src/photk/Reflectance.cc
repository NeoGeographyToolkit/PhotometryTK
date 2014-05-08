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
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <photk/Reflectance.h>
#include <photk/ShapeFromShading.h>
#include <photk/Shape.h>
#include <photk/Albedo.h>
#include <photk/Exposure.h>
#include <photk/Reconstruct.h>
#include <photk/Misc.h>
#include <photk/Weights.h>
using namespace photometry;
using namespace std;

// Read from file the 3D position of the sun or spacecraft in Moon-centric coordinates.
// We assume that each line in the file has the form:
// AS15-M-0073 -1353012.8795685 151848979.22428 689542.65035626
// We will use the 11-character string "AS15-M-0073" as a unique identifier.
void
photometry::ReadSunOrSpacecraftPosition(std::string const& filename,             // Input
                                        std::map<std::string, Vector3> & records // Output
                                        ){

  records.clear();
  
  std::ifstream infile( filename.c_str() );
  if ( !infile.is_open() ) {
    std::cerr << "ERROR: Could not read file: " << filename << std::endl;
    exit(1);
  }

  while ( infile.good() ){
    std::string line, key;
    Vector3 val;
    getline (infile, line);
    std::istringstream lh(line);
    if (line == "") continue;
    if (! ( lh >> key >> val[0] >> val[1] >> val[2] ) ){
      std::cerr << "ERROR: Unable to read from file: " << filename << " the line: '" << line << "'" << std::endl;
      exit(1);
    }

    if (records.find(key) != records.end()){
      std::cerr << "ERROR: Duplicate key: " << key << " in file: " << filename << std::endl;
      exit(1);
    }
    
    records[key] = val;
  }
  
  infile.close();

  return;
}

//computes the Lambertian reflectance model (cosine of the light direction and the normal to the Moon)
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float
photometry::computeLambertianReflectanceFromNormal(Vector3 sunPos,
                                                   Vector3 xyz,
                                                   Vector3 normal) {
  float reflectance;

  //Vector3 xyz = get_normal(lon_lat);
  //printf("xyz[0] = %f, xyz[1] = %f, xyz[2] = %f\n", xyz[0], xyz[1], xyz[2]);
  // sun coordinates relative to the xyz point on the Moon surface
  Vector3 sunDirection = normalize(sunPos-xyz);

  reflectance = sunDirection[0]*normal[0] + sunDirection[1]*normal[1] + sunDirection[2]*normal[2];

  return reflectance;
}


float
photometry::computeLunarLambertianReflectanceFromNormal(Vector3 const& sunPos, Vector3 const& viewPos, Vector3 const& xyz,
                                                        Vector3 const& normal, float phaseCoeffC1, float phaseCoeffC2,
                                                        float & alpha // output
                                                        ) {
  
  float reflectance;
  float L;

  double len = dot_prod(normal, normal);
  if (abs(len - 1.0) > 1.0e-4){
    std::cerr << "Error: Expecting unit normal in the reflectance computation, in "
              << __FILE__ << " at line " << __LINE__ << std::endl;
    exit(1);
  }

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  //sun coordinates relative to the xyz point on the Moon surface
  //Vector3 sunDirection = -normalize(sunPos-xyz);
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  double tol = 0.3;
  if (mu_0 < tol){
    // Sun is too low, reflectance is too close to 0, the albedo will be inaccurate
    return 0.0;
  }
  
  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  //viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  //compute the phase angle (alpha) between the viewing direction and the light source direction
  float deg_alpha;
  float cos_alpha;

  cos_alpha = dot_prod(sunDirection,viewDirection);
  if ((cos_alpha > 1)||(cos_alpha< -1)){
    printf("cos_alpha error\n");
  }

  alpha     = acos(cos_alpha);  // phase angle in radians
  deg_alpha = alpha*180.0/M_PI; // phase angle in degrees

  //printf("deg_alpha = %f\n", deg_alpha);

  //Bob Gaskell's model
  //L = exp(-deg_alpha/60.0);

  //Alfred McEwen's model
  float A = -0.019;
  float B =  0.000242;//0.242*1e-3;
  float C = -0.00000146;//-1.46*1e-6;

  L = 1.0 + A*deg_alpha + B*deg_alpha*deg_alpha + C*deg_alpha*deg_alpha*deg_alpha;

  //printf(" deg_alpha = %f, L = %f\n", deg_alpha, L);

  if (mu_0 < 0.0){
    return 0.0;
  }

  if (mu < 0.0){ //emission angle is > 90
    mu = 0.0;
  }

  if (mu_0 + mu == 0){
    //printf("negative reflectance\n");
    return 0.0;
  }
  else{
    reflectance = 2*L*mu_0/(mu_0+mu) + (1-L)*mu_0;
  }
  if (reflectance <= 0){
    //printf("negative reflectance\n");
    return 0.0;
  }

  // Attempt to compensate for points on the terrain being too bright
  // if the sun is behind the spacecraft as seen from those points.

  //reflectance *= std::max(0.4, exp(-alpha*alpha));
  reflectance *= ( exp(-phaseCoeffC1*alpha) + phaseCoeffC2 );
  
  return reflectance;
}

float
photometry::ComputeReflectance(Vector3 const& normal, Vector3 const& xyz,
                               ModelParams const& input_img_params,
                               GlobalParams const& globalParams,
                               float & phaseAngle // output
                               ) {
  float input_img_reflectance;

  switch ( globalParams.reflectanceType )
    {
    case LUNAR_LAMBERT:
      //printf("Lunar Lambert\n");
      input_img_reflectance = computeLunarLambertianReflectanceFromNormal(input_img_params.sunPosition,
                                                                          input_img_params.spacecraftPosition,
                                                                          xyz,  normal,
                                                                          globalParams.phaseCoeffC1,
                                                                          globalParams.phaseCoeffC2,
                                                                          phaseAngle // output
                                                                          );
      break;
    case LAMBERT:
      //printf("Lambert\n");
      input_img_reflectance = computeLambertianReflectanceFromNormal(input_img_params.sunPosition,
                                                                     xyz,  normal);
      break;

    default:
      //printf("No reflectance model\n");
      input_img_reflectance = 1;
    }

  return input_img_reflectance;

}

void photometry::computeReflectanceAtPixel(int x, int y,
                                           ImageView<PixelGray<float> > const& DEMTile,
                                           cartography::GeoReference const& DEMGeo,
                                           float noDEMDataValue,
                                           ModelParams const& input_img_params,
                                           GlobalParams const& globalParams,
                                           bool savePhaseAngle,
                                           float &outputReflectance,
                                           float & phaseAngle) {

  // We assume 1<= x and 1 <= y
  
  outputReflectance = 0.0;
  if (savePhaseAngle) phaseAngle = std::numeric_limits<double>::quiet_NaN();

  if (x <= 0 || y <= 0) return;
  
  Vector2 lonlat = DEMGeo.pixel_to_lonlat(Vector2(x, y));
  float   h      = DEMTile(x, y);
  if (h == noDEMDataValue) return;;
  Vector3 lonlat3(lonlat(0), lonlat(1), h);
  Vector3 base = DEMGeo.datum().geodetic_to_cartesian(lonlat3);

  lonlat = DEMGeo.pixel_to_lonlat(Vector2(x-1, y));
  h      = DEMTile(x-1, y);
  if (h == noDEMDataValue) return;;
  lonlat3 = Vector3(lonlat(0), lonlat(1), h);
  Vector3 x1 = DEMGeo.datum().geodetic_to_cartesian(lonlat3);

  lonlat = DEMGeo.pixel_to_lonlat(Vector2(x, y-1));
  h      = DEMTile(x, y-1);
  if (h == noDEMDataValue) return;;
  lonlat3 = Vector3(lonlat(0), lonlat(1), h);
  Vector3 y1 = DEMGeo.datum().geodetic_to_cartesian(lonlat3);

  Vector3 dx = base - x1;
  Vector3 dy = base - y1;
  Vector3 normal = -normalize(cross_prod(dx, dy));

  outputReflectance = ComputeReflectance(normal, base, input_img_params, globalParams, phaseAngle);
}

void photometry::computeReflectanceAux(ImageView<PixelGray<float> > const& DEMTile,
                                       cartography::GeoReference const& DEMGeo,
                                       float noDEMDataValue,
                                       ModelParams const& input_img_params,
                                       GlobalParams const& globalParams,
                                       ImageView<PixelMask<PixelGray<float> > >& outputReflectance,
                                       bool savePhaseAngle,
                                       ImageView<PixelMask<PixelGray<float> > >& phaseAngle) {
  
  //std::cout << "---sun        " << input_img_params.sunPosition << std::endl;
  //std::cout << "---spacecraft " << input_img_params.spacecraftPosition << std::endl;
  outputReflectance.set_size(DEMTile.cols(), DEMTile.rows());
  if (savePhaseAngle) phaseAngle.set_size(DEMTile.cols(), DEMTile.rows());

  for (int y = 1; y < (int)outputReflectance.rows(); y++) {
    for (int x = 1; x < (int)outputReflectance.cols(); x++) {

      float outputReflectanceVal, phaseAngleVal;
      computeReflectanceAtPixel(x, y, DEMTile, DEMGeo, noDEMDataValue, input_img_params, globalParams, savePhaseAngle, outputReflectanceVal, phaseAngleVal);
      
      outputReflectance(x, y) = outputReflectanceVal;
      if (savePhaseAngle) phaseAngle(x, y) = phaseAngleVal;
    }
  }
  
  return;
}

float photometry::actOnImage(std::vector<ImageRecord> & DEMTiles,
                             std::vector<ImageRecord> & albedoTiles,
                             std::vector<ImageRecord> & weightsSumTiles,
                             std::vector<int> & overlap,
                             ModelParams & input_img_params,
                             GlobalParams const& globalParams){

  // Perform one of the following operations for the given image:

  // 1. Compute the average reflectance of the current image (if
  // globalParams.initExposure is true).
  
  // 2. Update the exposure of that image (if
  // globalParams.updateExposure is true).

  // These operations use very similar logic. Namely, we will iterate
  // over all tiles overlapping with the image, and accumulate the
  // results.
  
  // We will keep in mind that the tiles partially overlap, so when we
  // accumulate the results we will ignore the padded region of the
  // tile (of size globalParams.pixelPadding).

  // Use double precision for the intermediate computations, even though
  // the input and output are float, as the float type is not precise
  // enough when accumulating a lot of numbers.
  
  std::string input_img_file = input_img_params.inputFilename;
  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  bool useWeights = (globalParams.useWeights != 0);
  if (useWeights){
    ReadWeightsParamsFromFile(input_img_params);
  }

  // Quantities needed for computing reflectance
  int count = 0;
  double reflectance_sum = 0.0;

  // Quantities needed for updating the exposure
  double numerator = 0.0, denominator = 0.0;
  
  // Iterate over the DEM tiles
  for (int i = 0; i < (int)overlap.size(); i++){

    std::string DEMTileFile = DEMTiles[overlap[i]].path;
    std::cout << "Reading file: "<< DEMTileFile << std::endl;
    DiskImageView<PixelGray<float> > DEMTile(DEMTileFile);
    GeoReference DEMGeo;
    read_georeference(DEMGeo, DEMTileFile);
    float noDEMDataValue;
    if ( !readNoDEMDataVal(DEMTileFile, noDEMDataValue)){
      std::cerr << "ERROR: Could not read the NoData Value from " << DEMTileFile << std::endl;
      exit(1);
    }

    //system("echo top1 is $(top -u $(whoami) -b -n 1|grep reconstruct)");

    // We need to keep in mind that the tile is padded with globalParams.pixelPadding pixels on each side.
    // To avoid double-counting pixels, skip them if they are part of the padding and not
    // of the tile proper.
    double min_tile_x, max_tile_x, min_tile_y, max_tile_y;
    getTileCornersWithoutPadding(// Inputs
                                 DEMTile.cols(),  DEMTile.rows(), DEMGeo,  
                                 globalParams.tileSize, globalParams.pixelPadding,  
                                 // Outputs
                                 min_tile_x, max_tile_x,  
                                 min_tile_y, max_tile_y
                                 );

    // Skip the current tile if there are no valid image pixels inside
    // of it. This optimization is necessary, but is a big of a hack
    // since it won't apply if we don't use weights (which is rare).
    if (useWeights){
      bool hasGoodPixels = false;
      Vector2 beg_pixel = input_img_geo.lonlat_to_pixel(Vector2(min_tile_x, max_tile_y));
      Vector2 end_pixel = input_img_geo.lonlat_to_pixel(Vector2(max_tile_x, min_tile_y));
      int beg_row = std::max(0,                (int)floor(beg_pixel(1)));
      int end_row = std::min(input_img.rows(), (int)ceil(end_pixel(1)));
      int beg_col = std::max(0,                (int)floor(beg_pixel(0)));
      int end_col = std::min(input_img.cols(), (int)ceil(end_pixel(0)));
      for (int k = beg_row; k < end_row; ++k) {
        if (hasGoodPixels) break;          
        for (int l = beg_col; l < end_col; ++l) {
          Vector2 input_img_pix(l,k);
          double weight = ComputeLineWeightsHV(input_img_pix, input_img_params);
          if (weight != 0) {
            hasGoodPixels = true;
            break;
          }
        }
      }
      //system("echo top2 is $(top -u $(whoami) -b -n 1|grep reconstruct)");
      if (!hasGoodPixels) continue;
    }
      
    std::string weightsSumFile = weightsSumTiles[overlap[i]].path;
    ImageView<float> weightsSum;
    if (globalParams.useNormalizedCostFun && globalParams.updateExposure){
      std::cout << "Reading " << weightsSumFile << std::endl;
      weightsSum = copy(DiskImageView<float>(weightsSumFile));
    }
    InterpolationView<EdgeExtensionView<ImageView<float>, ConstantEdgeExtension>, BilinearInterpolation>
      interp_weightsSum = interpolate(weightsSum,
                                      BilinearInterpolation(), ConstantEdgeExtension());
    
    ImageView<PixelMask<PixelGray<float> > > Reflectance;
    if (globalParams.reflectanceType != NO_REFL){
      // The reflectance type is non-zero. Compute the reflectance.
      // Declare dem_xyz and surface_normal in a block to de-allocate
      // fast the quantities dem_xyz and surface_normal which are
      // needed only temporarily. They take a lot of memory being
      // images of vectors.
      bool savePhaseAngle = false;
      ImageView<PixelMask<PixelGray<float> > > phaseAngle;
      computeReflectanceAux(DEMTile, DEMGeo, noDEMDataValue,
                            input_img_params, globalParams,  
                            Reflectance, // output
                            savePhaseAngle,
                            phaseAngle  // output  
                            );
      
      //system("echo top3 is $(top -u $(whoami) -b -n 1|grep reconstruct)");
    }else{
      // The reflectance is set to 1.
      Reflectance.set_size(DEMTile.cols(), DEMTile.rows());
      for (int row = 0; row < (int)DEMTile.rows(); row++){
        for (int col = 0; col < (int)DEMTile.cols(); col++){
          Reflectance(col, row) = 1.0;
          Reflectance(col, row).validate();
        }
      }
    }
    
    InterpolationView<EdgeExtensionView<ImageView< PixelMask<PixelGray<float> > >,
      ConstantEdgeExtension>, BilinearInterpolation>
      interp_reflectance = interpolate(Reflectance,
                                       BilinearInterpolation(), ConstantEdgeExtension());
    
    
    ImageView<PixelMask<PixelGray<float> > > albedoTile;
    if (globalParams.updateExposure){
      // Update the exposure. Need the current albedo.
      std::string albedoTileFile = albedoTiles[overlap[i]].path;
      std::ifstream aFile(albedoTileFile.c_str());
      if (!aFile) continue;
      std::cout << "Reading " << albedoTileFile << std::endl;
      albedoTile = copy(DiskImageView<PixelMask<PixelGray<uint8> > >(albedoTileFile));

      // Sanity check: if the user is not careful, the padding may have been stripped
      // from the albedo tile by now, which of course would cause a size mis-match.
      if (DEMTile.rows() != albedoTile.rows() || DEMTile.cols() != albedoTile.cols()){
        std::cout << "ERROR: We expect the DEM and albedo tiles to have the same number "
                  << "of rows/columns." << std::endl;
        exit(1);
      }
    }
    
    InterpolationView<EdgeExtensionView<ImageView< PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_albedo = interpolate(albedoTile,
                                  BilinearInterpolation(), ConstantEdgeExtension());

    //system("echo top4 is $(top -u $(whoami) -b -n 1|grep reconstruct)");

    // Iterate only over the portion of the image intersecting the DEM tile
    Vector2 beg_pixel = input_img_geo.lonlat_to_pixel(Vector2(min_tile_x, max_tile_y));
    Vector2 end_pixel = input_img_geo.lonlat_to_pixel(Vector2(max_tile_x, min_tile_y));
    int beg_row = std::max(0,                (int)floor(beg_pixel(1)));
    int end_row = std::min(input_img.rows(), (int)ceil(end_pixel(1)));
    int beg_col = std::max(0,                (int)floor(beg_pixel(0)));
    int end_col = std::min(input_img.cols(), (int)ceil(end_pixel(0)));
    for (int k = beg_row; k < end_row; ++k) {
      for (int l = beg_col; l < end_col; ++l) {

        Vector2 input_img_pix(l,k);

        Vector2 lonlat = input_img_geo.pixel_to_lonlat(input_img_pix);
        double lx = lonlat(0), ly = lonlat(1);
        bool isInTile = (min_tile_x <= lx && lx < max_tile_x && min_tile_y <= ly && ly < max_tile_y);
        if (!isInTile) continue;

        //if (k%2 != 0 && l%2 != 0) continue; // Skip some points when computing the exposure to save time
          
        // Check that all four grid points used for interpolation are valid
        Vector2 reflectance_pix = DEMGeo.lonlat_to_pixel(lonlat);
        double x = reflectance_pix[0];
        double y = reflectance_pix[1];
        bool isValidReflectance = (x>=0) && (x <= Reflectance.cols()-1) && (y>=0) && (y<= Reflectance.rows()-1) &&
          is_valid(Reflectance( floor(x), floor(y))) && ((double)Reflectance( floor(x), floor(y) ) != 0) &&
          is_valid(Reflectance( floor(x), ceil (y))) && ((double)Reflectance( floor(x), ceil(y)  ) != 0) &&
          is_valid(Reflectance( ceil (x), floor(y))) && ((double)Reflectance( ceil (x), floor(y) ) != 0) &&
          is_valid(Reflectance( ceil (x), ceil (y))) && ((double)Reflectance( ceil (x), ceil(y)  ) != 0);
        if (!isValidReflectance) continue;
        
        double R = interp_reflectance(x, y);
        double exposureRefl = R*input_img_params.exposureTime;
        
        // Check if the image is above the shadow threshold
        double t = getShadowThresh(globalParams, exposureRefl);
        bool isBlack = !( is_valid(input_img(l, k)) && (double)input_img(l, k) >= t );
        if (isBlack) continue;

        if (globalParams.initExposure){
          reflectance_sum += R;
          count++;
        }else if (globalParams.updateExposure){
          // Update the exposure
          // We assume that the DEM, reflectance, and albedo are on the same grid

          // Find the weight, and normalize it so that all weights at the current pixel sum up to 1.
          double weight = 1.0;
          if (useWeights) weight = ComputeLineWeightsHV(input_img_pix, input_img_params);
          if (globalParams.useNormalizedCostFun){
            double weights_sum = interp_weightsSum(x, y);
            if (weights_sum == 0) continue;
            weight /= weights_sum;
          }
                
          double A  = interp_albedo(x, y); // we assume albedo, reflectance, and DEM are on the same grid
          double RA = R*A;
          numerator   += ((double)input_img(l, k) - input_img_params.exposureTime*RA)*RA*weight;
          denominator += RA*RA*weight;
        }
      }
          
    }
    //system("echo top5 is $(top -u $(whoami) -b -n 1|grep reconstruct)");
  }

  if (globalParams.initExposure){
    // Compute the average reflectance
    double avg_reflectance = reflectance_sum/count;
    printf("avg_reflectance = %f\n", avg_reflectance);
    return avg_reflectance;
  }else if (globalParams.updateExposure){
    // Update the exposure
    double ratio = 0.0;
    if (denominator != 0) ratio = numerator/denominator;
    double exposure = input_img_params.exposureTime + ratio;
    printf("updated exposure is = %f\n", exposure);
    return exposure;
  }

  return 0;
}

// --------- Only old code below -------------------

//computes a reflectance image
//author: Ara Nefian
float photometry::computeImageReflectanceOld(ModelParams const& input_img_params,
                                             GlobalParams const& globalParams) {
  std::string input_img_file = input_img_params.inputFilename;
  std::string output_img_file = input_img_params.reliefFilename;;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  ImageView<PixelMask<PixelGray<float> > > output_img;
  float avg_reflectance =
    computeImageReflectanceNoWriteOld(input_img_params,
                                      globalParams,
                                      output_img);


  std::cout << "Writing " << output_img_file << std::endl;
  write_georeferenced_image
    (output_img_file,
     output_img,
     input_img_geo,
     TerminalProgressCallback("{Core}","Processing:"));
  
  return avg_reflectance;
}


//computes a reflectance image
//author: Ara Nefian
float photometry::computeImageReflectanceOld(ModelParams const& input_img_params,
                                             ModelParams const& overlap_img_params,
                                             GlobalParams const& globalParams) {
  int l, k;
  int count = 0;
  float avg_reflectance = 0.0;
  float overlap_avg_reflectance = 0.0;

  std::string input_img_file = input_img_params.inputFilename;
  std::string overlap_img_file = overlap_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;
  std::string overlap_shadow_file = overlap_img_params.shadowFilename;
  std::string output_img_file = input_img_params.reliefFilename;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

  ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());

  Vector3 xyz;
  Vector3 xyz_prior;
  int x, y;

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation>
    interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
                                               ConstantEdgeExtension()),
                                   BilinearInterpolation());



  // start
  DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_file);
  GeoReference overlap_geo;
  read_georeference(overlap_geo, overlap_img_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_shadow_file);
  
  ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                            ConstantEdgeExtension()),
                                                                                BilinearInterpolation());

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
                                                                                                  ConstantEdgeExtension()),
                                                                                      BilinearInterpolation());
  //end

  //initialize the nominator and denomitor images
  for (k = 0 ; k < (int)input_img.rows(); ++k) {
    for (l = 0; l < (int)input_img.cols(); ++l) {

      Vector2 input_image_pix(l,k);

      if ( is_valid(input_img(l,k)) ) {

        //get the corresponding DEM value

        Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
        Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));

        x = (int)input_dem_pix[0];
        y = (int)input_dem_pix[1];

        //check for valid DEM coordinates to be within the DEM boundaries
        if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){

          Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
          Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);

          Vector2 input_img_left_pix;
          input_img_left_pix(0) = l-1;
          input_img_left_pix(1) = k;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = l;
          input_img_top_pix(1) = k-1;

          //check for valid DEM pixel value and valid left and top coordinates
          if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != globalParams.noDEMDataValue)){

            //determine the 3D coordinates of the pixel left of the current pixel
            Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
            Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
            Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
            Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);

            //determine the 3D coordinates of the pixel top of the current pixel
            Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
            Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
            Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
            Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);

            Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);
            //This part is the only image depedent part - START

            float input_img_reflectance, phaseAngle;
            input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams, phaseAngle);

            output_img(l, k) = input_img_reflectance;

            //--------------------------------------------------------------
            //check for overlap between the output image and the input DEM image
            Vector2 input_img_pix(l,k);
            Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
            x = (int)overlap_pix[0];
            y = (int)overlap_pix[1];

            //image dependent part of the code  - START
            PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

            //check for valid overlap_img coordinates
            //remove shadow pixels in the overlap_img.
            if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && (interpOverlapShadowImage(x, y) == 0) && (shadowImage(l,k) == 0)) {

              if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                float overlap_img_reflectance, phaseAngle;
                overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params, globalParams, phaseAngle);
                if ((overlap_img_reflectance > 0) && (input_img_reflectance > 0)){
                  avg_reflectance = avg_reflectance + output_img(l,k);
                  overlap_avg_reflectance = overlap_avg_reflectance + overlap_img_reflectance;
                  count++;
                }
              }//if
            }//if
            //----------------------------
          }
        }
      }
    }
  }

  avg_reflectance = avg_reflectance/count;
  overlap_avg_reflectance = overlap_avg_reflectance/count;
  float reflectance_ratio = avg_reflectance/overlap_avg_reflectance;
  printf("avg_reflectance = %f, overlap_avg_reflectance = %f\n", avg_reflectance, overlap_avg_reflectance);
  write_georeferenced_image(output_img_file,
                            output_img,
                            input_img_geo, TerminalProgressCallback("{Core}","Processing:"));

  return reflectance_ratio;
}

//generates the normal of a point p1 from the 3D coordinates of p1, p2, p3
//this function is currently buggy. Needs to be fixed and this will replace
//computeNormalFrom3DPoints
Vector3
photometry::computeNormalFrom3DPointsGeneral(Vector3 p1,
                                             Vector3 p2,
                                             Vector3 p3) {
  return -normalize(cross_prod(p2-p1,p3-p1));
}

//generates the normal of a point p1 from the 3D coordinates of p1, p2, p3
//assumes unit dstance between adjacent pixels in x, y coordinates
//this function will be grandfathered
//author Ara Nefian
Vector3
photometry::computeNormalFrom3DPoints(Vector3 p1, Vector3 p2,
                                      Vector3 p3) {
  Vector3 normal;

  normal[0] = p1[2]-p2[2];
  normal[1] = p1[2]-p3[2];
  normal[2] = -1;

  return normalize(normal);
}

//computes the Lunar-Lambertian reflectance mode
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float
photometry::computeLunarLambertianReflectanceFromNormalOld(Vector3 sunPos, Vector3 viewPos, Vector3 xyz,  Vector3 normal, float B_0, float L) {
  float reflectance;

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  // sun coordinates relative to the xyz point on the Moon surface
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  // viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  reflectance = B_0*((2*L*mu_0/(mu_0+mu)) + (1-L)*mu_0 );
  return reflectance;
}

void photometry::computeXYZandSurfaceNormalOld(ImageView<PixelGray<float> > const& DEMTile,
                                               cartography::GeoReference const& DEMGeo,
                                               float noDEMDataValue,
                                               ImageView<Vector3> & dem_xyz,
                                               ImageView<Vector3> & surface_normal){


  // convert dem altitude to xyz cartesian pixels
  dem_xyz = geodetic_to_cartesian( dem_to_geodetic( DEMTile, DEMGeo ),
                                   DEMGeo.datum() );

  // transfer noDEMDataValue values
  for (int y=0; y < (int)dem_xyz.rows(); y++) {
    for (int x=0; x < (int)dem_xyz.cols(); x++) {
      if (DEMTile(x, y) == noDEMDataValue) {
        dem_xyz(x, y) = Vector3();
      }
    }
  }  

  // Init the surface normal to zero
  surface_normal.set_size(dem_xyz.cols(), dem_xyz.rows());
  for (int y=0; y < (int)surface_normal.rows(); y++) {
    for (int x=0; x < (int)surface_normal.cols(); x++) {
      surface_normal(x, y) = Vector3();
    }
  }

  for (int y=1; y < (int)surface_normal.rows(); y++) {
    for (int x=1; x < (int)surface_normal.cols(); x++) {
      
      Vector3 base = dem_xyz(x, y);
      if (base == Vector3()) {
        continue;
      }
      Vector3 x1 = dem_xyz(x-1, y);
      if (x1 == Vector3()) {
        continue;
      }
      Vector3 y1 = dem_xyz(x, y-1);
      if (y1 == Vector3()) {
        continue;
      }

      Vector3 dx = base - x1;
      Vector3 dy = base - y1;
      surface_normal(x, y) = -normalize(cross_prod(dx, dy));
    }
  }

  return;
}

void photometry::computeReflectanceAuxOld(ImageView<Vector3> const& dem_xyz,
                                          ImageView<Vector3> const& surface_normal,
                                          ModelParams const& input_img_params,
                                          GlobalParams const& globalParams,
                                          ImageView<PixelMask<PixelGray<float> > >& outputReflectance,
                                          bool savePhaseAngle,
                                          ImageView<PixelMask<PixelGray<float> > >& phaseAngle) {
  
  //std::cout << "---sun        " << input_img_params.sunPosition << std::endl;
  //std::cout << "---spacecraft " << input_img_params.spacecraftPosition << std::endl;
  outputReflectance.set_size(surface_normal.cols(), surface_normal.rows());
  if (savePhaseAngle) phaseAngle.set_size(surface_normal.cols(), surface_normal.rows());
  for (int y=0; y < (int)outputReflectance.rows(); y++) {
    for (int x=0; x < (int)outputReflectance.cols(); x++) {
      Vector3 normal = surface_normal(x, y);
      if (normal == Vector3()) {
        outputReflectance(x, y) = 0;
        outputReflectance(x, y).invalidate();
        if (savePhaseAngle) phaseAngle(x, y) = std::numeric_limits<double>::quiet_NaN();
      }else{
        Vector3 xyz = dem_xyz(x, y);
        float phaseAngleVal;
        //outputReflectance(x, y) = (ComputeReflectance(normal, xyz, input_img_params, globalParams, phaseAngleVal) != 0); // temporary!!!
        outputReflectance(x, y) = ComputeReflectance(normal, xyz, input_img_params, globalParams, phaseAngleVal);
        if (savePhaseAngle) phaseAngle(x, y) = phaseAngleVal;
      }
    }
  }
  
  return;
}

float photometry::computeImageReflectanceNoWriteOld(ModelParams const& input_img_params,
                                                    GlobalParams const& globalParams,
                                                    ImageView<PixelMask<PixelGray<float> > >& output_img) {
  
  std::string input_img_file = input_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  // warp dem to drg georef
  ImageView<PixelGray<double> > dem_in_drg_georef_0(input_img.cols(), input_img.rows());
  DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);
  dem_in_drg_georef_0 =
    crop
    (geo_transform
     (input_dem_image,
      input_dem_geo,
      input_img_geo,
      ConstantEdgeExtension(),
      BilinearInterpolation()),
     bounding_box(input_img)
     );
  
#if 0
  std::cout << std::endl;
  std::string tmp  = "./" + prefix_from_filename(suffix_from_filename(input_img_file)) + "_dem.tif";
  std::cout << "Writing " << tmp << std::endl;
  write_georeferenced_image(tmp,
                            dem_in_drg_georef_0,
                            input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
#endif
    
  // grow no-data areas to avoid bogus values created by interpolating
  // data with no-data
  ImageView<PixelGray<double> > dem_in_drg_georef = copy(dem_in_drg_georef_0);
  int nodata = globalParams.noDEMDataValue;
  for (int y=1; y < (int)dem_in_drg_georef.rows()-1; y++) {
    for (int x=1; x < (int)dem_in_drg_georef.cols()-1; x++) {

#define CHECK_PIX(dx,dy)                                \
      if (dem_in_drg_georef_0(x+dx, y+dy) == nodata) {  \
        dem_in_drg_georef(x, y) = nodata;               \
        continue;                                       \
      }

      CHECK_PIX(-1, -1);
      CHECK_PIX(-1,  0);
      CHECK_PIX(-1,  1);
      CHECK_PIX( 0, -1);
      CHECK_PIX( 0,  1);
      CHECK_PIX( 1, -1);
      CHECK_PIX( 1,  0);
      CHECK_PIX( 1,  1);
    }
  }

  ImageView<Vector3> dem_xyz;
  ImageView<Vector3> surface_normal;
  
  bool savePhaseAngle = false;
  ImageView<PixelMask<PixelGray<float> > > phaseAngle;
  computeXYZandSurfaceNormalOld(dem_in_drg_georef, input_img_geo, globalParams.noDEMDataValue, // Inputs 
                                dem_xyz, surface_normal                                        // Outputs
                                );
  computeReflectanceAuxOld(dem_xyz, surface_normal,  
                           input_img_params,
                           globalParams,  
                           output_img,
                           savePhaseAngle, phaseAngle
                           );
  
  // compute average reflectance
  int count = 0;
  float reflectance_sum = 0.0;
  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);
  for (int y=0; y < (int)output_img.rows(); y++) {
    for (int x=0; x < (int)output_img.cols(); x++) {
      float refl = output_img(x, y);
      if ((refl != 0.0) && (shadowImage(x, y) == 0)){ //valid non zero reflectance
        reflectance_sum += refl;
        count++;
      }
    }
  }
  float avg_reflectance = reflectance_sum / count;

  printf("avg_reflectance = %f\n", avg_reflectance);
  return avg_reflectance;
}

