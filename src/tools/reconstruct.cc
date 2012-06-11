// __BEGIN_LICENSE__
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

// To do:
// Good testcase for cutting: AS17-M-1982
// Add documentation about generated kmls
// Utility for  extracting sun and spacecraft position from a cube. Maybe Ara has it?
// Put info on the uint testcases in the documentation.
// Put the Apache license everywhere!
// Make the unit tests work
// Wipe all the code having to do with images.
// How to create kml tiles?
// Massive cleanup needed in reconstruct.cc. Remove all flow without tiles.
// Integrate all pieces calling actOnTile. Remove the loop over images.
// Test that in mosaic only mode the code can run without the DEM directory
// Fix a bug at 10 mpp. Zoom in a lot at 84E.
// When updating the exposure, some processes write exposure while
// some other processes at the same time read it. This can cause problems.
// Use normalized weights in shape-from-shading as well.
// Wipe the network machines code. Create a copy of reconstruct.sh for use
// on the supercomputer.
// Add tool to StereoPipeline to extract the sun and spacecraft
// position from cube. Then wipe reconstruct_aux.cc and everything
// having to do with ISIS from reconstruct.cc and reconstruct.sh.
// Add a blurb in the documentation about how to get images and
// sun/spacecraft position from isis cubes.
// Fix the bug in orthoproject with images going beyond 180 degrees. They show up
// with a huge black band.
// Simplify the script.
// More work in shape from shading: Strip padding at the last iteration.
// Bug: There is a shift in the AS17 mosaic.
// Validate the list of machines, the number of CPUs
// Validate the  sim box
// Check if  reconstruct, orthoproject, image2qtree, gdal, exist.
// To do: How to find the number of processors on a Mac
// Validate in advance if we have sun/spacecraft position for each image.
// Should images not having sun/spacecraft info be ignored?
// To do: Must regenerate the index files all the time, as they are too fragile
// Bad image: DIM_input_10mpp/AS17-M-0473.tif
// To do: Fix the logic for extracting sun/spacecraft position from cubes
// in both the .cc file and in the shell script.
// To do: Enable the flow of doing all the flow w/o the shell script.
// To do: Remove the file reconstruct_aux.cc which mostly duplicates orthoproject.cc.
// To do: Test this big image: AS17-M-0305
// To do: Test this as well. DIM_input_sub64_isis/AS17-M-0281.tif, has a lot of black.
// Note: We assume a certain convention about the isis cube file and corresponding
// isis adjust file.
// Note: The isis adjust file may not exist.
// To do: Move readDEMTilesIntersectingBox to Image/ and remove the DEM-specific language.
// To do: Note that the image extracted from cubes is uint16, not uint8.
// To do: Also note that those images are a bit darker than regular
// images we already have, even after converting them to uint8.
// To: The function ComputeGeoBoundary looks wrong. Use instead
// georef.bounding_box(img), which I think is accurate. Compare
// with gdalinfo.
// To do: Check the effect of pixel padding on final albedo. 
// To do: Implement the logic where one checks if a pixel is in the shadow,
// in one (inline) function, and call it wherever it is needed.
// To do: Copy the images from supercomp.
// To do: Fix the memory leaks where the weighs are read.
// Read the data from cubes, and remove all logic having to do with
// filters from reconstruct.sh.
// To do: Unit tests
// To do: Reorg the code which computes the reflectance and its
// derivative in ShapeFromShading.cc to only compute the
// derivative. Move all that code to Reflectance.cc. Convert all to
// double.
// Rename dem_out.tif to dem_mean.tif, and modelParams.outputFile to modelParams.albedoFile,
// inputFile to drgFile.
// Remove a lot of duplicate code related to overlaps
// There is only one image, rm the vector of images
// No need to initialize the albedo tiles on disk, create
//   them on the fly.
// Merge the imageRecord and modelParams classes
// See if to change the order of values in the corners vector
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <photk/Photometry.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace photometry;


int main( int argc, char *argv[] ) {

  for (int s = 0; s < argc; s++) std::cout << argv[s] << " ";
  std::cout << std::endl;

  std::vector<std::string> inputDRGFiles;
  std::string resDir, settingsFile, DRGInBoxList, albedoTilesList;
  std::string currImgOrTile;
  
  po::options_description general_options("Options");
  general_options.add_options()
    ("settings-file,s", po::value<std::string>(&settingsFile), "Settings file")
    ("results-directory,r", po::value<std::string>(&resDir),   "Results directory")
    ("images-list,f", po::value<std::string>(&DRGInBoxList),   "The list of images")
    ("tiles-list,t", po::value<std::string>(&albedoTilesList), "The list of albedo tiles")
    ("image-file,i", po::value<std::string>(&currImgOrTile),   "Current image or tile")
    ("initial-setup",            "Initial setup")
    ("save-weights",             "Save the weights")
    ("compute-weights-sum",      "Compute the sum of weights at each pixel")
    ("compute-shadow",           "Compute the shadow")
    ("init-dem",                 "Initialize the DEM")
    ("init-exposure",            "Initialize the exposure times")
    ("init-albedo",              "Initialize the albedo")
    ("update-exposure",          "Update the exposure times")
    ("update-tile-phase-coeffs", "Update the phase coefficients per tile")
    ("update-phase-coeffs",      "Update the phase coefficients by combining the results over all tiles")
    ("update-albedo",            "Update the albedo")
    ("update-height",            "Update the height (shape from shading)")
    ("compute-errors",           "Compute the errors in albedo")
    ("is-last-iter",             "Is this the last iteration")
    ("help,h",                   "Display this help message");

  po::options_description hidden_options("");

  hidden_options.add_options()
    ("inputDRGFiles", po::value<std::vector<std::string> >(&inputDRGFiles));
 
  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("inputDRGFiles", -1);

  std::ostringstream usage;
  usage << "Description: main code for albedo, mosaic, and shape reconstruction from multiple images" << std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch ( po::error const& e ) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if  (vm.count("help") || argc <= 1) {
    std::cerr << usage.str() << std::endl;
    return 1;
  }

  // Read the global parameters settings. Apply the command-line overrides.
  GlobalParams globalParams;
  ReadSettingsFile((char*)settingsFile.c_str(), globalParams);
  if ( vm.count("initial-setup"            ) ) globalParams.initialSetup          = true;
  if ( vm.count("save-weights"             ) ) globalParams.saveWeights           = true;
  if ( vm.count("compute-weights-sum"      ) ) globalParams.computeWeightsSum     = true;
  if ( vm.count("init-dem"                 ) ) globalParams.initDEM               = true;
  if ( vm.count("init-exposure"            ) ) globalParams.initExposure          = true;
  if ( vm.count("init-albedo"              ) ) globalParams.initAlbedo            = true;
  if ( vm.count("update-exposure"          ) ) globalParams.updateExposure        = true;
  if ( vm.count("update-tile-phase-coeffs" ) ) globalParams.updateTilePhaseCoeffs = true;
  if ( vm.count("update-phase-coeffs"      ) ) globalParams.updatePhaseCoeffs     = true;
  if ( vm.count("update-albedo"            ) ) globalParams.updateAlbedo          = true;
  PrintGlobalParams(globalParams);

  bool updateHeight  = vm.count("update-height");
  bool computeErrors = vm.count("compute-errors");
  bool isLastIter    = vm.count("is-last-iter");

  // Validation
  if ((int)globalParams.initialSetup
      + (int)globalParams.saveWeights
      + (int)globalParams.computeWeightsSum
      + (int)globalParams.initDEM
      + (int)globalParams.initExposure
      + (int)globalParams.initAlbedo
      + (int)globalParams.updateExposure
      + (int)globalParams.updateAlbedo
      + (int)globalParams.updateTilePhaseCoeffs
      + (int)globalParams.updatePhaseCoeffs
      + (int)updateHeight
      + (int)computeErrors > 1 ){
    std::cerr << "ERROR: Two or more actions requested which cannot be processed at the same time." << std::endl;
    exit(1);
  }
  
  // Double check to make sure all folders exist  
  if ( !fs::exists(resDir) )
    fs::create_directories(resDir);
  if ( !fs::exists(resDir+"/info") )
    fs::create_directories(resDir+"/info");
  if ( !fs::exists(resDir+"/reflectance") )
    fs::create_directories(resDir+"/reflectance");
  if ( !fs::exists(resDir+"/shadow") )
    fs::create_directories(resDir+"/shadow");
  if ( !fs::exists(resDir+"/exposure") )
    fs::create_directories(resDir+"/exposure");
  if ( !fs::exists(resDir+"/weight") )
    fs::create_directories(resDir+"/weight");
  std::string albedoDir     = resDir + "/albedo";
  std::string meanDEMDir    = resDir + "/DEM";
  std::string costFunDir    = resDir + "/costFun";
  std::string errorDir      = resDir + "/error";
  std::string weightsSumDir = resDir + "/weightsSum";
  std::string sfsDir        = resDir + "/DEM_sfs";
  std::string phaseDir      = resDir + "/phase";
  if ( !fs::exists(albedoDir    ) ) fs::create_directories(albedoDir    );
  if ( !fs::exists(meanDEMDir   ) ) fs::create_directories(meanDEMDir   );
  if ( !fs::exists(costFunDir   ) ) fs::create_directories(costFunDir   );
  if ( !fs::exists(errorDir     ) ) fs::create_directories(errorDir     );
  if ( !fs::exists(weightsSumDir) ) fs::create_directories(weightsSumDir);
  if ( !fs::exists(sfsDir       ) ) fs::create_directories(sfsDir       );
  if ( !fs::exists(phaseDir     ) ) fs::create_directories(phaseDir     );

  // This will overwrite the values of globalParams.phaseCoeffA1, globalParams.phaseCoeffA2
  // if that information is available on disk.
  ReadPhaseCoeffsFromFile(phaseDir, globalParams);
  if (globalParams.initialSetup == 1) AppendPhaseCoeffsToFile(globalParams);

  // The names of the files listing all DRGs and DEMs and the coordinates
  // of their corners. 
  std::string allDRGIndex  = globalParams.drgDir + "/index.txt";
  std::string allDEMIndex  = globalParams.demDir + "/index.txt";
  std::string DEMTilesList = resDir + "/DEMTilesList.txt";

  // Check if the input image exists
  std::ifstream iFile(currImgOrTile.c_str());
  if (!iFile) return 0;

  // Not using reflectance is the same as creating a mosaic instead of albedo
  bool useReflectance = (globalParams.reflectanceType != NO_REFL); 

  if (globalParams.initialSetup == 1) {

    // This block of code creates the list of images (DRGInBoxList). As such, it must be above
    // any code which reads the list of images.
    // Create the DRGInBoxList used in subsequent iterations.
    // Create the list of all DEM if not there yet.
    listDRGinBoxAndAllDEM(useReflectance,
                          allDRGIndex, allDEMIndex,
                          globalParams.simulationBox, globalParams.drgDir, globalParams.demDir, DRGInBoxList
                          );
    
    vw_out( VerboseDebugMessage, "photometry" ) << "Initializing the albedo tiles ... ";
    std::vector<ImageRecord> drgRecords;
    if (!readImagesFile(drgRecords, DRGInBoxList)) exit(1);
    std::string imageFile = currImgOrTile; // an image whose georef we will use
    createAlbedoTilesOverlappingWithDRG(globalParams.tileSize, globalParams.pixelPadding,
                                        imageFile, globalParams.simulationBox,
                                        drgRecords,
                                        DEMTilesList,    meanDEMDir,
                                        albedoTilesList, albedoDir
                                        );
  }

  std::vector<ImageRecord> drgRecords;
  if (!readImagesFile(drgRecords, DRGInBoxList)) exit(1);

  std::map<std::string, Vector3> sunPositions;
  std::map<std::string, Vector3> spacecraftPositions;
  if (useReflectance){
    ReadSunOrSpacecraftPosition(globalParams.sunPosFile, // Input
                                sunPositions             // Output
                                );
    ReadSunOrSpacecraftPosition(globalParams.spacecraftPosFile, // Input
                                spacecraftPositions             // Output
                                );
  }
  
  std::vector<ModelParams> modelParamsArray;

  //this will contain all the DRG files
  std::vector<std::string> DRGFiles;
  {
    int i = 0;
    for (unsigned int j = 0; j < drgRecords.size(); ++j) {
      if (!drgRecords[j].useImage) continue;
      
      DRGFiles.push_back(drgRecords[j].path);
      modelParamsArray.push_back(ModelParams());
      
      std::string temp = suffix_from_filename(DRGFiles[i]);
      modelParamsArray[i].exposureTime     = 1.0;
      modelParamsArray[i].hCenterLineDEM   = NULL;
      modelParamsArray[i].hMaxDistArrayDEM = NULL;
      modelParamsArray[i].inputFilename    = DRGFiles[i];//these filenames have full path
      
      std::string prefix = getFirstElevenCharsFromFileName(DRGFiles[i]);
      
      if (useReflectance){
        if ( sunPositions.find(prefix) == sunPositions.end()){
          std::cerr << "Could not find the sun position for the DRG file: " << DRGFiles[i] << std::endl;
          exit(1);
        }
        modelParamsArray[i].sunPosition = 1000*sunPositions[prefix];
        
        if (spacecraftPositions.find(prefix) == spacecraftPositions.end()){
          std::cerr << "Could not find the spacecraft position for the DRG file: " << DRGFiles[i] << std::endl;
          exit(1);
        }
        modelParamsArray[i].spacecraftPosition = 1000*spacecraftPositions[prefix];
      }
      
      modelParamsArray[i].infoFilename        = resDir + "/info/" + prefix_less3_from_filename(temp)+"info.txt";
      modelParamsArray[i].reliefFilename      = resDir + "/reflectance" + prefix_from_filename(temp) + "_reflectance.tif";
      modelParamsArray[i].shadowFilename      = resDir + "/shadow/" + prefix_from_filename(temp) + "_shadow.tif";
      modelParamsArray[i].errorFilename       = resDir + "/error" + prefix_from_filename(temp) + "_err.tif";
      modelParamsArray[i].outputFilename      = resDir + "/albedo" + prefix_from_filename(temp) + "_albedo.tif";
      modelParamsArray[i].sfsDEMFilename      = resDir + "/DEM_sfs" + prefix_less3_from_filename(temp) + "DEM_sfs.tif";
      modelParamsArray[i].errorHeightFilename = resDir + "/error" + prefix_from_filename(temp) + "_height_err.tif";
      modelParamsArray[i].weightFilename      = resDir + "/weight" + prefix_from_filename(temp) + "_weight.txt";
      modelParamsArray[i].exposureFilename    = resDir + "/exposure" + prefix_from_filename(temp) + "_exposure.txt";
      
      const ImageRecord& rec = drgRecords[j];
      modelParamsArray[i].corners = Vector4(rec.west, rec.east, rec.south, rec.north);
      
      modelParamsArray[i].hCenterLineDEM   = NULL;
      modelParamsArray[i].hMaxDistArrayDEM = NULL;
      modelParamsArray[i].vCenterLineDEM   = NULL;
      modelParamsArray[i].vMaxDistArrayDEM = NULL;
      
      vw_out( VerboseDebugMessage, "photometry" ) << modelParamsArray[i] << "\n";
      
      i++;
    }
  }
  
  if (globalParams.initialSetup == 1){
    // We performed all tasks, including validation, if we are in the initial setup.
    return 0;
  }
    
  vw_out() << "Number of Files = " << DRGFiles.size() << "\n";
  
  std::vector<int>  overlapIndicesArray = makeOverlapList(modelParamsArray, currImgOrTile);
  printOverlapList(overlapIndicesArray);

  std::vector<int> inputIndices = GetInputIndices(currImgOrTile, DRGFiles);
  
  // set up weights only for images that we want to process or that
  // overlap one of the files we want to process
  std::vector<int> relevantIndices;
  for (unsigned int i=0; i < inputIndices.size(); i++) {
    relevantIndices.push_back(inputIndices[i]);
  }
  for (unsigned int j=0; j < overlapIndicesArray.size(); j++) {
    relevantIndices.push_back(overlapIndicesArray[j]);
  }

  // Read the exposure information for the relevant indices
  for (unsigned int i=0; i < relevantIndices.size(); i++) {
    int j = relevantIndices[i];
    ReadExposureInfoFromFile(&(modelParamsArray[j]));
  }
  
  if (globalParams.useWeights == 1) {
    vw_out( VerboseDebugMessage, "photometry" ) << "Computing weights ... ";
    
    for (unsigned int i=0; i < relevantIndices.size(); i++) {
      int j = relevantIndices[i];

      // If we are not in weight saving mode, the weights should be on
      // disk already, so read them.
      if (globalParams.saveWeights != 1) continue;
      
      // We have globalParams.saveWeights == 1. Build the weights.
      if ((int)inputIndices.size() == 0){
        cerr << "Error: Could not find the image to process: " << currImgOrTile << " in the list of input images."
             << endl;
        exit(1);
      }
      
      if (j == inputIndices[0]){
        // Compute and save the weights only for the current image,
        // not for all images overlapping with it.

        ComputeImageCenterLines(modelParamsArray[j]);

        if (globalParams.saveWeights == 1) SaveWeightsParamsToFile(modelParamsArray[j]);
        
      }
      
    }
    
    vw_out( VerboseDebugMessage, "photometry" ) << "Done.\n";
  }
  
  if ( globalParams.initDEM == 1 && useReflectance ){
    // Initialize the DEM tiles
    std::string albedoTileFile = currImgOrTile;
    std::string DEMTileFile    = meanDEMDir + suffix_from_filename(albedoTileFile);
    std::vector<ImageRecord> DEMImages;
    if (!readImagesFile(DEMImages, allDEMIndex)) exit(1);
    std::vector<int> overlap = makeOverlapList(DEMImages, albedoTileFile);
    InitMeanDEMTile(albedoTileFile, DEMTileFile, 
                    DEMImages, overlap, globalParams
                    );
  }
  
  if ( useReflectance &&
       ( globalParams.initExposure == 1 || globalParams.updateExposure == 1 )
       ){

    // Init/update reflectance

    if ((int)inputIndices.size() == 0){
        cerr << "Error: Could not find the image to process: " << currImgOrTile << " in the list of input images."
             << endl;
        exit(1);
    }
    
    std::string currDRG = currImgOrTile;
    std::vector<ImageRecord> DEMTiles, albedoTiles, weightsSumTiles;
    std::vector<int> overlap;
    if (!readImagesFile(DEMTiles,    DEMTilesList))    exit(1);
    if (!readImagesFile(albedoTiles, albedoTilesList)) exit(1);
    weightsSumTiles = albedoTiles;
    for (int s = 0; s < (int)weightsSumTiles.size(); s++){
      // From the list of albedo tiles get the list of weights sum tiles by simply
      // substituting the right directory name.
      weightsSumTiles[s].path = weightsSumDir + suffix_from_filename(weightsSumTiles[s].path);
    }
    overlap = makeOverlapList(DEMTiles, currDRG);
    float val = actOnImage(DEMTiles, albedoTiles, weightsSumTiles,
                           overlap,
                           modelParamsArray[inputIndices[0]],
                           globalParams);
    
    if (globalParams.initExposure){
      // exposure time = TRConst/avgReflectance
      modelParamsArray[inputIndices[0]].exposureTime = globalParams.TRConst/val;
    }else if (globalParams.updateExposure){
      modelParamsArray[inputIndices[0]].exposureTime = val;
    }
    
    AppendExposureInfoToFile(modelParamsArray[inputIndices[0]]);
  }
  
  if ( useReflectance && globalParams.updatePhaseCoeffs){
    // Combine the components of the phase coefficients over all tiles
    double A1_num = 0.0, A1_den = 0.0, A2_num = 0.0, A2_den = 0.0;
    std::vector<ImageRecord> albedoTiles;
    if (!readImagesFile(albedoTiles, albedoTilesList)) exit(1);
    for (int s = 0; s < (int)albedoTiles.size(); s++){
      // From the list of albedo tiles get the list of phase coeffs per tile by simply
      // substituting the right directory name.
      //std::string phaseTileFile = phaseDir + suffix_from_filename(albedoTiles[s].path);
      std::string phaseTileFile = phaseDir + prefix_from_filename(suffix_from_filename(albedoTiles[s].path)) + ".txt";
      phaseCoeffsData PCD;
      PCD.readFromFile(phaseTileFile);
      A1_num += PCD.phaseCoeffA1_num; A1_den += PCD.phaseCoeffA1_den;
      A2_num += PCD.phaseCoeffA2_num; A2_den += PCD.phaseCoeffA2_den;
    }
    if (A1_den != 0.0) globalParams.phaseCoeffA1 += A1_num/A1_den;
    if (A2_den != 0.0) globalParams.phaseCoeffA2 += A2_num/A2_den;
    AppendPhaseCoeffsToFile(globalParams);
  }
  
  if (globalParams.initAlbedo || globalParams.updateAlbedo || computeErrors ||
      (globalParams.useNormalizedCostFun && globalParams.computeWeightsSum) ||
      (useReflectance && globalParams.updateTilePhaseCoeffs)
      ){

    // Perform one of the several operations on the given tile
    
    if (isLastIter && globalParams.updateTilePhaseCoeffs){
      std::cout << "ERROR: Cannot update the tile phase coefficients in the last iteration." << std::endl;
      exit(1);
    }

    std::vector<ModelParams> overlapParamsArray(overlapIndicesArray.size());
    for (unsigned int j = 0; j < overlapIndicesArray.size(); j++){
      overlapParamsArray[j] = modelParamsArray[overlapIndicesArray[j]];
    }
    
    std::string albedoTileFile = currImgOrTile;
    std::string DEMTileFile    = meanDEMDir    + suffix_from_filename(albedoTileFile);
    std::string errorTileFile  = errorDir      + suffix_from_filename(albedoTileFile);
    std::string weightsSumFile = weightsSumDir + suffix_from_filename(albedoTileFile);
    std::string phaseTileFile  = phaseDir      + prefix_from_filename(suffix_from_filename(albedoTileFile)) + ".txt";
    phaseCoeffsData PCD;
    
    // If this is not the last albedo iteration, we must
    // init/update the exposure and albedo for all
    // tiles. Otherwise, do the update only if the tile overlaps
    // with the sim box, and wipe that tile altogether if it does
    // not.
    Vector4 tileCorners = getImageCorners(albedoTileFile);
    if ( isLastIter && !boxesOverlap(tileCorners, globalParams.simulationBox)){
      std::cout << "Skipping/removing tile: "
                << albedoTileFile << " as it does not overlap with the simulation box." << std::endl;
      try { boost::filesystem::remove(albedoTileFile); } catch(...){}
      return 0;
    }
    double costFunVal = actOnTile(isLastIter, computeErrors,
                                  DEMTileFile, albedoTileFile, errorTileFile, weightsSumFile,
                                  overlapParamsArray, globalParams, PCD);
    if (globalParams.updateTilePhaseCoeffs) PCD.writeToFile(phaseTileFile);
    if (!globalParams.initAlbedo && !globalParams.computeWeightsSum){
      std::string costFunFile = costFunDir
        + prefix_from_filename(suffix_from_filename(albedoTileFile)) + ".txt";
      AppendCostFunToFile(costFunVal, costFunFile);
    }
  }

  //re-estimate the height map  - shape from shading
  if (useReflectance && (updateHeight == 1)){

    std::vector<ModelParams> overlapParamsArray(overlapIndicesArray.size());
    for (unsigned int j = 0; j < overlapIndicesArray.size(); j++){
      overlapParamsArray[j] = modelParamsArray[overlapIndicesArray[j]];
    }
    
    std::string albedoTileFile = currImgOrTile;
    std::string DEMTileFile    = meanDEMDir + suffix_from_filename(albedoTileFile);
    std::string sfsTileFile    = sfsDir     + suffix_from_filename(albedoTileFile);
    UpdateHeightMapTiles(DEMTileFile, albedoTileFile, sfsTileFile,
                         overlapParamsArray, globalParams);
  }

  return 0;
}

