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

// To do:
// Copy the DIM data to lou
// Simplify the script.
// Use normalized weights in shape-from-shading as well.
// More work in shape from shading: Strip padding at the last iteration.
// Must regenerate the index files all the time, as they are too fragile
// The isis adjust file may not exist.
// Copy the images from supercomp.
// Rename dem_out.tif to dem_mean.tif, and modelParams.outputFile to modelParams.albedoFile,
// inputFile to drgFile.
// Merge the imageRecord and modelParams classes
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
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

  time_t Start_t = time(NULL);
 
  std::vector<std::string> inputDRGFiles;
  std::string resDir, settingsFile, imagesList, albedoTilesList;
  std::string currImgOrTile;
  
  po::options_description general_options("Options");
  general_options.add_options()
    ("settings-file,s", po::value<std::string>(&settingsFile), "Settings file")
    ("results-directory,r", po::value<std::string>(&resDir),   "Results directory")
    ("images-list,f", po::value<std::string>(&imagesList),     "The list of images")
    ("tiles-list,t", po::value<std::string>(&albedoTilesList), "The list of albedo tiles")
    ("image-file,i", po::value<std::string>(&currImgOrTile),   "Current image or tile")
    ("initial-setup",            "Initial setup")
    ("save-weights",             "Save the weights")
    ("compute-weights-sum",      "Compute the sum of weights at each pixel")
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

  // Not using reflectance is the same as creating a mosaic instead of albedo
  bool useReflectance = (globalParams.reflectanceType != NO_REFL); 

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
  if ( useReflectance &&
       !fs::exists(meanDEMDir   ) ) fs::create_directories(meanDEMDir   );
  if ( !fs::exists(costFunDir   ) ) fs::create_directories(costFunDir   );
  if ( computeErrors &&
       !fs::exists(errorDir     ) ) fs::create_directories(errorDir     );
  if ( globalParams.useNormalizedCostFun &&
       !fs::exists(weightsSumDir) ) fs::create_directories(weightsSumDir);
  if ( updateHeight &&
       !fs::exists(sfsDir ) )       fs::create_directories(sfsDir       );
  if ( !fs::exists(phaseDir     ) ) fs::create_directories(phaseDir     );

  // This will overwrite the values of globalParams.phaseCoeffC1, globalParams.phaseCoeffC2
  // if that information is available on disk.
  ReadPhaseCoeffsFromFile(phaseDir, globalParams);
  if (globalParams.initialSetup == 1) AppendPhaseCoeffsToFile(phaseDir, globalParams);

  // The names of the files listing all DRGs and DEMs and the coordinates
  // of their corners. 
  std::string allDRGIndex  = globalParams.drgDir + "/index.txt";
  std::string allDEMIndex  = globalParams.demDir + "/index.txt";
  std::string DEMTilesList = resDir              + "/DEMTilesList.txt";

  // A tile which will hold the georef info for all images and tiles
  std::string sampleTileFile = resDir + "/sampleTile.tif";
  
  if (globalParams.initialSetup == 1) {

    // This block of code creates the list of images (imagesList). As such, it must be above
    // any code which reads the list of images.
    // Create the imagesList used in subsequent iterations.
    // Create the list of all DEM if not there yet.
    listDRGinBoxAndAllDEM(useReflectance,
                          allDRGIndex, allDEMIndex,
                          globalParams.simulationBox, globalParams.drgDir, globalParams.demDir, imagesList
                          );

    vw_out( VerboseDebugMessage, "photometry" ) << "Initializing the albedo tiles ... ";
    std::vector<ImageRecord> drgRecords;
    if (!readImagesFile(drgRecords, imagesList)) exit(1);
    std::string imageFile = currImgOrTile; // an image whose georef we will use
    listAlbedoTilesOverlappingWithDRG(globalParams.tileSize, globalParams.pixelPadding,
                                      imageFile, globalParams.simulationBox,
                                      drgRecords,
                                      DEMTilesList,    meanDEMDir,
                                      albedoTilesList, albedoDir,
                                      sampleTileFile
                                      );
  }

  std::vector<ImageRecord> drgRecords;
  if (!readImagesFile(drgRecords, imagesList)) exit(1);

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
  for (unsigned int j = 0; j < drgRecords.size(); ++j) {
      
    DRGFiles.push_back(drgRecords[j].path);
    modelParamsArray.push_back(ModelParams());
      
    std::string temp = suffix_from_filename(DRGFiles[j]);
    modelParamsArray[j].exposureTime     = 1.0;
    modelParamsArray[j].hCenterLineDEM   = NULL;
    modelParamsArray[j].hMaxDistArrayDEM = NULL;
    modelParamsArray[j].inputFilename    = DRGFiles[j];//these filenames have full path
      
    std::string prefix = getFirstElevenCharsFromFileName(DRGFiles[j]);
      
    if (useReflectance){
      if ( sunPositions.find(prefix) == sunPositions.end()){
        std::cerr << "Could not find the sun position for the DRG file: " << DRGFiles[j] << std::endl;
        exit(1);
      }
      modelParamsArray[j].sunPosition = 1000*sunPositions[prefix];
        
      if (spacecraftPositions.find(prefix) == spacecraftPositions.end()){
        std::cerr << "Could not find the spacecraft position for the DRG file: " << DRGFiles[j] << std::endl;
        exit(1);
      }
      // Go from kilometers to meters
      modelParamsArray[j].spacecraftPosition = 1000*spacecraftPositions[prefix];
    }

     modelParamsArray[j].infoFilename        = resDir + "/info/"       + prefix_less3_from_filename(temp) + "info.txt";
     modelParamsArray[j].reliefFilename      = resDir + "/reflectance" + prefix_from_filename(temp)       + "_reflectance.tif";
     modelParamsArray[j].errorFilename       = resDir + "/error"       + prefix_from_filename(temp)       + "_err.tif";
     modelParamsArray[j].outputFilename      = resDir + "/albedo"      + prefix_from_filename(temp)       + "_albedo.tif";
     modelParamsArray[j].sfsDEMFilename      = resDir + "/DEM_sfs"     + prefix_less3_from_filename(temp) + "DEM_sfs.tif";
     modelParamsArray[j].errorHeightFilename = resDir + "/error"       + prefix_from_filename(temp)       + "_height_err.tif";
     modelParamsArray[j].weightFilename      = resDir + "/weight"      + prefix_from_filename(temp)       + "_weight.txt";
     modelParamsArray[j].exposureFilename    = resDir + "/exposure"    + prefix_from_filename(temp)       + "_exposure.txt";
      
    const ImageRecord& rec = drgRecords[j];
    modelParamsArray[j].corners = Vector4(rec.west, rec.east, rec.south, rec.north);
      
    modelParamsArray[j].hCenterLineDEM   = NULL;
    modelParamsArray[j].hMaxDistArrayDEM = NULL;
    modelParamsArray[j].vCenterLineDEM   = NULL;
    modelParamsArray[j].vMaxDistArrayDEM = NULL;
      
    vw_out( VerboseDebugMessage, "photometry" ) << modelParamsArray[j] << "\n";
  }
  
  if (globalParams.initialSetup == 1){
    // We performed all tasks, including validation, if we are in the initial setup.
    return 0;
  }
    
  vw_out() << "Number of Files = " << DRGFiles.size() << "\n";

  // Get the corners for currImgOrTile from the list
  // of images or tiles
  std::vector<int> inputIndices = GetInputIndices(currImgOrTile, drgRecords);
  ImageRecord currImgOrTileCorners;
  if (inputIndices.size() != 0){
    // currImgOrTile is an image
    currImgOrTileCorners = drgRecords[inputIndices[0]];
  }else{
    // currImgOrTile is a tile
    std::vector<ImageRecord> albedoTiles;
    if (!readImagesFile(albedoTiles, albedoTilesList)) exit(1);
    std::vector<int> inputIndicesLocal = GetInputIndices(currImgOrTile, albedoTiles);
    if (inputIndicesLocal.size() >= 1){
      currImgOrTileCorners = albedoTiles[inputIndicesLocal[0]];
    }else if(drgRecords.size() >= 1){
      // In this case we don't really use currImgOrTileCorners but do this for  consistency
      currImgOrTileCorners = drgRecords[0];
    }else{
      cerr << "ERROR: No images." << endl;
      exit(1);
    }
  }

  // The images overlapping with currImgOrTile
  std::vector<int>  overlapIndicesArray = makeOverlapList(drgRecords, currImgOrTileCorners);
  printOverlapList(overlapIndicesArray);
  
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
    std::string currTile        = currImgOrTile;
    ImageRecord currTileCorners = currImgOrTileCorners;
    std::string DEMTileFile     = meanDEMDir + suffix_from_filename(currTile);
    std::vector<ImageRecord> DEMImages;
    if (!readImagesFile(DEMImages, allDEMIndex)) exit(1);
    std::vector<int> overlap = makeOverlapList(DEMImages, currTileCorners);
    InitMeanDEMTile(sampleTileFile, currTileCorners, DEMTileFile, 
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
    ImageRecord currImgCorners = currImgOrTileCorners;
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
    overlap = makeOverlapList(DEMTiles, currImgCorners);
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
      A1_num += PCD.phaseCoeffC1_num; A1_den += PCD.phaseCoeffC1_den;
      A2_num += PCD.phaseCoeffC2_num; A2_den += PCD.phaseCoeffC2_den;
    }
    if (A1_den != 0.0) globalParams.phaseCoeffC1 += A1_num/A1_den;
    if (A2_den != 0.0) globalParams.phaseCoeffC2 += A2_num/A2_den;
    AppendPhaseCoeffsToFile(phaseDir, globalParams);
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
    
    ImageRecord albedoTileCorners = currImgOrTileCorners;
    std::string albedoTileFile    = currImgOrTile;
    std::string DEMTileFile       = meanDEMDir    + suffix_from_filename(albedoTileFile);
    std::string errorTileFile     = errorDir      + suffix_from_filename(albedoTileFile);
    std::string weightsSumFile    = weightsSumDir + suffix_from_filename(albedoTileFile);
    std::string phaseTileFile     = phaseDir      + prefix_from_filename(suffix_from_filename(albedoTileFile)) + ".txt";
    phaseCoeffsData PCD;
    
    // If this is not the last albedo iteration, we must
    // init/update the exposure and albedo for all
    // tiles. Otherwise, do the update only if the tile overlaps
    // with the sim box, and wipe that tile altogether if it does
    // not.
    Vector4 cornersVec = Vector4(albedoTileCorners.west, albedoTileCorners.east,
                                 albedoTileCorners.south, albedoTileCorners.north);
    if ( isLastIter && !boxesOverlap(cornersVec, globalParams.simulationBox)){
      std::cout << "Removing tile: "
                << albedoTileFile << " as it does not overlap with the simulation box." << std::endl;
      try { boost::filesystem::remove(albedoTileFile); } catch(...){}
      return 0;
    }
    double costFunVal = actOnTile(isLastIter, computeErrors,
                                  sampleTileFile, albedoTileCorners,
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

  time_t End_t = time(NULL);
  double time_task = difftime(End_t, Start_t);
  ostringstream cmd;
  cmd << "echo \"\n\n Job " << currImgOrTile << " at $(date) on machine $(uname -n) took " << setw(3) << time_task << " seconds.\n\n\"" << endl;
  //std::cout << cmd.str() << std::endl;
  system(cmd.str().c_str());
  
  return 0;
}

