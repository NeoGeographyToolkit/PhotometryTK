// //__BEGIN_LICENSE__
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


#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <photk/Photometry.h>
#include <test/Helpers.h>

using namespace vw;
using namespace photometry;
using namespace std;

template<class ImageT>
bool compareFilesWithTol(string file1, string file2, int tol){
  // Compare if two uint8 tif files have their pixels within tol of each other.
  std::cout << "Comparing files: " << file1 << ' ' << file2 << std::endl;
  DiskImageView<ImageT>  img1(file1);
  DiskImageView<ImageT>  img2(file2);
  if (img1.rows() != img2.rows() || img1.cols() != img2.cols()){
    std::cout << "The following images have different sizes: " << file1  << ' ' << file2 << std::endl;
    return false;
  }
  
  for (int k = 0 ; k < img1.rows(); ++k) {
    for (int l = 0; l < img1.cols(); ++l) {
      if ( abs( (double)img1(l, k) - (double)img2(l, k) ) > tol ) return false;
    }
  }
  return true;
}

bool compareDirsWithTol(string dir1, string dir2, int tol){

  // Given two output directories of running albedo, see if the albedo tiles
  // in those directories are similar enough
  
  vector<string> tifsInDir;
  listTifsInDir(dir1 + "/albedo", tifsInDir);
  for (int s = 0; s < tifsInDir.size(); s++){
    string file1 = tifsInDir[s];
    string file2 = dir2 + "/albedo" + suffix_from_filename(file1);
    
    std::cout << "Comparing: " << file1  << ' ' << file2 << std::endl;

    ifstream h1(file1.c_str());
    if (!h1) {
      std::cerr << "Missing file: " << file1 << std::endl;
      return false;
    }
    ifstream h2(file2.c_str());
    if (!h2) {
      std::cerr << "Missing file: " << file2 << std::endl;
      return false;
    }

    bool flag = compareFilesWithTol< PixelMask<PixelGray<uint8> > >(file1, file2, tol);
    if (!flag) return false;
  }

  return true;
}

class TestAlbedo : public ::testing::Test {
protected:
  TestAlbedo() {}
  virtual void SetUp() {}
};

// Function for expansion of a pre-processor macro into a string.
#define xstr(a) str(a)
#define str(a) #a

TEST_F( TestAlbedo, FullCycle ) {

  std::cout << std::endl;
  string cmd = "echo Running the tests in $(pwd)";
  system(cmd.c_str());
  
  // Make symbolic links in the current directory to the script and data necessary to run this test.
  string files[] = {
    "meta",
    "DIM_input_1280mpp_masked", "DIM_input_2560mpp",
    "DEM_tiles_sub64",
    "photometry_settings_1.txt", "albedo_gold_1",
    "photometry_settings_2.txt", "albedo_gold_2",
    "AS15-M-1134.lev1.cub", "apollo-DEM.tif",
    "demimg_iter1_int_gold.tif", "demrefl_iter1_int_gold.tif"
  };

  int tol = 2; // image values go from 0 to 255, use 2 as the tolerance
  
  std::string paths = xstr(PHOTK_SOURCE_DIR) + std::string("/build ") + xstr(VISIONWORKBENCH_ROOT);
  
  for (int s = 0 ; s < sizeof(files)/sizeof(string); s++){
    string cmd = std::string("rm -f ") + files[s] +
      "; ln -s ../../../../src/photk/tests/" + files[s] + " .";
    std::cout << cmd << std::endl;
    system(cmd.c_str());
  }
  cmd = "rm -f reconstruct.sh; ln -s ../../../../src/tools/reconstruct.sh .";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  
  // Run test 1
  std::cout << "\nRunning albedo test 1" << std::endl;
  system("echo Run directory is $(pwd)");
  cmd="./reconstruct.sh photometry_settings_1.txt test_1 " + paths + " > output1.txt 2>&1";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  int flag = compareDirsWithTol("albedo_gold_1", "albedo_test_1", tol);
  EXPECT_EQ(flag, 1);
  
  // Run test 2
  std::cout << "\nRunning albedo test 2" << std::endl;
  system("echo Run directory is $(pwd)");
  cmd="./reconstruct.sh photometry_settings_2.txt test_2 " + paths + " > output2.txt 2>&1";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  flag = compareDirsWithTol("albedo_gold_2", "albedo_test_2", tol);
  EXPECT_EQ(flag, 1);

#ifdef ENABLE_SFS
  
  std::cout << "\nRunning SfS test" << std::endl;
  system("echo Run directory is $(pwd)");

  if (!getenv("ISIS3DATA")){
    std::cerr << "Must set ISIS3DATA before running SFS" << std::endl;
    EXPECT_EQ(1, 0);
  }

  std::string cube = "AS15-M-1134.lev1.cub";

  cmd = std::string("export ISISROOT=") + xstr(ISIS_ROOT) + "; " + xstr(ISIS_ROOT)
    + "/bin/spiceinit from = " + cube;
  std::cout << cmd << std::endl;
  system(cmd.c_str());

  cmd = std::string("export ISISROOT=") + xstr(ISIS_ROOT) + "; "
    + xstr(PHOTK_SOURCE_DIR) + std::string("/build/src/tools/sfs ") + cube
    + " apollo-DEM.tif meta 1";
  std::cout << cmd << std::endl;
  system(cmd.c_str());

  flag = compareFilesWithTol<uint8>("demimg_iter1_int.tif",  "demimg_iter1_int_gold.tif", tol);
  EXPECT_EQ(flag, 1);
  
  flag = compareFilesWithTol<uint8>("demrefl_iter1_int.tif", "demrefl_iter1_int_gold.tif", tol);
  EXPECT_EQ(flag, 1);
  
#endif
  
}

