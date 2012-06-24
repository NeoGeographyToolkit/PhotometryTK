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

bool compareDirsWithTol(string dir1, string dir2){

  // Given two output directories of running albedo, see if the albedo tiles
  // in those directories are similar enough
  
  int tol = 2; // image values go from 0 to 255, use 2 as the tolerance
  
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

    DiskImageView<PixelMask<PixelGray<uint8> > >  img1(file1);
    DiskImageView<PixelMask<PixelGray<uint8> > >  img2(file2);
    if (img1.rows() != img2.rows() || img1.cols() != img2.cols()){
      std::cout << "The following images have different sizes: " << file1  << ' ' << file2 << std::endl;
      return false;
    }

    for (int k = 0 ; k < img1.rows(); ++k) {
      for (int l = 0; l < img1.cols(); ++l) {
        if ( abs( (double)img1(l, k) - (double)img2(l, k) ) > tol ) return false;
      }
    }
    
  }

  return true;
}

class TestAlbedo : public ::testing::Test {
protected:
  TestAlbedo() {}
  virtual void SetUp() {}
};

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
    "photometry_settings_2.txt", "albedo_gold_2"
  };
  
  for (int s = 0 ; s < sizeof(files)/sizeof(string); s++){
    string cmd = "ln -s ../../../../src/photk/tests/" + files[s] + " . > /dev/null 2>&1";
    //std::cout << cmd << std::endl;
    system(cmd.c_str());
  }
  cmd = "ln -s ../../../../src/tools/reconstruct.sh . > /dev/null 2>&1";
  //std::cout << cmd << std::endl;
  int flag = system(cmd.c_str());
  
  // Run test 1
  std::cout << "\nRunning test 1" << std::endl;
  cmd="./reconstruct.sh photometry_settings_1.txt test_1 > output1.txt";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  flag = compareDirsWithTol("albedo_gold_1", "albedo_test_1");
  EXPECT_EQ(flag, 1);
  
  // Run test 2
  std::cout << "\nRunning test 2" << std::endl;
  cmd="./reconstruct.sh photometry_settings_2.txt test_2 > output2.txt";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  flag = compareDirsWithTol("albedo_gold_2", "albedo_test_2");
  EXPECT_EQ(flag, 1);
}

