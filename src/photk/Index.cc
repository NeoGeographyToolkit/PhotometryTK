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
#include <fstream>
#include <sys/stat.h>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

//using namespace std;
//using namespace vw;
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions/gamma.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <photk/Reconstruct.h>
#include <photk/Index.h>

template <class ViewT>
std::vector<float> build_histogram(ImageViewBase<ViewT> const& view) {
        typedef typename ViewT::pixel_accessor pixel_accessor;
        std::vector<float> result(DYNAMIC_RANGE);

        int num_valid = 0;
        pixel_accessor row_acc = view.impl().origin();
        for( int32 row=0; row < view.impl().rows(); ++row ) {
                pixel_accessor col_acc = row_acc;
                for( int32 col=0; col < view.impl().cols(); ++col ) {
                        if ( is_valid(*col_acc) ) {
                                result[ (*col_acc)[0] ] += 1.0;
                                ++num_valid;
                        }
                        col_acc.next_col();
                }
                row_acc.next_row();
        }

        // Renormalize the histogram values
        for (size_t i=0; i< result.size(); ++i)
                result[i] /= num_valid;

        return result;
}

