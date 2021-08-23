/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFixedPointVolumeRayCastCompositeShadeHelper.cxx
  Language:  C++

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkFixedPointVolumeRayCastCompositeShadeHelper.h"

#include "vtkCommand.h"
#include "vtkDataArray.h"
#include "vtkFixedPointRayCastImage.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkRectilinearGrid.h"
#include "vtkRenderWindow.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include <typeinfo>
#include <bitset>

#include <fstream>
#include <chrono>
// #include "../../../mfa/include/mfa/mfa.hpp"
// #include "../../../mfa/include/diy/include/diy/master.hpp"
// #include "../../../mfa/include/diy/include/diy/io/block.hpp"

// #include    "/Users/jianxinsun/research/intern/VTK/temp/opts.h"
// #include    "/Users/jianxinsun/research/intern/VTK/temp/block.hpp"

// #include <mfa/mfa.hpp>
#include "/Users/jianxinsun/research/intern/mfa/include/mfa/mfa.hpp"
// #include <diy/master.hpp>
#include "/Users/jianxinsun/research/intern/mfa/include/diy/include/diy/master.hpp"

#include "/Users/jianxinsun/research/intern/VTK/external/opts.h"
#include "/Users/jianxinsun/research/intern/VTK/external/block.hpp"

#include <cmath>

// sinc 20-11
// #define LENGTH 19.0

// sinc 200-20
#define LENGTH 199.0

#define MIN -2.0
#define MAX 9.0

vtkStandardNewMacro(vtkFixedPointVolumeRayCastCompositeShadeHelper);

// Construct a new vtkFixedPointVolumeRayCastCompositeShadeHelper with default values
vtkFixedPointVolumeRayCastCompositeShadeHelper::vtkFixedPointVolumeRayCastCompositeShadeHelper() =
  default;

// Destruct a vtkFixedPointVolumeRayCastCompositeShadeHelper - clean up any memory used
vtkFixedPointVolumeRayCastCompositeShadeHelper::~vtkFixedPointVolumeRayCastCompositeShadeHelper() =
  default;

// This method is used when the interpolation type is nearest neighbor and
// the data has one component and scale == 1.0 and shift == 0.0. In the inner
// loop we get the data value as an unsigned short, and use this index to
// lookup a color and opacity for this sample. We then composite this into
// the color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleShadeNN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val = static_cast<unsigned short>(((*dptr)));
    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp);
    if (tmp[3])
    {
      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading(diffuseShadingTable[0], specularShadingTable[0], normal, tmp);
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has one component. In the inner loop we get the data value as
// an unsigned short using the scale/shift, and use this index to lookup
// a color and opacity for this sample. We then composite this into the
// color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageOneNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleShadeNN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val = static_cast<unsigned short>(((*dptr) + shift[0]) * scale[0]);
    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp);
    if (tmp[3])
    {
      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading(diffuseShadingTable[0], specularShadingTable[0], normal, tmp);
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has two components which are not considered independent. In the
// inner loop we compute the two unsigned short index values from the data
// values (using the scale/shift). We use the first index to lookup a color,
// and we use the second index to look up the opacity. We then composite
// the color into the color computed so far along this ray, and check to
// see if we can terminate here (if the opacity accumulated exceed some
// threshold). Finally we move to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageTwoDependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleShadeNN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val[2];
    val[0] = static_cast<unsigned short>(((*(dptr)) + shift[0]) * scale[0]);
    val[1] = static_cast<unsigned short>(((*(dptr + 1)) + shift[1]) * scale[1]);

    tmp[3] = scalarOpacityTable[0][val[1]];
    if (tmp[3])
    {
      tmp[0] = static_cast<unsigned short>(
        (colorTable[0][3 * val[0]] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
      tmp[1] = static_cast<unsigned short>(
        (colorTable[0][3 * val[0] + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
      tmp[2] = static_cast<unsigned short>(
        (colorTable[0][3 * val[0] + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading(diffuseShadingTable[0], specularShadingTable[0], normal, tmp);
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has four components which are not considered independent . This
// means that the first three components directly represent color, and this
// data must be of unsigned char type. In the inner loop we directly access
// the four data values (no scale/shift is needed). The first three are the
// color of this sample and the fourth is used to look up an opacity in the
// scalar opacity transfer function. We then composite this color into the
// color we have accumulated so far along the ray, and check if we can
// terminate here (if our accumulated opacity has exceed some threshold).
// Finally we move onto the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageFourDependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleShadeNN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val[4];
    val[0] = *(dptr);
    val[1] = *(dptr + 1);
    val[2] = *(dptr + 2);
    val[3] = static_cast<unsigned short>(((*(dptr + 3)) + shift[3]) * scale[3]);

    tmp[3] = scalarOpacityTable[0][val[3]];
    if (tmp[3])
    {
      tmp[0] = (val[0] * tmp[3] + 0x7f) >> (8);
      tmp[1] = (val[1] * tmp[3] + 0x7f) >> (8);
      tmp[2] = (val[2] * tmp[3] + 0x7f) >> (8);

      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading(diffuseShadingTable[0], specularShadingTable[0], normal, tmp);
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has more than one component and the components are considered to
// be independent. In the inner loop we access each component value, using
// the scale/shift to turn the data value into an unsigned short index. We
// then lookup the color/opacity for each component and combine them according
// to the weighting value for each component. We composite this resulting
// color into the color already accumulated for this ray, and we check
// whether we can terminate here (if the accumulated opacity exceeds some
// threshold). Finally we increment to the next sample on the ray.
//
// TODO: short circuit calculations when opacity is 0
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageIndependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeMultiNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleShadeNN();
    }

    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short normal[4];
    for (c = 0; c < components; c++)
    {
      val[c] = static_cast<unsigned short>(((*(dptr + c)) + shift[c]) * scale[c]);
      normal[c] = *(dirPtr + c);
    }

    VTKKWRCHelper_LookupAndCombineIndependentColorsShadeUS(colorTable, scalarOpacityTable,
      diffuseShadingTable, specularShadingTable, val, normal, weights, components, tmp);

    if (tmp[3])
    {
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

int count = 0;
int count1 = 0;
int minn = 255;
int maxx = 0;
long int time0 = 0;
long int time1 = 0;
long int time2 = 0;
long int time3 = 0;
long int time4 = 0;
long int time5 = 0;
long int time6 = 0;
long int time7 = 0;

long int mfa_time0 = 0;
long int mfa_time1 = 0;
long int mfa_time2 = 0;
long int mfa_time3 = 0;
long int mfa_time4 = 0;
long int mfa_time5 = 0;
long int mfa_time6 = 0;
long int mfa_time7 = 0;

long int mfa_counter0 = 0;
long int mfa_counter1 = 0;
long int mfa_counter2 = 0;
long int mfa_counter3 = 0;
long int mfa_counter4 = 0;
long int mfa_counter5 = 0;
long int mfa_counter6 = 0;
long int mfa_counter7 = 0;

long int test_count = 0;

// This method is used when the interpolation type is linear and the data
// has one component and scale = 1.0 and shift = 0.0. In the inner loop we
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (the range must be right and we don't need the
// scale/shift). We compute our weights within the cell according to our
// fractional position within the cell, apply trilinear interpolation to
// compute the index, and use this index to lookup a color and opacity for
// this sample. We then composite this into the color computed so far along
// the ray, and check if we can terminate at this point (if the accumulated
// opacity is higher than some threshold). Finally we move on to the next
// sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol) // handling batch of image pixels assigned to this thread
{
 
  /* Retrieve MFA reference pointer */ 
  Block<real_t>* b = (Block<real_t>*)(mapper->GetMfaBlock());
  mfa::MFA_Data<real_t>& mfa_data = *(b->vars[0].mfa_data);
  TensorProduct<real_t>&  t = mfa_data.tmesh.tensor_prods[0];
  mfa::Decoder<real_t> decoder(mfa_data, 0);
  mfa::FastDecodeInfo<real_t> di(decoder);
  /* Retrieve dataset size and update scalar for POS */
  int size =  mapper->GetMfaSize() - 1;
  std::string dataset =  mapper->GetDataset();
  double recip = 1 / ((size*32767.0) + 0.5);

  /* Time-stamp initialization */
  std::chrono::time_point<std::chrono::system_clock> before;
  std::chrono::time_point<std::chrono::system_clock> after;
  std::chrono::time_point<std::chrono::system_clock> mfa_before;
  std::chrono::time_point<std::chrono::system_clock> mfa_after;
  if (threadID == 0) {
    time0 = 0;
    mfa_time0 = 0;
    mfa_counter0 = 0;
  }
  if (threadID == 1) {
    time1 = 0;
    mfa_time1 = 0;
    mfa_counter1 = 0;
  }
  if (threadID == 2) {
    time2 = 0;
    mfa_time2 = 0;
    mfa_counter2 = 0;
  }
  if (threadID == 3) {
    time3 = 0;
    mfa_time3 = 0;
    mfa_counter3 = 0;
  }
  if (threadID == 4) {
    time4 = 0;
    mfa_time4 = 0;
    mfa_counter4 = 0;
  }
  if (threadID == 5) {
    time5 = 0;
    mfa_time5 = 0;
    mfa_counter5 = 0;
  }
  if (threadID == 6) {
    time6 = 0;
    mfa_time6 = 0;
    mfa_counter6 = 0;
  }
  if (threadID == 7) {
    time7 = 0;
    mfa_time7 = 0;
    mfa_counter7 = 0;
  }



  VectorX<double> param(3); 
  VectorX<double> out_pt(1);
  VectorX<double> out_pt4(4);


  // cerr << "AAA" << threadID << endl;
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  
  int skip_factor = 1; 
  for (k = 0; k < numSteps; k++)
  {
    if (k%skip_factor != 0) {
        continue;
    }

    if (k)
    {
      mapper->FixedPointIncrement(pos, dir); // increment a dir on pos
      if (threadID == 0) {
      }
    }
    // VTKKWRCHelper_SpaceLeapCheck(); // contine to next iteration if check fails
    VTKKWRCHelper_CroppingCheckTrilin(pos); // contine if check cropping

    mapper->ShiftVectorDown(pos, spos); // only keep the first 15 bits of the 32bit int
        if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
        {
          oldSPos[0] = spos[0];
          oldSPos[1] = spos[1];
          oldSPos[2] = spos[2];

          dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
          VTKKWRCHelper_GetCellScalarValuesSimple(dptr); // update 8 corners of the cell
          dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
          dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
          needToSampleDirection = 1;
        }

    if (!mapper->GetUseMfa()) {
        before = std::chrono::system_clock::now();
        VTKKWRCHelper_ComputeWeights(pos);
        VTKKWRCHelper_InterpolateScalar(val); // val is unsigned short, update val, val is index

        /* benchmark time used */
        after = std::chrono::system_clock::now();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(after - before);
        if (threadID == 0) {
            time0 += nanoseconds.count();
        }
        if (threadID == 1) {
            time1 += nanoseconds.count();
        }
        if (threadID == 2) {
            time2 += nanoseconds.count();
        }
        if (threadID == 3) {
            time3 += nanoseconds.count();
        }
        if (threadID == 4) {
            time4 += nanoseconds.count();
        }
        if (threadID == 5) {
            time5 += nanoseconds.count();
        }
        if (threadID == 6) {
            time6 += nanoseconds.count();
        }
        if (threadID == 7) {
            time7 += nanoseconds.count();
        }
    } else {
        mfa_before = std::chrono::system_clock::now();
        VTKKWRCHelper_ComputeWeights(pos);
        param(0) = pos[0] * recip;
        param(1) = pos[1] * recip;
        param(2) = pos[2] * recip;
        // param(0) = 0.5;
        // param(1) = 0.5;
        // param(2) = 0.5;
        /* Do MFA decoding */
        if (threadID == 0) {
            // cerr << "before decoder" << endl;
        }
        decoder.FastVolPt(param, out_pt, di, t);
        // b->my_decode_point(param, out_pt4);
        if (threadID == 0) {
            // cerr << "after decoder" << endl;
        }
        /* Normalize value */
        double v = out_pt[0];
        // double v = out_pt[3];
        double ratio;
        double bounds_lower = b->bounds_mins[3];
        // double bounds_lower = -14.6001;
        double bounds_upper = b->bounds_maxs[3];
        // double bounds_upper = 130.085;
        // /*
        if (threadID == 0) {
            // cerr << "val" << v << endl;
            // cerr << "min: " << bounds_lower << endl;
            // cerr << "max: " << bounds_upper << endl;
            if (test_count == 0) {
                cerr << param << endl;
                cerr << "val: " << v << endl;
            }
            test_count++;
        }
        // */
        ratio = (v - bounds_lower)/(bounds_upper - bounds_lower);
        /*
        if (dataset == "sinc") {
            ratio = (v + 2.0)/11.0625; // for sinc dataset
            // ratio = (v + 2.16934)/(9.98008 + 2.16934); // for sinc 200-20 dataset
            // ratio = (v + 1.85646)/(8.00917 + 1.85646); // for sinc 20-11 dataset
        }
        if (dataset == "nek5000") {
            ratio = v/135.111; // for nek5000 dataset
        }
        */
        if (ratio > 1) {
            ratio = 1;
        }
        if (ratio < 0) {
            ratio = 0;
        }
        unsigned short new_val = (unsigned short)(255*ratio);
        val = new_val;

        /* benchmark time used */
        mfa_after = std::chrono::system_clock::now();
        auto mfa_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(mfa_after - mfa_before);
        if (threadID == 0) {
            mfa_time0 += mfa_nanoseconds.count();
            mfa_counter0++;
        }
        if (threadID == 1) {
            mfa_time1 += mfa_nanoseconds.count();
            mfa_counter1++;
        }
        if (threadID == 2) {
            mfa_time2 += mfa_nanoseconds.count();
            mfa_counter2++;
        }
        if (threadID == 3) {
            mfa_time3 += mfa_nanoseconds.count();
            mfa_counter3++;
        }
        if (threadID == 4) {
            mfa_time4 += mfa_nanoseconds.count();
            mfa_counter4++;
        }
        if (threadID == 5) {
            mfa_time5 += mfa_nanoseconds.count();
            mfa_counter5++;
        }
        if (threadID == 6) {
            mfa_time6 += mfa_nanoseconds.count();
            mfa_counter6++;
        }
        if (threadID == 7) {
            mfa_time7 += mfa_nanoseconds.count();
            mfa_counter7++;
        }
    }
    /* Do Compositing */
    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp); // update tmp, tmp[0,1,2] unsigned short
    if (needToSampleDirection) // true
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }
    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp); // update tmp again
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity); // break for loop if terminates
    unsigned short* imagePtr;                                                                        \
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
 
  /* Save Time-stamp */ 
  if (threadID == 0) {
      std::ofstream filetime0;
      filetime0.open("time0.txt");
      filetime0 << time0;
      filetime0.close();

      std::ofstream mfa_filetime0;
      mfa_filetime0.open("mfa_time0.txt");
      mfa_filetime0 << mfa_time0;
      mfa_filetime0.close();

      std::ofstream mfa_filecount0;
      mfa_filecount0.open("mfa_count0.txt");
      mfa_filecount0 << mfa_counter0;
      mfa_filecount0.close();
  }
  if (threadID == 1) {
      std::ofstream filetime1;
      filetime1.open("time1.txt");
      filetime1 << time1;
      filetime1.close();

      std::ofstream mfa_filetime1;
      mfa_filetime1.open("mfa_time1.txt");
      mfa_filetime1 << mfa_time1;
      mfa_filetime1.close();

      std::ofstream mfa_filecount1;
      mfa_filecount1.open("mfa_count1.txt");
      mfa_filecount1 << mfa_counter1;
      mfa_filecount1.close();
  }
  if (threadID == 2) {
      std::ofstream filetime2;
      filetime2.open("time2.txt");
      filetime2 << time2;
      filetime2.close();

      std::ofstream mfa_filetime2;
      mfa_filetime2.open("mfa_time2.txt");
      mfa_filetime2 << mfa_time2;
      mfa_filetime2.close();

      std::ofstream mfa_filecount2;
      mfa_filecount2.open("mfa_count2.txt");
      mfa_filecount2 << mfa_counter2;
      mfa_filecount2.close();
  }
  if (threadID == 3) {
      std::ofstream filetime3;
      filetime3.open("time3.txt");
      filetime3 << time3;
      filetime3.close();
      // cerr << "time3: " << time3 << endl;

      std::ofstream mfa_filetime3;
      mfa_filetime3.open("mfa_time3.txt");
      mfa_filetime3 << mfa_time3;
      mfa_filetime3.close();
      // cerr << "mfa_time3: " << mfa_time3 << endl;

      std::ofstream mfa_filecount3;
      mfa_filecount3.open("mfa_count3.txt");
      mfa_filecount3 << mfa_counter3;
      mfa_filecount3.close();
      // cerr << "mfa_count3: " << mfa_counter3 << endl;
  }
  if (threadID == 4) {
      std::ofstream filetime4;
      filetime4.open("time4.txt");
      filetime4 << time4;
      filetime4.close();

      std::ofstream mfa_filetime4;
      mfa_filetime4.open("mfa_time4.txt");
      mfa_filetime4 << mfa_time4;
      mfa_filetime4.close();

      std::ofstream mfa_filecount4;
      mfa_filecount4.open("mfa_count4.txt");
      mfa_filecount4 << mfa_counter4;
      mfa_filecount4.close();
  }
  if (threadID == 5) {
      std::ofstream filetime5;
      filetime5.open("time5.txt");
      filetime5 << time5;
      filetime5.close();

      std::ofstream mfa_filetime5;
      mfa_filetime5.open("mfa_time5.txt");
      mfa_filetime5 << mfa_time5;
      mfa_filetime5.close();

      std::ofstream mfa_filecount5;
      mfa_filecount5.open("mfa_count5.txt");
      mfa_filecount5 << mfa_counter5;
      mfa_filecount5.close();
  }
  if (threadID == 6) {
      std::ofstream filetime6;
      filetime6.open("time6.txt");
      filetime6 << time6;
      filetime6.close();

      std::ofstream mfa_filetime6;
      mfa_filetime6.open("mfa_time6.txt");
      mfa_filetime6 << mfa_time6;
      mfa_filetime6.close();

      std::ofstream mfa_filecount6;
      mfa_filecount6.open("mfa_count6.txt");
      mfa_filecount6 << mfa_counter6;
      mfa_filecount6.close();
  }
  if (threadID == 7) {
      std::ofstream filetime7;
      filetime7.open("time7.txt");
      filetime7 << time7;
      filetime7.close();

      std::ofstream mfa_filetime7;
      mfa_filetime7.open("mfa_time7.txt");
      mfa_filetime7 << mfa_time7;
      mfa_filetime7.close();

      std::ofstream mfa_filecount7;
      mfa_filecount7.open("mfa_count7.txt");
      mfa_filecount7 << mfa_counter7;
      mfa_filecount7.close();
  }
}

// This method is used when the interpolation type is linear and the data
// has one component and scale != 1.0 or shift != 0.0. In the inner loop we
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (we use the scale/shift to ensure the correct range).
// We compute our weights within the cell according to our fractional position
// within the cell, apply trilinear interpolation to compute the index, and use
// this index to lookup a color and opacity for this sample. We then composite
// this into the color computed so far along the ray, and check if we can
// terminate at this point (if the accumulated opacity is higher than some
// threshold). Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageOneTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValues(dptr, scale[0], shift[0]);
      dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
      dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
      needToSampleDirection = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);

    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp);
    if (needToSampleDirection)
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }
    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp);
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// two components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for both components as an unsigned shorts (we use the
// scale/shift to ensure the correct range). We compute our weights within
// the cell according to our fractional position within the cell, and apply
// trilinear interpolation to compute the two index value. We use the first
// index to lookup a color and the second to look up an opacity for this sample.
// We then composite this into the color computed so far along the ray, and
// check if we can terminate at this point (if the accumulated opacity is
// higher than some threshold). Finally we move on to the next sample along
// the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageTwoDependentTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 0, scale[0], shift[0]);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 1, scale[1], shift[1]);

      dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
      dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
      needToSampleDirection = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, 2);

    tmp[3] = scalarOpacityTable[0][val[1]];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleDirection)
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }

    tmp[0] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0]] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp);
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// four components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we don't have to
// use the scale/shift because only unsigned char data is supported for four
// component data when the components are not independent). We compute our
// weights within the cell according to our fractional position within the cell,
// and apply trilinear interpolation to compute a value for each component. We
// use the first three directly as the color of the sample, and the fourth is
// used to look up an opacity for this sample. We then composite this into the
// color computed so far along the ray, and check if we can terminate at this
// point (if the accumulated opacity is higher than some threshold). Finally we
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageFourDependentTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 0);

      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 1);

      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 2);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);

      dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
      dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
      needToSampleDirection = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, 4);

    tmp[3] = scalarOpacityTable[0][val[3]];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleDirection)
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }

    tmp[0] = (val[0] * tmp[3] + 0x7f) >> 8;
    tmp[1] = (val[1] * tmp[3] + 0x7f) >> 8;
    tmp[2] = (val[2] * tmp[3] + 0x7f) >> 8;

    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp);
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// more than one component and the components are considered independent. In
// the inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we have to use the
// scale/shift to ensure that we obtained unsigned short indices) We compute our
// weights within the cell according to our fractional position within the cell,
// and apply trilinear interpolation to compute a value for each component. We
// look up a color/opacity for each component and blend them according to the
// component weights. We then composite this resulting color into the
// color computed so far along the ray, and check if we can terminate at this
// point (if the accumulated opacity is higher than some threshold). Finally we
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageIndependentTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeMultiShadeTrilin();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 0, scale[0], shift[0]);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 1, scale[1], shift[1]);

      if (components > 2)
      {
        dptr++;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, 2, scale[2], shift[2]);
        if (components > 3)
        {
          dptr++;
          VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);
        }
      }

      dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
      dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
      VTKKWRCHelper_GetCellComponentDirectionValues(dirPtrABCD, dirPtrEFGH, 0);

      dirPtrABCD++;
      dirPtrEFGH++;
      VTKKWRCHelper_GetCellComponentDirectionValues(dirPtrABCD, dirPtrEFGH, 1);

      if (components > 2)
      {
        dirPtrABCD++;
        dirPtrEFGH++;
        VTKKWRCHelper_GetCellComponentDirectionValues(dirPtrABCD, dirPtrEFGH, 2);
        if (components > 3)
        {
          dirPtrABCD++;
          dirPtrEFGH++;
          VTKKWRCHelper_GetCellComponentDirectionValues(dirPtrABCD, dirPtrEFGH, 3);
        }
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);

    VTKKWRCHelper_LookupAndCombineIndependentColorsInterpolateShadeUS(colorTable,
      scalarOpacityTable, diffuseShadingTable, specularShadingTable, val, weights, components, tmp);

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

void vtkFixedPointVolumeRayCastCompositeShadeHelper::GenerateImage(
  int threadID, int threadCount, vtkVolume* vol, vtkFixedPointVolumeRayCastMapper* mapper)
{
  void* data = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType(); // VTK_UNSIGNED_CHAR   3

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol)) // false
  {
    // cerr << "use Nearest Neighbor Interpolation" << endl;
    // One component data
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      else
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
    }
    // More that one independent components
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageIndependentNN(
          static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent (color) components
    else
    {
      // Two components - the first specifies color (through a lookup table)
      // and the second specified opacity (through a lookup table)
      if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 2)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageTwoDependentNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup
      // table)
      else
      {
        if (scalarType == VTK_UNSIGNED_CHAR)
        {
          vtkFixedPointCompositeShadeHelperGenerateImageFourDependentNN(
            static_cast<unsigned char*>(data), threadID, threadCount, mapper, vol);
        }
        else
        {
          vtkErrorMacro("Four component dependent data must be unsigned char");
        }
      }
    }
  }
  // Trilinear Interpolation
  else
  {
    // cerr << "use Trilinear Interpolation" << endl;
    // cerr << "mapper->GetCurrentScalars()->GetNumberOfComponents() == 1" << mapper->GetCurrentScalars()->GetNumberOfComponents() << endl;
    // One component
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1) // true
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      // if (threadID == 0) {
      //   cerr << "mapper->GetTableScale()[0]" << mapper->GetTableScale()[0] << endl; 
      //   cerr << "mapper->GetTableShift()[0]" << mapper->GetTableShift()[0] << endl;
      // }
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0) // true
      {
        // cerr << "in if" << endl;
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
      {
        // cerr << "in else" << endl;
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
    }
    // Independent components (more than one)
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageIndependentTrilin(
          static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent components
    else
    {
      // Two components - the first specifies color (through a lookup table)
      // and the second specified opacity (through a lookup table)
      if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 2)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageTwoDependentTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup
      // table)
      else
      {
        if (scalarType == VTK_UNSIGNED_CHAR)
        {
          vtkFixedPointCompositeShadeHelperGenerateImageFourDependentTrilin(
            static_cast<unsigned char*>(data), threadID, threadCount, mapper, vol);
        }
        else
        {
          vtkErrorMacro("Four component dependent data must be unsigned char");
        }
      }
    }
  }
}




// Print method for vtkFixedPointVolumeRayCastCompositeShadeHelper
void vtkFixedPointVolumeRayCastCompositeShadeHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
