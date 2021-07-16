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

vtkStandardNewMacro(vtkFixedPointVolumeRayCastCompositeShadeHelper);

// Construct a new vtkFixedPointVolumeRayCastCompositeShadeHelper with default values
vtkFixedPointVolumeRayCastCompositeShadeHelper::vtkFixedPointVolumeRayCastCompositeShadeHelper() =
  default;

// Destruct a vtkFixedPointVolumeRayCastCompositeShadeHelper - clean up any memory used
vtkFixedPointVolumeRayCastCompositeShadeHelper::~vtkFixedPointVolumeRayCastCompositeShadeHelper() =
  default;

double getMFAValue(unsigned int *pos, float range)
{
    // cerr << "MFA=======" << endl; 
    return 1;
}

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
  /*
  Block<real_t>* b = (Block<real_t>*)(mapper->GetMfaBlock());
  VectorXi no_derivs;
  mfa::MFA_Data<real_t>&  mfa_data = *(b->vars[0].mfa_data);
  TensorProduct<real_t>&  t = mfa_data.tmesh.tensor_prods[0];
  mfa::Decoder<real_t>    decoder(mfa_data, 0);
  mfa::DecodeInfo<real_t> decode_info(mfa_data, no_derivs);
  */

  Block<real_t>* b = (Block<real_t>*)(mapper->GetMfaBlock());
  // mfa::MFA_Data<T>& mfa_data = *(b->vars[0].mfa_data);
  mfa::MFA_Data<real_t>& mfa_data = *(b->vars[0].mfa_data);
  TensorProduct<real_t>&  t = mfa_data.tmesh.tensor_prods[0];
  // mfa::Decoder<T> decoder(mfa_data, 0);
  mfa::Decoder<real_t> decoder(mfa_data, 0);
  // mfa::FastDecodeInfo<T> di(decoder);
  mfa::FastDecodeInfo<real_t> di(decoder);

  constexpr double recip = 1 / 2211772.5;

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

  // cerr << "I am MFA Data: " << mapper->GetMfaTest() << endl;
  // Block<real_t>* b = (Block<real_t>*)(mapper->GetMfaBlock());
  // cerr << "use MFA? " << mapper->GetUseMfa() << endl;

  
  /* 
  if (threadID == 0) { 
    VectorX<double> param(3); 
    VectorX<double> out_pt(4);
    param(0) = 0.5;
    param(1) = 0.5;
    param(2) = 0.5;
    b->my_decode_point(param, out_pt);
    cerr << "========thread " << threadID << "===========evaluation parameters: [" << param(0) << ", " << param(1) << ", " << param(2) << "]" << endl;
    cerr << "-decoded point: [" << out_pt[0] << ", " <<  out_pt[1] << ", " << out_pt[2] << ", " << out_pt[3] << "]" << endl;
  }
 */

#if 0
  // read in MFA model
  // diy::mpi::environment  env(argc, argv);             // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
  diy::mpi::environment  env();             // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
  diy::mpi::communicator world;                       // equivalent of MPI_COMM_WORLD
  string      infile  = "approx.out";                 // diy input file
  // initialize DIY
  diy::FileStorage storage("./DIY.XXXXXX");
  diy::Master      master(world,
          1,
          -1,
          &Block<real_t>::create,
          &Block<real_t>::destroy);
  diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()


  // Read a previously-solved MFA from disk into a DIY block
  diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
  std::cout << master.size() << " blocks read from file "<< infile << "\n\n";
#endif

  // std::ofstream myfile;
  // myfile.open("save.txt");
  // std::ofstream myval;
  // myval.open("myval.txt");
  // std::ofstream myxyz;
  // myxyz.open("myxyz.txt");
  // std::ofstream mymfaval;
  // mymfaval.open("mymfaval.txt");

  // cerr << "AAA" << threadID << endl;
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  // if (threadID == 0) {
    // cerr << "-------numSteps: " << numSteps << " count: " << count << endl;
    // count++;
    // cerr << "numstep = 0" << endl;
  // }
  // if (threadID == 2) {
  //   cerr << "+++++++numSteps: " << numSteps << " count1: " << count1 << endl;
  //   count1++;
  // }
  double maxMFAValueOnRay = 0;
  // myfile << j << "," << i << "," << numSteps << ",";
  // if (threadID == 0) { 
    // cerr << "numSteps: " << numSteps << endl;
  // }
  
  int skip_factor = 1; 
  for (k = 0; k < numSteps; k++)
  {
    if (k%skip_factor != 0) {
        continue;
    }
    /*
    if(threadID == 0) {
        if (k == numSteps -1 ){
            myfile << pos[0] <<"," << pos[1] << "," << pos[2] << "\n";
        } else {
            myfile << pos[0] <<"," << pos[1] << "," << pos[2] << ",";
        }
    }
    */
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir); // increment a dir on pos
      if (threadID == 0) {
      //   cerr << "pos " << pos[0] <<"; " << pos[1] << "; " << pos[2] << endl;
      //   cerr << "dir " << dir[0] <<"; " << dir[1] << "; " << dir[2] << endl;
      }
    }
    // VTKKWRCHelper_SpaceLeapCheck(); // contine to next iteration if check fails
    VTKKWRCHelper_CroppingCheckTrilin(pos); // contine if check cropping

    mapper->ShiftVectorDown(pos, spos); // only keep the first 15 bits of the 32bit int
    /*
    if (threadID == 0) {
        cerr << "not skiped " << endl;
        cerr << "pos " << pos[0] <<"; " << pos[1] << "; " << pos[2] << endl;
        cerr << "spos " << spos[0] <<"; " << spos[1] << "; " << spos[2] << endl;
        cerr << "data" << typeid(data).name() << endl;
        cerr << "data print" << static_cast<unsigned int>(*(data)) << endl;
    }
    */
    // start ray traversal
    // if (!mapper->GetUseMfa()) {
       //  before = std::chrono::system_clock::now();
        if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
        {
          oldSPos[0] = spos[0];
          oldSPos[1] = spos[1];
          oldSPos[2] = spos[2];

          dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
          /*
          if (threadID == 0 ) { 
              cerr << "inc: " << inc[0] << "; " << inc[1] << "; " << inc[2] << endl;
              cerr << "dptr" << static_cast<unsigned int>(*(dptr)) << endl;
              printf("data u: %u \n", static_cast<unsigned int>(*(data)));
              printf("dptr u: %u \n", static_cast<unsigned int>(*(dptr)));
          }
          */
          VTKKWRCHelper_GetCellScalarValuesSimple(dptr); // update 8 corners of the cell
          dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
          dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
          needToSampleDirection = 1;
        }

    if (!mapper->GetUseMfa()) {
        before = std::chrono::system_clock::now();
        VTKKWRCHelper_ComputeWeights(pos);
        /*
        if (threadID == 0) {
            printf("val before %u\n", (unsigned int)val);
        }
        */
        VTKKWRCHelper_InterpolateScalar(val); // val is unsigned short, update val, val is index

        after = std::chrono::system_clock::now();
        // auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(after - before);
        if (threadID == 0) {
            time0 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 1) {
            time1 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 2) {
            time2 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 3) {
            time3 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 4) {
            time4 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 5) {
            time5 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 6) {
            time6 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        if (threadID == 7) {
            time7 += nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << nanoseconds.count() << endl;
        }
        // Replace val by value retrieved by querying MFA model, val in range 0~255
    } else {
        mfa_before = std::chrono::system_clock::now();
        // get position x, y, z in range [0, 1]
        /*
        double x = (double)pos[0] / 2211772.5;
        double y = (double)pos[1] / 2211772.5;
        double z = (double)pos[2] / 2211772.5;
        param(0) = x;
        param(1) = y;
        param(2) = z;
        */
        param(0) = pos[0] * recip;
        param(1) = pos[1] * recip;
        param(2) = pos[2] * recip;
        // do MFA decoding
        // b->my_decode_point(param, out_pt4);
        // decoder.VolPt(param, out_pt, decode_info, t);
        decoder.FastVolPt(param, out_pt, di, t);
        // cerr << "========thread " << threadID << "===========evaluation parameters: [" << param(0) << ", " << param(1) << ", " << param(2) << "]" << endl;
        // cerr << "-decoded point: [" << out_pt[0] << ", " <<  out_pt[1] << ", " << out_pt[2] << ", " << out_pt[3] << "]" << endl;
        // mymfaval << out_pt[0] << "\n";
        // if (threadID == 0)
        //     cerr << "+++++++++++mymfaval: " << out_pt[0] << "," << out_pt4[3] << "\n";

        double v = out_pt[0];
        double ratio = (v + 2.0)/11.0; // for sinc dataset
        // double ratio = v/136.0; // for nek5000 dataset
        if (ratio > 1) {
            ratio = 1;
        }
        if (ratio < 0) {
            ratio = 0;
        }
        unsigned short new_val = (unsigned short)(255*ratio);
        val = new_val;
        mfa_after = std::chrono::system_clock::now();
        auto mfa_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(mfa_after - mfa_before);
        if (threadID == 0) {
            mfa_time0 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter0++;
        }
        if (threadID == 1) {
            mfa_time1 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter1++;
        }
        if (threadID == 2) {
            mfa_time2 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter2++;
        }
        if (threadID == 3) {
            mfa_time3 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter3++;
        }
        if (threadID == 4) {
            mfa_time4 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter4++;
        }
        if (threadID == 5) {
            mfa_time5 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter5++;
        }
        if (threadID == 6) {
            mfa_time6 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter6++;
        }
        if (threadID == 7) {
            mfa_time7 += mfa_nanoseconds.count();
            // cerr << "time0: " << time0 << endl;
            // cerr << "time0: " << mfa_nanoseconds.count() << endl;
            mfa_counter7++;
        }
    }

    // cerr << "new val: " << val << endl;
    
    /*
    if (threadID == 0) {
        printf("val after %u\n", (unsigned int)val);
        if (val > maxx) { maxx = val; }
        if (val < minn) { minn = val; }
        cerr << "min: " << minn << "; max: " << maxx << endl;
    }
    */
    // if (threadID == 0) {val = 1;} 
    // val = 40; // range 0~255 higher -> darker in final image
    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp); // update tmp, tmp[0,1,2] unsigned short
    // if (threadID == 0) { cerr << "1 tmp[0, 1, 2]: " << tmp[0] << "; " << tmp[1] << "; " << tmp[2] << endl; }
    // if (threadID == 0) { cerr << "needToSampleDirection" << needToSampleDirection << endl; }
    if (needToSampleDirection) // true
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }

    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp); // update tmp again
    // if (threadID == 0) { cerr << "2 tmp[0, 1, 2]: " << tmp[0] << "; " << tmp[1] << "; " << tmp[2] << endl; }
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity); // break for loop if terminates
    // if (threadID == 0) { cerr << "color[0, 1, 2]: " << color[0] << "; " << color[1] << "; " << color[2] << endl; }
    //
    //
    unsigned short* imagePtr;                                                                        \
    double currMFAValue;
    currMFAValue = getMFAValue(pos, 68.0);
    if (currMFAValue > maxMFAValueOnRay) {
        maxMFAValueOnRay = currMFAValue;
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  if (threadID == 0) {
      // cerr << "imagePtr[0, 1, 2, 3]: " << imagePtr[0] << "; " << imagePtr[1] << "; " << imagePtr[2] << "; " << imagePtr[3] << endl; // 0~32767
      // imagePtr[0] = 32767; // R
      // imagePtr[1] = 0; // G
      // imagePtr[2] = 0; // B
      // imagePtr[3] = 32767; // Alpha
  }
  VTKKWRCHelper_IncrementAndLoopEnd();
  // myfile.close();
  // myval.close();
  // myxyz.close();
  // mymfaval.close();
  
  if (threadID == 0) {
      // std::ofstream tffile;
      // tffile.open("tffile.txt");
      // for (int cc = 0; cc < 32768; cc++) {
      //   tffile << scalarOpacityTable[0][cc] << "," ;
      // }
      // tffile.close();
  }

  if (threadID == 0) {
      std::ofstream filetime0;
      filetime0.open("time0.txt");
      filetime0 << time0;
      filetime0.close();
      // cerr << "time0: " << time0 << endl;

      std::ofstream mfa_filetime0;
      mfa_filetime0.open("mfa_time0.txt");
      mfa_filetime0 << mfa_time0;
      mfa_filetime0.close();
      // cerr << "mfa_time0: " << mfa_time0 << endl;

      std::ofstream mfa_filecount0;
      mfa_filecount0.open("mfa_count0.txt");
      mfa_filecount0 << mfa_counter0;
      mfa_filecount0.close();
      // cerr << "mfa_count0: " << mfa_counter0 << endl;

  }
  if (threadID == 1) {
      std::ofstream filetime1;
      filetime1.open("time1.txt");
      filetime1 << time1;
      filetime1.close();
      // cerr << "time1: " << time1 << endl;

      std::ofstream mfa_filetime1;
      mfa_filetime1.open("mfa_time1.txt");
      mfa_filetime1 << mfa_time1;
      mfa_filetime1.close();
      // cerr << "mfa_time1: " << mfa_time1 << endl;

      std::ofstream mfa_filecount1;
      mfa_filecount1.open("mfa_count1.txt");
      mfa_filecount1 << mfa_counter1;
      mfa_filecount1.close();
      // cerr << "mfa_count1: " << mfa_counter1 << endl;
  }
  if (threadID == 2) {
      std::ofstream filetime2;
      filetime2.open("time2.txt");
      filetime2 << time2;
      filetime2.close();
      // cerr << "time2: " << time2 << endl;

      std::ofstream mfa_filetime2;
      mfa_filetime2.open("mfa_time2.txt");
      mfa_filetime2 << mfa_time2;
      mfa_filetime2.close();
      // cerr << "mfa_time2: " << mfa_time2 << endl;

      std::ofstream mfa_filecount2;
      mfa_filecount2.open("mfa_count2.txt");
      mfa_filecount2 << mfa_counter2;
      mfa_filecount2.close();
      // cerr << "mfa_count2: " << mfa_counter2 << endl;
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
      // cerr << "time4: " << time4 << endl;

      std::ofstream mfa_filetime4;
      mfa_filetime4.open("mfa_time4.txt");
      mfa_filetime4 << mfa_time4;
      mfa_filetime4.close();
      // cerr << "mfa_time4: " << mfa_time4 << endl;

      std::ofstream mfa_filecount4;
      mfa_filecount4.open("mfa_count4.txt");
      mfa_filecount4 << mfa_counter4;
      mfa_filecount4.close();
      // cerr << "mfa_count4: " << mfa_counter4 << endl;
  }
  if (threadID == 5) {
      std::ofstream filetime5;
      filetime5.open("time5.txt");
      filetime5 << time5;
      filetime5.close();
      // cerr << "time5: " << time5 << endl;

      std::ofstream mfa_filetime5;
      mfa_filetime5.open("mfa_time5.txt");
      mfa_filetime5 << mfa_time5;
      mfa_filetime5.close();
      // cerr << "mfa_time5: " << mfa_time5 << endl;

      std::ofstream mfa_filecount5;
      mfa_filecount5.open("mfa_count5.txt");
      mfa_filecount5 << mfa_counter5;
      mfa_filecount5.close();
      // cerr << "mfa_count5: " << mfa_counter5 << endl;
  }
  if (threadID == 6) {
      std::ofstream filetime6;
      filetime6.open("time6.txt");
      filetime6 << time6;
      filetime6.close();
      // cerr << "time6: " << time6 << endl;

      std::ofstream mfa_filetime6;
      mfa_filetime6.open("mfa_time6.txt");
      mfa_filetime6 << mfa_time6;
      mfa_filetime6.close();
      // cerr << "mfa_time6: " << mfa_time6 << endl;

      std::ofstream mfa_filecount6;
      mfa_filecount6.open("mfa_count6.txt");
      mfa_filecount6 << mfa_counter6;
      mfa_filecount6.close();
      // cerr << "mfa_count6: " << mfa_counter6 << endl;
  }
  if (threadID == 7) {
      std::ofstream filetime7;
      filetime7.open("time7.txt");
      filetime7 << time7;
      filetime7.close();
      // cerr << "time7: " << time7 << endl;

      std::ofstream mfa_filetime7;
      mfa_filetime7.open("mfa_time7.txt");
      mfa_filetime7 << mfa_time7;
      mfa_filetime7.close();
      // cerr << "mfa_time7: " << mfa_time7 << endl;

      std::ofstream mfa_filecount7;
      mfa_filecount7.open("mfa_count7.txt");
      mfa_filecount7 << mfa_counter7;
      mfa_filecount7.close();
      // cerr << "mfa_count7: " << mfa_counter7 << endl;
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
