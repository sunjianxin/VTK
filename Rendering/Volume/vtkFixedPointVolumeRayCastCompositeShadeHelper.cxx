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

#define RECORD_TIMING

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

void getDiffuseSpecularCoefficients(double viewDirection[3],
                                    double lightDirection[3],
                                    double lightAmbientColor[3],
                                    double lightDiffuseColor[3],
                                    double lightSpecularColor[3],
                                    double lightIntensity,
                                    double material[4],
                                    double gx,
                                    double gy,
                                    double gz,
                                    // double* diffuseCoefficients,
                                    // double* specularCoefficients)
                                    unsigned short* diffuseCoefficients,
                                    unsigned short* specularCoefficients)

{
  double half_x, half_y, half_z;
  double g_dot_l, g_dot_h, g_dot_v;
  double Ka, Es, Kd_intensity, Ks_intensity;
  double mag, specular_value;
  float dr, dg, db, sr, sg, sb;
  
  // cerr << "out estimator: _viewDirection: " << viewDirection[0] << ", " << viewDirection[1] << ", " << viewDirection[2] << endl;
  // cerr << "out estimator: _lightDirection: " << lightDirection[0] << ", " << lightDirection[1] << ", " << lightDirection[2] << endl;
  // cerr << "out estimator: _lightAmbientColor " << lightAmbientColor[0] << ", " << lightAmbientColor[1] << ", " << lightAmbientColor[2] << endl;
  // cerr << "out estimator: _lightDiffuseColor " << lightDiffuseColor[0] << ", " << lightDiffuseColor[1] << ", " << lightDiffuseColor[2] << endl;
  // cerr << "out estimator: _lightDiffuseColor " << lightSpecularColor[0] << ", " << lightSpecularColor[1] << ", " << lightSpecularColor[2] << endl;
  // cerr << "out estimator: _lightIntensity: " << lightIntensity << endl;
  // cerr << "out estimator: _material " << material[0] << ", " << material[1] << ", " << material[2] << ", " << material[3] << endl;

  half_x = lightDirection[0] - viewDirection[0];
  half_y = lightDirection[1] - viewDirection[1];
  half_z = lightDirection[2] - viewDirection[2];

  mag = sqrt(static_cast<double>(half_x * half_x + half_y * half_y + half_z * half_z));

  if (mag != 0.0)
  {
    half_x /= mag;
    half_y /= mag;
    half_z /= mag;
  }

  Ka = material[0] * lightIntensity;
  Es = material[3];
  Kd_intensity = material[1] * lightIntensity;
  Ks_intensity = material[2] * lightIntensity;
      
  g_dot_l = gx * lightDirection[0] + gy * lightDirection[1] + gz * lightDirection[2];
  g_dot_h = gx * half_x + gy * half_y + gz * half_z;
  g_dot_v = gx * viewDirection[0] + gy * viewDirection[1] + gz * viewDirection[2];
  if (g_dot_v > 0.0) {
    g_dot_l = -g_dot_l;
    g_dot_h = -g_dot_h;
  }

  dr = static_cast<float>(Ka*lightAmbientColor[0]);
  dg = static_cast<float>(Ka*lightAmbientColor[1]);
  db = static_cast<float>(Ka*lightAmbientColor[2]);
  sr = 0;
  sg = 0;
  sb = 0;
  if (g_dot_l > 0)
  {
    dr += static_cast<float>(Kd_intensity * g_dot_l * lightDiffuseColor[0]);
    dg += static_cast<float>(Kd_intensity * g_dot_l * lightDiffuseColor[1]);
    db += static_cast<float>(Kd_intensity * g_dot_l * lightDiffuseColor[2]);

    if (g_dot_h > 0.001)
    {
      specular_value =
        Ks_intensity * pow(static_cast<double>(g_dot_h), static_cast<double>(Es));
      sr += static_cast<float>(specular_value * lightSpecularColor[0]);
      sg += static_cast<float>(specular_value * lightSpecularColor[1]);
      sb += static_cast<float>(specular_value * lightSpecularColor[2]);
    }
  }

  // diffuseCoefficients[0] = dr;
  // diffuseCoefficients[1] = dg;
  // diffuseCoefficients[2] = db;
  // specularCoefficients[0] = sr;
  // specularCoefficients[1] = sg;
  // specularCoefficients[2] = sb;
    
  
  diffuseCoefficients[0] = static_cast<unsigned short>(dr * VTKKW_FP_SCALE + 0.5);
  diffuseCoefficients[1] = static_cast<unsigned short>(dg * VTKKW_FP_SCALE + 0.5);
  diffuseCoefficients[2] = static_cast<unsigned short>(db * VTKKW_FP_SCALE + 0.5);
  specularCoefficients[0] = static_cast<unsigned short>(sr * VTKKW_FP_SCALE + 0.5);
  specularCoefficients[1] = static_cast<unsigned short>(sg * VTKKW_FP_SCALE + 0.5);
  specularCoefficients[2] = static_cast<unsigned short>(sb * VTKKW_FP_SCALE + 0.5);
}

long int value_dis_time0 = 0;
long int value_dis_time1 = 0;
long int value_dis_time2 = 0;
long int value_dis_time3 = 0;
long int value_dis_time4 = 0;
long int value_dis_time5 = 0;
long int value_dis_time6 = 0;
long int value_dis_time7 = 0;
long int shade_dis_time0 = 0;
long int shade_dis_time1 = 0;
long int shade_dis_time2 = 0;
long int shade_dis_time3 = 0;
long int shade_dis_time4 = 0;
long int shade_dis_time5 = 0;
long int shade_dis_time6 = 0;
long int shade_dis_time7 = 0;
long int dis_counter0 = 0;
long int dis_counter1 = 0;
long int dis_counter2 = 0;
long int dis_counter3 = 0;
long int dis_counter4 = 0;
long int dis_counter5 = 0;
long int dis_counter6 = 0;
long int dis_counter7 = 0;
long int _dis_counter0 = 0;
long int _dis_counter1 = 0;
long int _dis_counter2 = 0;
long int _dis_counter3 = 0;
long int _dis_counter4 = 0;
long int _dis_counter5 = 0;
long int _dis_counter6 = 0;
long int _dis_counter7 = 0;
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
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  /* Time-stamp initialization */
  std::chrono::time_point<std::chrono::system_clock> before;
  std::chrono::time_point<std::chrono::system_clock> after;
#if defined(RECORD_TIMING) // record timing
  if (threadID == 0) {
    value_dis_time0 = 0;
    shade_dis_time0 = 0;
    dis_counter0 = 0;
    _dis_counter0 = 0;
  }
  if (threadID == 1) {
    value_dis_time1 = 0;
    shade_dis_time1 = 0;
    dis_counter1 = 0;
    _dis_counter1 = 0;
  }
  if (threadID == 2) {
    value_dis_time2 = 0;
    shade_dis_time2 = 0;
    dis_counter2 = 0;
    _dis_counter2 = 0;
  }
  if (threadID == 3) {
    value_dis_time3 = 0;
    shade_dis_time3 = 0;
    dis_counter3 = 0;
    _dis_counter3 = 0;
  }
  if (threadID == 4) {
    value_dis_time4 = 0;
    shade_dis_time4 = 0;
    dis_counter4 = 0;
    _dis_counter4 = 0;
  }
  if (threadID == 5) {
    value_dis_time5 = 0;
    shade_dis_time5 = 0;
    dis_counter5 = 0;
    _dis_counter5 = 0;
  }
  if (threadID == 6) {
    value_dis_time6 = 0;
    shade_dis_time6 = 0;
    dis_counter6 = 0;
    _dis_counter6 = 0;
  }
  if (threadID == 7) {
    value_dis_time7 = 0;
    shade_dis_time7 = 0;
    dis_counter7 = 0;
    _dis_counter7 = 0;
  }
#endif

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
      VTKKWRCHelper_GetCellScalarValuesSimple(dptr);
      dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1];
      dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1];
      needToSampleDirection = 1;
    }

    before = std::chrono::system_clock::now();
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);
    after = std::chrono::system_clock::now();
    auto value_dis_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(after - before);
#if defined(RECORD_TIMING) // record timing
    if (threadID == 0) {
        value_dis_time0 += value_dis_nanoseconds.count();
        dis_counter0++;
    }
    if (threadID == 1) {
        value_dis_time1 += value_dis_nanoseconds.count();
        dis_counter1++;
    }
    if (threadID == 2) {
        value_dis_time2 += value_dis_nanoseconds.count();
        dis_counter2++;
    }
    if (threadID == 3) {
        value_dis_time3 += value_dis_nanoseconds.count();
        dis_counter3++;
    }
    if (threadID == 4) {
        value_dis_time4 += value_dis_nanoseconds.count();
        dis_counter4++;
    }
    if (threadID == 5) {
        value_dis_time5 += value_dis_nanoseconds.count();
        dis_counter5++;
    }
    if (threadID == 6) {
        value_dis_time6 += value_dis_nanoseconds.count();
        dis_counter6++;
    }
    if (threadID == 7) {
        value_dis_time7 += value_dis_nanoseconds.count();
        dis_counter7++;
    }
#endif


    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp);
    
    before = std::chrono::system_clock::now();
    if (needToSampleDirection)
    {
      VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH);
      needToSampleDirection = 0;
    }

    VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp);
    after = std::chrono::system_clock::now();
    auto shade_dis_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(after - before);
#if defined(RECORD_TIMING) // record timing
    if (threadID == 0) {
        shade_dis_time0 += shade_dis_nanoseconds.count();
        _dis_counter0++;
    }
    if (threadID == 1) {
        shade_dis_time1 += shade_dis_nanoseconds.count();
        _dis_counter1++;
    }
    if (threadID == 2) {
        shade_dis_time2 += shade_dis_nanoseconds.count();
        _dis_counter2++;
    }
    if (threadID == 3) {
        shade_dis_time3 += shade_dis_nanoseconds.count();
        _dis_counter3++;
    }
    if (threadID == 4) {
        shade_dis_time4 += shade_dis_nanoseconds.count();
        _dis_counter4++;
    }
    if (threadID == 5) {
        shade_dis_time5 += shade_dis_nanoseconds.count();
        _dis_counter5++;
    }
    if (threadID == 6) {
        shade_dis_time6 += shade_dis_nanoseconds.count();
        _dis_counter6++;
    }
    if (threadID == 7) {
        shade_dis_time7 += shade_dis_nanoseconds.count();
        _dis_counter7++;
    }
#endif
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();

#if defined(RECORD_TIMING) // record timing
  /* Save Time-stamp */ 
  if (threadID == 0) {
      std::ofstream value_dis_time0_file;
      value_dis_time0_file.open("shade/value_dis_time0.txt");
      value_dis_time0_file << value_dis_time0;
      value_dis_time0_file.close();

      std::ofstream shade_dis_time0_file;
      shade_dis_time0_file.open("shade/shade_dis_time0.txt");
      shade_dis_time0_file << shade_dis_time0;
      shade_dis_time0_file.close();

      std::ofstream dis_count0_file;
      dis_count0_file.open("shade/dis_count0.txt");
      dis_count0_file << dis_counter0;
      dis_count0_file.close();

      std::ofstream _dis_count0_file;
      _dis_count0_file.open("shade/_dis_count0.txt");
      _dis_count0_file << _dis_counter0;
      _dis_count0_file.close();
  }
  if (threadID == 1) {
      std::ofstream value_dis_time1_file;
      value_dis_time1_file.open("shade/value_dis_time1.txt");
      value_dis_time1_file << value_dis_time1;
      value_dis_time1_file.close();

      std::ofstream shade_dis_time1_file;
      shade_dis_time1_file.open("shade/shade_dis_time1.txt");
      shade_dis_time1_file << shade_dis_time1;
      shade_dis_time1_file.close();

      std::ofstream dis_count1_file;
      dis_count1_file.open("shade/dis_count1.txt");
      dis_count1_file << dis_counter1;
      dis_count1_file.close();

      std::ofstream _dis_count1_file;
      _dis_count1_file.open("shade/_dis_count1.txt");
      _dis_count1_file << _dis_counter1;
      _dis_count1_file.close();
  }
  if (threadID == 2) {
      std::ofstream value_dis_time2_file;
      value_dis_time2_file.open("shade/value_dis_time2.txt");
      value_dis_time2_file << value_dis_time2;
      value_dis_time2_file.close();

      std::ofstream shade_dis_time2_file;
      shade_dis_time2_file.open("shade/shade_dis_time2.txt");
      shade_dis_time2_file << shade_dis_time2;
      shade_dis_time2_file.close();

      std::ofstream dis_count2_file;
      dis_count2_file.open("shade/dis_count2.txt");
      dis_count2_file << dis_counter2;
      dis_count2_file.close();

      std::ofstream _dis_count2_file;
      _dis_count2_file.open("shade/_dis_count2.txt");
      _dis_count2_file << _dis_counter2;
      _dis_count2_file.close();
  }
  if (threadID == 3) {
      std::ofstream value_dis_time3_file;
      value_dis_time3_file.open("shade/value_dis_time3.txt");
      value_dis_time3_file << value_dis_time3;
      value_dis_time3_file.close();

      std::ofstream shade_dis_time3_file;
      shade_dis_time3_file.open("shade/shade_dis_time3.txt");
      shade_dis_time3_file << shade_dis_time3;
      shade_dis_time3_file.close();

      std::ofstream dis_count3_file;
      dis_count3_file.open("shade/dis_count3.txt");
      dis_count3_file << dis_counter3;
      dis_count3_file.close();

      std::ofstream _dis_count3_file;
      _dis_count3_file.open("shade/_dis_count3.txt");
      _dis_count3_file << _dis_counter3;
      _dis_count3_file.close();
  }
  if (threadID == 4) {
      std::ofstream value_dis_time4_file;
      value_dis_time4_file.open("shade/value_dis_time4.txt");
      value_dis_time4_file << value_dis_time4;
      value_dis_time4_file.close();

      std::ofstream shade_dis_time4_file;
      shade_dis_time4_file.open("shade/shade_dis_time4.txt");
      shade_dis_time4_file << shade_dis_time4;
      shade_dis_time4_file.close();

      std::ofstream dis_count4_file;
      dis_count4_file.open("shade/dis_count4.txt");
      dis_count4_file << dis_counter4;
      dis_count4_file.close();

      std::ofstream _dis_count4_file;
      _dis_count4_file.open("shade/_dis_count4.txt");
      _dis_count4_file << _dis_counter4;
      _dis_count4_file.close();
  }
  if (threadID == 5) {
      std::ofstream value_dis_time5_file;
      value_dis_time5_file.open("shade/value_dis_time5.txt");
      value_dis_time5_file << value_dis_time5;
      value_dis_time5_file.close();

      std::ofstream shade_dis_time5_file;
      shade_dis_time5_file.open("shade/shade_dis_time5.txt");
      shade_dis_time5_file << shade_dis_time5;
      shade_dis_time5_file.close();

      std::ofstream dis_count5_file;
      dis_count5_file.open("shade/dis_count5.txt");
      dis_count5_file << dis_counter5;
      dis_count5_file.close();

      std::ofstream _dis_count5_file;
      _dis_count5_file.open("shade/_dis_count5.txt");
      _dis_count5_file << _dis_counter5;
      _dis_count5_file.close();
  }
  if (threadID == 6) {
      std::ofstream value_dis_time6_file;
      value_dis_time6_file.open("shade/value_dis_time6.txt");
      value_dis_time6_file << value_dis_time6;
      value_dis_time6_file.close();

      std::ofstream shade_dis_time6_file;
      shade_dis_time6_file.open("shade/shade_dis_time6.txt");
      shade_dis_time6_file << shade_dis_time6;
      shade_dis_time6_file.close();

      std::ofstream dis_count6_file;
      dis_count6_file.open("shade/dis_count6.txt");
      dis_count6_file << dis_counter6;
      dis_count6_file.close();

      std::ofstream _dis_count6_file;
      _dis_count6_file.open("shade/_dis_count6.txt");
      _dis_count6_file << _dis_counter6;
      _dis_count6_file.close();
  }
  if (threadID == 7) {
      std::ofstream value_dis_time7_file;
      value_dis_time7_file.open("shade/value_dis_time7.txt");
      value_dis_time7_file << value_dis_time7;
      value_dis_time7_file.close();

      std::ofstream shade_dis_time7_file;
      shade_dis_time7_file.open("shade/shade_dis_time7.txt");
      shade_dis_time7_file << shade_dis_time7;
      shade_dis_time7_file.close();

      std::ofstream dis_count7_file;
      dis_count7_file.open("shade/dis_count7.txt");
      dis_count7_file << dis_counter7;
      dis_count7_file.close();

      std::ofstream _dis_count7_file;
      _dis_count7_file.open("shade/_dis_count7.txt");
      _dis_count7_file << _dis_counter7;
      _dis_count7_file.close();
  }
#endif
}

long int value_mfa_time0 = 0;
long int value_mfa_time1 = 0;
long int value_mfa_time2 = 0;
long int value_mfa_time3 = 0;
long int value_mfa_time4 = 0;
long int value_mfa_time5 = 0;
long int value_mfa_time6 = 0;
long int value_mfa_time7 = 0;
long int shade_mfa_time0 = 0;
long int shade_mfa_time1 = 0;
long int shade_mfa_time2 = 0;
long int shade_mfa_time3 = 0;
long int shade_mfa_time4 = 0;
long int shade_mfa_time5 = 0;
long int shade_mfa_time6 = 0;
long int shade_mfa_time7 = 0;
long int shade_mfa_counter0 = 0;
long int shade_mfa_counter1 = 0;
long int shade_mfa_counter2 = 0;
long int shade_mfa_counter3 = 0;
long int shade_mfa_counter4 = 0;
long int shade_mfa_counter5 = 0;
long int shade_mfa_counter6 = 0;
long int shade_mfa_counter7 = 0;

long int _shade_mfa_counter0 = 0;
long int _shade_mfa_counter1 = 0;
long int _shade_mfa_counter2 = 0;
long int _shade_mfa_counter3 = 0;
long int _shade_mfa_counter4 = 0;
long int _shade_mfa_counter5 = 0;
long int _shade_mfa_counter6 = 0;
long int _shade_mfa_counter7 = 0;

template <class T>
void vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleTrilinMfa(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  // cerr << "AAA" << threadID << endl;
  /* Retrieve MFA reference pointer */ 
  Block<real_t>* b = (Block<real_t>*)(mapper->GetMfaBlock());
  mfa::MFA_Data<real_t>& mfa_data = *(b->vars[0].mfa_data);
  TensorProduct<real_t>&  t = mfa_data.tmesh.tensor_prods[0];
  mfa::Decoder<real_t> decoder(mfa_data, 0);
  mfa::FastDecodeInfo<real_t> di(decoder);
  di.ResizeDers(1);
  mfa::DecodeInfo<real_t> di2(mfa_data, VectorXi::Ones(3));
  /* Retrieve dataset size and update scalar for POS */
  int size =  mapper->GetMfaSize() - 1;
  std::string dataset =  mapper->GetDataset();
  double recip = 1 / ((size*32767.0) + 0.5); // 1 / (size*VTKKW_FP_SCALE + 0.5)
  if (threadID == 0) {
    cerr << "size: " << size << endl;
    cerr << "dataset: " << dataset << endl;
  }
  
  /* Place holder for getting value */
  VectorX<double> param(3); 
  VectorX<double> out_pt(1);
  VectorX<double> out_pt4(4);
  /* Place holder for getting derivative */
  real_t          pt_value = 0;
  VectorX<real_t> gradient(3);
  VectorX<real_t> out_dx(4);       // will store x-derivative at the point
  VectorX<real_t> out_dy(4);       // will store y-deriv…
  VectorX<real_t> out_dz(4);       // will store z-deriv… 
  VectorX<real_t> extents = b->bounds_maxs - b->bounds_mins;                 // for rescaling the derivative
  real_t          value_extent_recip = 1 / (b->bounds_maxs(3) - b->bounds_mins(3)); // for calculating intensity      
  /* Place holder for getting diffuse and specular coefficients */
  unsigned short diffuseCoefficients[3];
  unsigned short specularCoefficients[3];

  /* Time-stamp initialization */
  std::chrono::time_point<std::chrono::system_clock> mfa_before;
  std::chrono::time_point<std::chrono::system_clock> mfa_after;
#if defined(RECORD_TIMING) // record timing
  if (threadID == 0) {
    value_mfa_time0 = 0;
    shade_mfa_time0 = 0;
    shade_mfa_counter0 = 0;
    _shade_mfa_counter0 = 0;
  }
  if (threadID == 1) {
    value_mfa_time1 = 0;
    shade_mfa_time1 = 0;
    shade_mfa_counter1 = 0;
    _shade_mfa_counter1 = 0;
  }
  if (threadID == 2) {
    value_mfa_time2 = 0;
    shade_mfa_time2 = 0;
    shade_mfa_counter2 = 0;
    _shade_mfa_counter2 = 0;
  }
  if (threadID == 3) {
    value_mfa_time3 = 0;
    shade_mfa_time3 = 0;
    shade_mfa_counter3 = 0;
    _shade_mfa_counter3 = 0;
  }
  if (threadID == 4) {
    value_mfa_time4 = 0;
    shade_mfa_time4 = 0;
    shade_mfa_counter4 = 0;
    _shade_mfa_counter4 = 0;
  }
  if (threadID == 5) {
    value_mfa_time5 = 0;
    shade_mfa_time5 = 0;
    shade_mfa_counter5 = 0;
    _shade_mfa_counter5 = 0;
  }
  if (threadID == 6) {
    value_mfa_time6 = 0;
    shade_mfa_time6 = 0;
    shade_mfa_counter6 = 0;
    _shade_mfa_counter6 = 0;
  }
  if (threadID == 7) {
    value_mfa_time7 = 0;
    shade_mfa_time7 = 0;
    shade_mfa_counter7 = 0;
    _shade_mfa_counter7 = 0;
  }
#endif


  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  // VTKKWRCHelper_InitializationAndLoopStartTrilin();
  // VTKKWRCHelper_InitializeCompositeOneTrilin();
  // VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0; //
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
      VTKKWRCHelper_GetCellScalarValuesSimple(dptr);
      // dirPtrABCD = gradientDir[spos[2]] + spos[0] * dInc[0] + spos[1] * dInc[1]; // get base address for ABCD
      // dirPtrEFGH = gradientDir[spos[2] + 1] + spos[0] * dInc[0] + spos[1] * dInc[1]; // get base address for EFGH
      // needToSampleDirection = 1; //
    }

    // VTKKWRCHelper_ComputeWeights(pos);
    // VTKKWRCHelper_InterpolateScalar(val); // step of value query: dis or mfa

    mfa_before = std::chrono::system_clock::now();
    // VTKKWRCHelper_ComputeWeights(pos);
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
    decoder.FastGrad(param, di, t, gradient, &pt_value);
    // decoder.FastVolPt(param, out_pt, di, t);
    // b->my_decode_point(param, out_pt4);
    if (threadID == 0) {
        // cerr << "after decoder" << endl;
    }
    /* Normalize value */
    // double v = out_pt[0];
    // double v = out_pt[3];
    double ratio;
    double bounds_lower = b->bounds_mins[3];
    double bounds_upper = b->bounds_maxs[3];
    ratio = (pt_value - bounds_lower) * value_extent_recip;
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
    auto value_mfa_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(mfa_after - mfa_before);
#if defined(RECORD_TIMING) // record timing
    if (threadID == 0) {
        value_mfa_time0 += value_mfa_nanoseconds.count();
        shade_mfa_counter0++;
    }
    if (threadID == 1) {
        value_mfa_time1 += value_mfa_nanoseconds.count();
        shade_mfa_counter1++;
    }
    if (threadID == 2) {
        value_mfa_time2 += value_mfa_nanoseconds.count();
        shade_mfa_counter2++;
    }
    if (threadID == 3) {
        value_mfa_time3 += value_mfa_nanoseconds.count();
        shade_mfa_counter3++;
    }
    if (threadID == 4) {
        value_mfa_time4 += value_mfa_nanoseconds.count();
        shade_mfa_counter4++;
    }
    if (threadID == 5) {
        value_mfa_time5 += value_mfa_nanoseconds.count();
        shade_mfa_counter5++;
    }
    if (threadID == 6) {
        value_mfa_time6 += value_mfa_nanoseconds.count();
        shade_mfa_counter6++;
    }
    if (threadID == 7) {
        value_mfa_time7 += value_mfa_nanoseconds.count();
        shade_mfa_counter7++;
    }
#endif
    VTKKWRCHelper_LookupColorUS(colorTable[0], scalarOpacityTable[0], val, tmp);

    //  if (needToSampleDirection) // 
    // {
    //   VTKKWRCHelper_GetCellDirectionValues(dirPtrABCD, dirPtrEFGH); // get normal for A, B, C, D, E, F, G, H
    //   needToSampleDirection = 0; //
    // }

    // VTKKWRCHelper_InterpolateShading(diffuseShadingTable[0], specularShadingTable[0], tmp); //diffuseShadingTable[0]/specularShadingTable[0] unsigned short pointer, starting address of the shading tables, range:[0, 65535]
#if 1
    mfa_before = std::chrono::system_clock::now();
    /* getting normal from MFA */
    // b->differentiate_point_simple(param, 1, 0, -1, out_dx); // first derivative, direction 0
    // b->differentiate_point_simple(param, 1, 1, -1, out_dy); // first derivative, direction 1
    // b->differentiate_point_simple(param, 1, 2, -1, out_dz); // first derivative, direction 2
    // out_dx(3) *= extents(0);
    // out_dy(3) *= extents(1);
    // out_dz(3) *= extents(2);
    gradient(0) *= extents(0);
    gradient(1) *= extents(1);
    gradient(2) *= extents(2);
    // normalize surface normal
    gradient.normalize();
    // double norm = sqrt(pow(out_dx(3), 2.0) + pow(out_dy(3), 2.0) + pow(out_dz(3), 2.0));
    // out_dx(3) /= norm;
    // out_dy(3) /= norm;
    // out_dz(3) /= norm;

    /* getting view direction */
    // ren->GetActiveCamera()->GetPosition(cameraPosition);
    if (threadID == 0) {
      // cerr << "in mapper: _viewDirection: " << mapper->GetViewDirection()[0] << ", " << mapper->GetViewDirection()[1] << ", " << mapper->GetViewDirection()[2] << endl;
      // cerr << "in mapper: _lightDirection: " << mapper->GetLightDirection()[0] << ", " << mapper->GetLightDirection()[1] << ", " << mapper->GetLightDirection()[2] << endl;
      // cerr << "in mapper: _lightAmbientColor " << mapper->GetLightAmbientColor()[0] << ", " << mapper->GetLightAmbientColor()[1] << ", " << mapper->GetLightAmbientColor()[2] << endl;
      // cerr << "in mapper: _lightDiffuseColor " << mapper->GetLightDiffuseColor()[0] << ", " << mapper->GetLightDiffuseColor()[1] << ", " << mapper->GetLightDiffuseColor()[2] << endl;
      // cerr << "in mapper: _lightDiffuseColor " << mapper->GetLightSpecularColor()[0] << ", " << mapper->GetLightSpecularColor()[1] << ", " << mapper->GetLightSpecularColor()[2] << endl;
      // cerr << "in mapper: _lightIntensity: " << mapper->GetLightIntensity() << endl;
      // cerr << "in mapper: _material " << mapper->GetMaterial()[0] << ", " << mapper->GetMaterial()[1] << ", " << mapper->GetMaterial()[2] << ", " << mapper->GetMaterial()[3] << endl;
      // cerr << "normal: " << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << endl;
    }
    getDiffuseSpecularCoefficients(mapper->GetViewDirection(),
                                   mapper->GetLightDirection(),
                                   mapper->GetLightAmbientColor(),
                                   mapper->GetLightDiffuseColor(),
                                   mapper->GetLightSpecularColor(),
                                   mapper->GetLightIntensity(),
                                   mapper->GetMaterial(),
                                   gradient(0),
                                   gradient(1),
                                   gradient(2),
                                   diffuseCoefficients,
                                   specularCoefficients);
    VTKKWRCHelper_InterpolateShadingMfa(diffuseCoefficients, specularCoefficients, tmp); //diffuseShadingTable[0]/specularShadingTable[0] unsigned short pointer, starting address of the shading tables, range:[0, 65535]
    mfa_after = std::chrono::system_clock::now();
    auto shade_mfa_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(mfa_after - mfa_before);
#if defined(RECORD_TIMING) // record timing
    if (threadID == 0) {
        shade_mfa_time0 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter0++;
    }
    if (threadID == 1) {
        shade_mfa_time1 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter1++;
    }
    if (threadID == 2) {
        shade_mfa_time2 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter2++;
    }
    if (threadID == 3) {
        shade_mfa_time3 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter3++;
    }
    if (threadID == 4) {
        shade_mfa_time4 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter4++;
    }
    if (threadID == 5) {
        shade_mfa_time5 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter5++;
    }
    if (threadID == 6) {
        shade_mfa_time6 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter6++;
    }
    if (threadID == 7) {
        shade_mfa_time7 += shade_mfa_nanoseconds.count();
        _shade_mfa_counter7++;
    }
#endif
#endif
    /* Save Debug Info */
#if 0 
    if (threadID == 0) {
      std::ofstream filedebug0;
      filedebug0.open("debug0.txt", std::ofstream::out | std::ofstream::app);
      filedebug0 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug0.close();
    }
    if (threadID == 1) {
      std::ofstream filedebug1;
      filedebug1.open("debug1.txt", std::ofstream::out | std::ofstream::app);
      filedebug1 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug1.close();
    }
    if (threadID == 2) {
      std::ofstream filedebug2;
      filedebug2.open("debug2.txt", std::ofstream::out | std::ofstream::app);
      filedebug2 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug2.close();
    }
    if (threadID == 3) {
      std::ofstream filedebug3;
      filedebug3.open("debug3.txt", std::ofstream::out | std::ofstream::app);
      filedebug3 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug3.close();
    }
    if (threadID == 4) {
      std::ofstream filedebug4;
      filedebug4.open("debug4.txt", std::ofstream::out | std::ofstream::app);
      filedebug4 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug4.close();
    }
    if (threadID == 5) {
      std::ofstream filedebug5;
      filedebug5.open("debug5.txt", std::ofstream::out | std::ofstream::app);
      filedebug5 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug5.close();
    }
    if (threadID == 6) {
      std::ofstream filedebug6;
      filedebug6.open("debug6.txt", std::ofstream::out | std::ofstream::app);
      filedebug6 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug6.close();
    }
    if (threadID == 7) {
      std::ofstream filedebug7;
      filedebug7.open("debug7.txt", std::ofstream::out | std::ofstream::app);
      filedebug7 << j << "," << i << ","
                 << pos[0] << "," << pos[1] << "," << pos[2] << ","
                 // << pos[0] << "," << pos[1] << "," << pos[2] << "," << A << "," << B << "," << C << "," << D << "," << E << "," << F << "," << G << "," << H << "," << val << ","
                 << _tmpDColor[0] << "," << _tmpDColor[1] << "," << _tmpDColor[2] << ","
                 << _tmpSColor[0] << "," << _tmpSColor[1] << "," << _tmpSColor[2] << ","
                 << diffuseCoefficients[0] << "," << diffuseCoefficients[1] << "," << diffuseCoefficients[2] << ","
                 << specularCoefficients[0] << "," << specularCoefficients[1] << "," << specularCoefficients[2] << ","
                 << out_dx(3) << ", " << out_dy(3) << ", " << out_dz(3) << "\n";
      filedebug7.close();
    }
    // ------- test end ------- //
#endif

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();

#if defined(RECORD_TIMING) // record timing
  /* Save Time-stamp */ 
  if (threadID == 0) {
      std::ofstream value_mfa_time0_file;
      value_mfa_time0_file.open("shade/value_mfa_time0.txt");
      value_mfa_time0_file << value_mfa_time0;
      value_mfa_time0_file.close();

      std::ofstream shade_mfa_time0_file;
      shade_mfa_time0_file.open("shade/shade_mfa_time0.txt");
      shade_mfa_time0_file << shade_mfa_time0;
      shade_mfa_time0_file.close();

      std::ofstream mfa_count0_file;
      mfa_count0_file.open("shade/mfa_count0.txt");
      mfa_count0_file << shade_mfa_counter0;
      mfa_count0_file.close();

      std::ofstream _mfa_count0_file;
      _mfa_count0_file.open("shade/_mfa_count0.txt");
      _mfa_count0_file << _shade_mfa_counter0;
      _mfa_count0_file.close();
  }
  if (threadID == 1) {
      std::ofstream value_mfa_time1_file;
      value_mfa_time1_file.open("shade/value_mfa_time1.txt");
      value_mfa_time1_file << value_mfa_time1;
      value_mfa_time1_file.close();

      std::ofstream shade_mfa_time1_file;
      shade_mfa_time1_file.open("shade/shade_mfa_time1.txt");
      shade_mfa_time1_file << shade_mfa_time1;
      shade_mfa_time1_file.close();

      std::ofstream mfa_count1_file;
      mfa_count1_file.open("shade/mfa_count1.txt");
      mfa_count1_file << shade_mfa_counter1;
      mfa_count1_file.close();

      std::ofstream _mfa_count1_file;
      _mfa_count1_file.open("shade/_mfa_count1.txt");
      _mfa_count1_file << _shade_mfa_counter1;
      _mfa_count1_file.close();
  }
  if (threadID == 2) {
      std::ofstream value_mfa_time2_file;
      value_mfa_time2_file.open("shade/value_mfa_time2.txt");
      value_mfa_time2_file << value_mfa_time2;
      value_mfa_time2_file.close();

      std::ofstream shade_mfa_time2_file;
      shade_mfa_time2_file.open("shade/shade_mfa_time2.txt");
      shade_mfa_time2_file << shade_mfa_time2;
      shade_mfa_time2_file.close();

      std::ofstream mfa_count2_file;
      mfa_count2_file.open("shade/mfa_count2.txt");
      mfa_count2_file << shade_mfa_counter2;
      mfa_count2_file.close();

      std::ofstream _mfa_count2_file;
      _mfa_count2_file.open("shade/_mfa_count2.txt");
      _mfa_count2_file << _shade_mfa_counter2;
      _mfa_count2_file.close();
  }
  if (threadID == 3) {
      std::ofstream value_mfa_time3_file;
      value_mfa_time3_file.open("shade/value_mfa_time3.txt");
      value_mfa_time3_file << value_mfa_time3;
      value_mfa_time3_file.close();

      std::ofstream shade_mfa_time3_file;
      shade_mfa_time3_file.open("shade/shade_mfa_time3.txt");
      shade_mfa_time3_file << shade_mfa_time3;
      shade_mfa_time3_file.close();

      std::ofstream mfa_count3_file;
      mfa_count3_file.open("shade/mfa_count3.txt");
      mfa_count3_file << shade_mfa_counter3;
      mfa_count3_file.close();

      std::ofstream _mfa_count3_file;
      _mfa_count3_file.open("shade/_mfa_count3.txt");
      _mfa_count3_file << _shade_mfa_counter3;
      _mfa_count3_file.close();
  }
  if (threadID == 4) {
      std::ofstream value_mfa_time4_file;
      value_mfa_time4_file.open("shade/value_mfa_time4.txt");
      value_mfa_time4_file << value_mfa_time4;
      value_mfa_time4_file.close();

      std::ofstream shade_mfa_time4_file;
      shade_mfa_time4_file.open("shade/shade_mfa_time4.txt");
      shade_mfa_time4_file << shade_mfa_time4;
      shade_mfa_time4_file.close();

      std::ofstream mfa_count4_file;
      mfa_count4_file.open("shade/mfa_count4.txt");
      mfa_count4_file << shade_mfa_counter4;
      mfa_count4_file.close();

      std::ofstream _mfa_count4_file;
      _mfa_count4_file.open("shade/_mfa_count4.txt");
      _mfa_count4_file << _shade_mfa_counter4;
      _mfa_count4_file.close();
  }
  if (threadID == 5) {
      std::ofstream value_mfa_time5_file;
      value_mfa_time5_file.open("shade/value_mfa_time5.txt");
      value_mfa_time5_file << value_mfa_time5;
      value_mfa_time5_file.close();

      std::ofstream shade_mfa_time5_file;
      shade_mfa_time5_file.open("shade/shade_mfa_time5.txt");
      shade_mfa_time5_file << shade_mfa_time5;
      shade_mfa_time5_file.close();

      std::ofstream mfa_count5_file;
      mfa_count5_file.open("shade/mfa_count5.txt");
      mfa_count5_file << shade_mfa_counter5;
      mfa_count5_file.close();

      std::ofstream _mfa_count5_file;
      _mfa_count5_file.open("shade/_mfa_count5.txt");
      _mfa_count5_file << _shade_mfa_counter5;
      _mfa_count5_file.close();
  }
  if (threadID == 6) {
      std::ofstream value_mfa_time6_file;
      value_mfa_time6_file.open("shade/value_mfa_time6.txt");
      value_mfa_time6_file << value_mfa_time6;
      value_mfa_time6_file.close();

      std::ofstream shade_mfa_time6_file;
      shade_mfa_time6_file.open("shade/shade_mfa_time6.txt");
      shade_mfa_time6_file << shade_mfa_time6;
      shade_mfa_time6_file.close();

      std::ofstream mfa_count6_file;
      mfa_count6_file.open("shade/mfa_count6.txt");
      mfa_count6_file << shade_mfa_counter6;
      mfa_count6_file.close();

      std::ofstream _mfa_count6_file;
      _mfa_count6_file.open("shade/_mfa_count6.txt");
      _mfa_count6_file << _shade_mfa_counter6;
      _mfa_count6_file.close();
  }
  if (threadID == 7) {
      std::ofstream value_mfa_time7_file;
      value_mfa_time7_file.open("shade/value_mfa_time7.txt");
      value_mfa_time7_file << value_mfa_time7;
      value_mfa_time7_file.close();

      std::ofstream shade_mfa_time7_file;
      shade_mfa_time7_file.open("shade/shade_mfa_time7.txt");
      shade_mfa_time7_file << shade_mfa_time7;
      shade_mfa_time7_file.close();

      std::ofstream mfa_count7_file;
      mfa_count7_file.open("shade/mfa_count7.txt");
      mfa_count7_file << shade_mfa_counter7;
      mfa_count7_file.close();

      std::ofstream _mfa_count7_file;
      _mfa_count7_file.open("shade/_mfa_count7.txt");
      _mfa_count7_file << _shade_mfa_counter7;
      _mfa_count7_file.close();
  }
#endif
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

    // VTKKWRCHelper_ComputeWeights(pos);
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
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol))
  {
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
    // One component
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0)
      {
        if (!mapper->GetUseMfa()) {
          switch (scalarType)
          {
            vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleTrilin(static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
          }
        } else {
          switch (scalarType)
          {
            vtkTemplateMacro(vtkFixedPointCompositeShadeHelperGenerateImageOneSimpleTrilinMfa(static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
          }
        }
      }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
      {
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
