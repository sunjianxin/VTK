/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSMPToolsInternal.h.in

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMPToolsInternal_h
#define vtkSMPToolsInternal_h

#include "vtkCommonCoreModule.h"       // For export macro
#include "vtkSMPToolsInternalCommon.h" // For common vtk smp class

#include <algorithm> // For std::sort()

#ifndef __VTK_WRAP__
namespace vtk
{
namespace detail
{
namespace smp
{

typedef void (*ExecuteFunctorPtrType)(void*, vtkIdType, vtkIdType, vtkIdType);

int VTKCOMMONCORE_EXPORT GetNumberOfThreads();
void VTKCOMMONCORE_EXPORT vtkSMPTools_Impl_For_STD(vtkIdType first, vtkIdType last, vtkIdType grain,
  ExecuteFunctorPtrType functorExecuter, void* functor);

template <typename FunctorInternal>
void ExecuteFunctor(void* functor, vtkIdType from, vtkIdType grain, vtkIdType last)
{
  const vtkIdType to = std::min(from + grain, last);

  FunctorInternal& fi = *reinterpret_cast<FunctorInternal*>(functor);
  fi.Execute(from, to);
}

template <typename FunctorInternal>
void vtkSMPTools_Impl_For(vtkIdType first, vtkIdType last, vtkIdType grain, FunctorInternal& fi)
{
  vtkIdType n = last - first;
  if (n <= 0)
  {
    return;
  }

  if (grain >= n)
  {
    fi.Execute(first, last);
  }
  else
  {
    vtkSMPTools_Impl_For_STD(first, last, grain, ExecuteFunctor<FunctorInternal>, &fi);
  }
}

//--------------------------------------------------------------------------------
template <typename InputIt, typename OutputIt, typename Functor>
void vtkSMPTools_Impl_Transform(
  InputIt inBegin, InputIt inEnd, OutputIt outBegin, Functor transform)
{
  auto size = std::distance(inBegin, inEnd);

  UnaryTransformCall<InputIt, OutputIt, Functor> exec(inBegin, outBegin, transform);
  vtkSMPTools_Impl_For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <typename InputIt1, typename InputIt2, typename OutputIt, typename Functor>
static void vtkSMPTools_Impl_Transform(
  InputIt1 inBegin1, InputIt1 inEnd, InputIt2 inBegin2, OutputIt outBegin, Functor transform)
{
  auto size = std::distance(inBegin1, inEnd);

  BinaryTransformCall<InputIt1, InputIt2, OutputIt, Functor> exec(
    inBegin1, inBegin2, outBegin, transform);
  vtkSMPTools_Impl_For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <typename Iterator, typename T>
void vtkSMPTools_Impl_Fill(Iterator begin, Iterator end, const T& value)
{
  auto size = std::distance(begin, end);

  FillFunctor<T> fill(value);
  UnaryTransformCall<Iterator, Iterator, FillFunctor<T>> exec(begin, begin, fill);
  vtkSMPTools_Impl_For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <typename RandomAccessIterator>
void vtkSMPTools_Impl_Sort(RandomAccessIterator begin, RandomAccessIterator end)
{
  std::sort(begin, end);
}

//--------------------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare>
void vtkSMPTools_Impl_Sort(RandomAccessIterator begin, RandomAccessIterator end, Compare comp)
{
  std::sort(begin, end, comp);
}

} // namespace smp
} // namespace detail
} // namespace vtk

#endif // __VTK_WRAP__

#endif
