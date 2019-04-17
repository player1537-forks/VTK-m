//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/worklet/NDimsHistMarginalization.h>
#include <vtkm/worklet/NDimsHistogram.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/testing/Testing.h>

namespace
{
// Make testing dataset with three fields(variables), each one has 1000 values
vtkm::cont::DataSet MakeTestDataSet()
{
  vtkm::cont::DataSet dataSet;

  const int xVerts = 20;
  const int yVerts = 50;
  const int nVerts = xVerts * yVerts;

  vtkm::Float32 fieldA[nVerts] = {
    8,  10, 9,  8,  14, 11, 12, 9,  19, 7,  8,  11, 7,  10, 11, 11, 11, 6,  8,  8,  7,  15, 9,  7,
    8,  10, 9,  10, 10, 12, 7,  6,  14, 10, 14, 10, 7,  11, 13, 9,  13, 11, 10, 10, 12, 12, 7,  12,
    10, 11, 12, 8,  13, 9,  5,  12, 11, 9,  5,  9,  12, 9,  6,  10, 11, 9,  9,  11, 9,  7,  7,  18,
    16, 13, 12, 8,  10, 11, 9,  8,  17, 3,  15, 15, 9,  10, 10, 8,  10, 9,  7,  9,  8,  10, 13, 9,
    7,  11, 7,  10, 13, 10, 11, 9,  10, 7,  10, 6,  12, 6,  9,  7,  6,  12, 12, 9,  12, 12, 11, 6,
    1,  12, 8,  13, 14, 8,  8,  10, 7,  7,  6,  7,  5,  11, 6,  11, 13, 8,  13, 5,  9,  12, 7,  11,
    10, 15, 11, 9,  7,  12, 15, 7,  8,  7,  12, 8,  21, 16, 13, 11, 10, 14, 12, 11, 12, 14, 7,  11,
    7,  12, 16, 8,  10, 8,  9,  7,  8,  7,  13, 13, 11, 15, 7,  7,  6,  11, 7,  12, 12, 13, 14, 11,
    13, 11, 11, 9,  15, 8,  6,  11, 12, 10, 11, 7,  6,  14, 11, 10, 12, 5,  8,  9,  11, 15, 11, 10,
    17, 14, 9,  10, 10, 12, 11, 13, 13, 12, 11, 7,  8,  10, 7,  11, 10, 5,  8,  10, 13, 13, 12, 6,
    10, 7,  13, 8,  11, 7,  10, 7,  8,  7,  14, 16, 9,  11, 8,  11, 9,  15, 11, 10, 10, 12, 7,  7,
    11, 7,  5,  17, 9,  11, 11, 11, 10, 17, 10, 15, 7,  11, 12, 16, 9,  8,  11, 14, 9,  22, 8,  8,
    8,  13, 12, 12, 1,  14, 15, 6,  15, 8,  11, 16, 14, 8,  6,  9,  8,  9,  9,  10, 8,  6,  13, 8,
    6,  12, 11, 12, 13, 8,  6,  6,  5,  6,  10, 9,  11, 12, 14, 12, 10, 11, 10, 10, 8,  13, 8,  11,
    7,  13, 13, 12, 12, 13, 15, 4,  9,  16, 7,  9,  8,  10, 6,  9,  11, 12, 6,  7,  14, 6,  4,  15,
    5,  18, 9,  9,  11, 12, 9,  5,  6,  7,  15, 6,  11, 14, 8,  12, 6,  9,  5,  9,  14, 9,  12, 6,
    9,  14, 11, 12, 12, 13, 15, 9,  8,  7,  13, 12, 7,  13, 6,  9,  10, 10, 10, 9,  11, 5,  9,  13,
    16, 9,  10, 8,  9,  6,  13, 12, 8,  12, 9,  12, 17, 8,  11, 10, 8,  7,  11, 7,  13, 13, 10, 14,
    11, 9,  6,  6,  14, 16, 5,  9,  13, 11, 12, 7,  4,  6,  9,  11, 11, 10, 12, 9,  7,  13, 8,  8,
    12, 5,  10, 7,  11, 11, 10, 10, 14, 6,  8,  8,  3,  12, 16, 11, 11, 7,  6,  12, 11, 5,  9,  12,
    9,  13, 7,  8,  9,  9,  12, 7,  9,  8,  12, 11, 6,  10, 6,  7,  6,  11, 10, 8,  9,  8,  4,  19,
    12, 6,  10, 9,  6,  12, 9,  14, 7,  8,  11, 7,  7,  12, 13, 9,  13, 12, 8,  6,  10, 17, 19, 10,
    10, 13, 5,  11, 8,  10, 8,  16, 12, 6,  6,  7,  10, 9,  12, 8,  5,  10, 7,  18, 9,  12, 10, 4,
    9,  9,  15, 15, 6,  7,  7,  11, 12, 4,  8,  18, 5,  12, 12, 11, 10, 14, 9,  9,  10, 8,  10, 8,
    10, 9,  9,  4,  10, 12, 5,  13, 6,  9,  7,  5,  12, 8,  11, 10, 9,  17, 9,  9,  8,  11, 18, 11,
    10, 9,  4,  13, 10, 15, 5,  10, 9,  7,  7,  8,  10, 6,  6,  19, 10, 16, 7,  7,  9,  10, 10, 13,
    10, 10, 14, 13, 12, 8,  7,  13, 12, 11, 13, 12, 9,  8,  6,  8,  10, 3,  8,  8,  12, 12, 13, 13,
    10, 5,  10, 7,  13, 7,  9,  5,  13, 7,  10, 8,  13, 11, 17, 9,  6,  14, 10, 10, 13, 9,  15, 8,
    15, 9,  12, 11, 12, 8,  3,  9,  8,  10, 12, 8,  14, 13, 12, 11, 12, 9,  18, 10, 13, 7,  4,  4,
    11, 8,  3,  7,  9,  10, 12, 7,  11, 21, 9,  7,  8,  9,  10, 10, 11, 9,  15, 13, 21, 12, 8,  11,
    9,  10, 11, 9,  17, 8,  9,  8,  14, 6,  13, 9,  8,  11, 12, 12, 12, 11, 6,  13, 7,  9,  11, 15,
    17, 17, 11, 10, 7,  8,  11, 8,  6,  9,  13, 7,  9,  6,  5,  10, 7,  16, 16, 9,  7,  6,  14, 8,
    13, 16, 7,  7,  10, 11, 6,  10, 9,  9,  8,  14, 11, 9,  11, 9,  10, 11, 9,  8,  14, 11, 7,  12,
    11, 8,  9,  9,  10, 11, 11, 10, 9,  6,  6,  11, 16, 10, 7,  6,  6,  13, 18, 8,  12, 11, 14, 13,
    8,  8,  10, 17, 17, 6,  6,  10, 18, 5,  8,  11, 6,  6,  14, 10, 9,  6,  11, 6,  13, 12, 10, 6,
    9,  9,  9,  13, 7,  17, 10, 14, 10, 9,  10, 10, 11, 10, 11, 15, 13, 6,  12, 19, 10, 12, 12, 15,
    13, 10, 10, 13, 11, 13, 13, 17, 6,  5,  6,  7,  6,  9,  13, 11, 8,  12, 9,  6,  10, 16, 11, 12,
    5,  12, 14, 13, 13, 16, 11, 6,  12, 12, 15, 8,  7,  11, 8,  5,  10, 8,  9,  11, 9,  12, 10, 5,
    12, 11, 9,  6,  14, 12, 10, 11, 9,  6,  7,  12, 8,  12, 8,  15, 9,  8,  7,  9,  3,  6,  14, 7,
    8,  11, 9,  10, 12, 9,  10, 9,  8,  6,  12, 11, 6,  8,  9,  8,  15, 11, 7,  18, 12, 11, 10, 13,
    11, 11, 10, 7,  9,  8,  8,  11, 11, 13, 6,  12, 13, 16, 11, 11, 5,  12, 14, 15, 9,  14, 15, 6,
    8,  7,  6,  8,  9,  19, 7,  12, 11, 8,  14, 12, 10, 9,  3,  7
  };

  vtkm::Float32 fieldB[nVerts] = {
    24, 19, 28, 19, 25, 28, 25, 22, 27, 26, 35, 26, 30, 28, 24, 23, 21, 31, 20, 11, 21, 22, 14, 25,
    20, 24, 24, 21, 24, 29, 26, 21, 32, 29, 23, 28, 31, 25, 23, 30, 18, 24, 22, 25, 33, 24, 22, 23,
    21, 17, 20, 28, 30, 18, 20, 32, 25, 24, 32, 15, 27, 24, 27, 19, 30, 27, 17, 24, 29, 23, 22, 19,
    24, 19, 28, 24, 25, 24, 25, 30, 24, 31, 30, 27, 25, 25, 25, 15, 29, 23, 29, 29, 21, 25, 35, 24,
    28, 10, 31, 23, 22, 22, 22, 33, 29, 27, 18, 27, 27, 24, 20, 20, 21, 29, 23, 31, 23, 23, 22, 23,
    30, 27, 28, 31, 16, 29, 25, 19, 33, 28, 25, 24, 15, 27, 37, 29, 15, 19, 14, 19, 24, 23, 30, 29,
    35, 22, 19, 26, 26, 14, 24, 30, 32, 23, 30, 29, 26, 27, 25, 23, 17, 26, 32, 29, 20, 17, 21, 23,
    22, 20, 36, 12, 26, 23, 15, 29, 24, 22, 26, 33, 24, 23, 20, 26, 22, 17, 26, 26, 34, 22, 26, 17,
    23, 18, 29, 27, 21, 29, 28, 29, 24, 25, 28, 19, 18, 21, 23, 23, 27, 25, 24, 25, 24, 25, 21, 25,
    21, 27, 23, 20, 29, 15, 28, 30, 24, 27, 17, 23, 16, 21, 25, 17, 27, 28, 21, 13, 19, 27, 16, 30,
    31, 25, 30, 17, 17, 25, 26, 22, 21, 17, 24, 17, 25, 22, 27, 14, 27, 24, 27, 25, 26, 31, 21, 23,
    30, 30, 22, 19, 23, 22, 23, 25, 24, 25, 24, 28, 26, 30, 18, 25, 30, 37, 27, 34, 28, 34, 25, 10,
    25, 22, 35, 30, 24, 32, 24, 34, 19, 29, 26, 16, 27, 17, 26, 23, 27, 25, 26, 21, 31, 21, 28, 15,
    32, 24, 23, 23, 18, 15, 22, 25, 16, 25, 31, 26, 25, 28, 24, 26, 23, 25, 33, 20, 27, 28, 24, 29,
    32, 20, 24, 20, 19, 32, 24, 6,  24, 21, 26, 18, 15, 30, 19, 26, 22, 30, 35, 23, 22, 30, 20, 22,
    18, 30, 28, 25, 16, 25, 27, 30, 18, 24, 30, 28, 20, 19, 20, 28, 21, 24, 15, 33, 20, 18, 20, 36,
    30, 26, 25, 18, 28, 27, 31, 31, 15, 26, 16, 22, 27, 14, 17, 27, 27, 22, 32, 30, 22, 34, 22, 25,
    20, 22, 26, 29, 28, 33, 18, 23, 20, 20, 27, 24, 28, 21, 25, 27, 25, 19, 19, 25, 19, 32, 29, 27,
    23, 21, 28, 33, 23, 23, 28, 26, 31, 19, 21, 29, 21, 27, 23, 32, 24, 26, 21, 28, 28, 24, 17, 31,
    27, 21, 19, 32, 28, 23, 30, 23, 29, 15, 26, 26, 15, 20, 25, 26, 27, 31, 21, 23, 23, 33, 28, 19,
    23, 22, 22, 25, 27, 17, 23, 17, 25, 28, 26, 30, 32, 31, 19, 25, 25, 19, 23, 29, 27, 23, 34, 22,
    13, 21, 32, 10, 20, 33, 21, 17, 29, 31, 14, 24, 23, 19, 19, 22, 17, 26, 37, 26, 22, 26, 38, 29,
    29, 27, 30, 20, 31, 14, 32, 32, 24, 23, 23, 18, 21, 31, 24, 20, 28, 15, 21, 25, 25, 20, 30, 25,
    22, 21, 21, 25, 24, 25, 18, 23, 28, 30, 20, 27, 27, 19, 10, 32, 24, 20, 29, 26, 25, 20, 25, 29,
    28, 24, 32, 26, 22, 19, 23, 27, 27, 29, 20, 25, 21, 30, 28, 31, 24, 19, 23, 19, 19, 18, 30, 18,
    16, 24, 20, 20, 30, 25, 29, 25, 31, 21, 28, 31, 24, 26, 27, 21, 24, 23, 26, 18, 32, 26, 28, 26,
    24, 26, 29, 30, 22, 20, 24, 28, 25, 29, 20, 21, 22, 15, 30, 27, 33, 26, 22, 32, 30, 31, 20, 19,
    24, 26, 27, 31, 17, 17, 33, 27, 16, 27, 27, 22, 27, 19, 24, 21, 17, 24, 28, 23, 26, 24, 19, 26,
    20, 24, 22, 19, 22, 21, 21, 28, 29, 39, 19, 16, 25, 29, 31, 22, 22, 29, 26, 22, 22, 22, 26, 23,
    23, 23, 30, 25, 25, 25, 27, 29, 18, 33, 21, 12, 22, 29, 12, 20, 35, 22, 34, 28, 18, 29, 21, 20,
    24, 33, 24, 26, 23, 34, 31, 25, 31, 22, 35, 21, 20, 29, 27, 22, 30, 22, 27, 23, 22, 32, 16, 19,
    27, 22, 24, 27, 21, 33, 25, 25, 19, 28, 20, 27, 21, 25, 28, 20, 27, 22, 21, 20, 26, 30, 33, 23,
    20, 24, 17, 23, 28, 35, 14, 23, 22, 28, 28, 26, 25, 18, 20, 28, 28, 22, 13, 24, 22, 20, 30, 26,
    26, 18, 22, 20, 23, 24, 20, 27, 34, 28, 18, 24, 34, 33, 25, 33, 37, 21, 20, 31, 19, 23, 29, 22,
    21, 24, 19, 27, 19, 32, 25, 23, 33, 26, 33, 27, 29, 30, 19, 22, 30, 19, 18, 24, 25, 17, 31, 19,
    31, 26, 22, 23, 28, 28, 25, 24, 19, 19, 27, 28, 23, 21, 29, 26, 31, 22, 22, 25, 16, 29, 21, 22,
    23, 25, 22, 21, 22, 19, 27, 26, 28, 30, 22, 21, 24, 22, 23, 26, 28, 22, 18, 25, 23, 27, 31, 19,
    15, 29, 20, 19, 27, 25, 21, 29, 22, 24, 25, 17, 36, 29, 22, 22, 24, 28, 27, 22, 26, 31, 29, 31,
    18, 25, 23, 16, 37, 27, 21, 31, 25, 24, 20, 23, 28, 33, 24, 21, 26, 20, 18, 31, 20, 24, 23, 19,
    27, 17, 23, 23, 20, 26, 28, 23, 26, 31, 25, 31, 19, 32, 26, 18, 19, 29, 20, 21, 15, 25, 27, 29,
    22, 22, 22, 26, 23, 22, 23, 29, 28, 20, 21, 22, 20, 22, 27, 25, 23, 32, 23, 20, 31, 20, 27, 26,
    34, 20, 22, 36, 21, 29, 25, 20, 21, 22, 29, 29, 25, 22, 24, 22
  };

  vtkm::Float32 fieldC[nVerts] = {
    3,  1,  4,  6,  5,  4,  8,  7,  2,  9,  2,  0,  0,  4,  3,  2,  5,  2,  3,  6,  3,  8,  3,  4,
    3,  3,  2,  7,  2,  10, 9,  6,  1,  1,  4,  7,  3,  3,  1,  4,  4,  3,  9,  4,  4,  7,  3,  2,
    4,  7,  3,  3,  2,  10, 1,  6,  2,  2,  3,  8,  3,  3,  6,  9,  4,  1,  4,  3,  16, 7,  0,  1,
    8,  7,  13, 3,  5,  0,  3,  8,  10, 3,  5,  5,  1,  5,  2,  1,  3,  2,  5,  3,  4,  3,  3,  3,
    3,  1,  13, 2,  3,  1,  2,  7,  3,  4,  1,  2,  5,  4,  4,  4,  2,  6,  3,  2,  7,  8,  1,  3,
    4,  1,  2,  0,  1,  6,  1,  8,  8,  1,  1,  4,  2,  1,  4,  3,  5,  4,  6,  4,  2,  3,  8,  8,
    3,  3,  3,  4,  5,  8,  8,  16, 7,  12, 4,  3,  14, 8,  3,  12, 5,  0,  5,  3,  5,  2,  9,  2,
    9,  4,  1,  0,  0,  4,  4,  6,  3,  4,  11, 2,  4,  7,  4,  2,  1,  9,  4,  3,  2,  5,  1,  5,
    3,  8,  2,  8,  1,  8,  0,  4,  1,  3,  2,  1,  2,  3,  2,  1,  8,  5,  4,  1,  9,  9,  1,  3,
    5,  0,  1,  6,  10, 8,  3,  12, 3,  4,  4,  7,  1,  3,  6,  4,  4,  6,  1,  4,  7,  5,  6,  11,
    6,  5,  2,  7,  2,  5,  3,  5,  6,  3,  6,  2,  1,  10, 8,  3,  7,  0,  2,  6,  9,  3,  11, 3,
    2,  5,  1,  4,  6,  10, 9,  1,  4,  3,  7,  12, 3,  10, 0,  2,  11, 2,  1,  0,  4,  1,  2,  16,
    5,  17, 7,  8,  2,  10, 10, 3,  1,  3,  2,  2,  4,  8,  4,  3,  2,  4,  4,  6,  8,  6,  2,  3,
    2,  4,  2,  4,  7,  10, 5,  3,  5,  2,  4,  6,  9,  3,  1,  1,  1,  1,  4,  2,  2,  7,  4,  9,
    2,  3,  5,  6,  2,  5,  1,  6,  5,  7,  8,  3,  7,  2,  2,  8,  6,  2,  10, 2,  1,  4,  5,  1,
    1,  1,  5,  6,  1,  1,  4,  5,  4,  2,  4,  3,  2,  7,  19, 4,  7,  2,  7,  5,  2,  5,  3,  8,
    4,  6,  7,  2,  0,  0,  2,  12, 6,  2,  2,  3,  5,  9,  4,  9,  2,  2,  7,  8,  3,  3,  10, 6,
    3,  2,  1,  6,  2,  4,  6,  3,  5,  8,  2,  3,  6,  14, 0,  3,  6,  5,  2,  7,  0,  3,  8,  5,
    3,  2,  2,  5,  1,  3,  12, 11, 16, 2,  1,  3,  7,  3,  1,  6,  4,  3,  12, 5,  1,  3,  1,  4,
    9,  1,  3,  3,  4,  4,  6,  7,  7,  5,  2,  4,  2,  3,  2,  2,  6,  4,  2,  2,  3,  5,  1,  4,
    9,  1,  0,  7,  6,  4,  3,  3,  7,  3,  3,  6,  2,  7,  9,  3,  1,  16, 5,  4,  3,  6,  3,  2,
    5,  2,  2,  4,  3,  1,  3,  3,  6,  3,  5,  9,  1,  10, 1,  7,  2,  2,  6,  7,  3,  5,  3,  7,
    2,  2,  2,  2,  6,  4,  3,  2,  5,  5,  3,  15, 4,  2,  7,  7,  4,  3,  3,  5,  1,  2,  9,  0,
    5,  7,  12, 2,  4,  8,  5,  7,  8,  3,  2,  2,  18, 1,  7,  2,  2,  1,  3,  3,  3,  7,  1,  9,
    8,  4,  3,  7,  6,  4,  5,  2,  0,  5,  1,  5,  10, 4,  2,  8,  2,  2,  0,  5,  6,  4,  5,  0,
    1,  5,  11, 3,  3,  4,  4,  2,  3,  5,  1,  6,  5,  7,  2,  2,  5,  7,  4,  8,  4,  1,  1,  7,
    2,  3,  9,  6,  13, 1,  5,  4,  6,  2,  4,  11, 2,  5,  5,  1,  4,  1,  4,  7,  1,  5,  8,  3,
    1,  10, 9,  13, 1,  7,  2,  9,  4,  3,  3,  10, 12, 2,  0,  4,  6,  5,  5,  1,  4,  7,  2,  12,
    7,  6,  5,  0,  6,  4,  4,  12, 1,  3,  10, 1,  9,  2,  2,  2,  1,  5,  5,  6,  9,  6,  4,  1,
    11, 6,  9,  3,  2,  7,  1,  7,  4,  3,  0,  3,  1,  12, 17, 2,  1,  6,  4,  4,  2,  1,  5,  5,
    3,  2,  2,  4,  6,  5,  4,  6,  11, 3,  12, 6,  3,  6,  3,  0,  6,  3,  7,  4,  8,  5,  14, 5,
    1,  9,  4,  6,  5,  3,  9,  3,  1,  1,  0,  3,  7,  3,  5,  1,  6,  2,  2,  6,  2,  12, 1,  0,
    6,  3,  3,  5,  4,  7,  2,  2,  15, 7,  3,  10, 4,  2,  6,  3,  4,  8,  3,  1,  5,  5,  5,  4,
    3,  7,  3,  4,  5,  5,  2,  4,  2,  5,  1,  12, 5,  6,  3,  2,  8,  5,  2,  3,  11, 11, 6,  5,
    0,  3,  3,  9,  4,  2,  11, 1,  5,  3,  5,  6,  3,  6,  4,  2,  4,  10, 11, 3,  3,  4,  1,  1,
    1,  3,  5,  5,  1,  1,  4,  1,  5,  1,  6,  8,  6,  4,  6,  7,  6,  3,  5,  3,  6,  6,  6,  4,
    0,  6,  3,  1,  2,  4,  2,  6,  1,  1,  1,  2,  2,  4,  7,  2,  6,  2,  5,  7,  6,  4,  6,  3,
    1,  4,  5,  1,  4,  6,  2,  3,  0,  6,  11, 2,  9,  2,  6,  4,  5,  6,  2,  19, 2,  10, 4,  2,
    3,  3,  11, 7,  3,  3,  1,  5,  3,  6,  4,  3,  0,  6,  6,  6,  4,  2,  5,  2,  2,  2,  6,  10,
    4,  9,  3,  7,  7,  0,  6,  8,  5,  2,  3,  2,  3,  3,  3,  1,  6,  1,  8,  2,  5,  3,  6,  11,
    5,  7,  2,  6,  7,  3,  4,  1,  0,  5,  8,  3,  2,  9,  3,  1,  2,  3,  3,  9,  5,  6,  5,  1,
    4,  5,  6,  7,  6,  1,  5,  1,  6,  6,  2,  6,  7,  2,  4,  6
  };

  vtkm::cont::ArrayHandleUniformPointCoordinates coordinates(vtkm::Id3(xVerts, yVerts, 1));
  dataSet.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coordinates", coordinates));

  // Set point scalars
  dataSet.AddField(vtkm::cont::make_Field(
    "fieldA", vtkm::cont::Field::Association::POINTS, fieldA, nVerts, vtkm::CopyFlag::On));
  dataSet.AddField(vtkm::cont::make_Field(
    "fieldB", vtkm::cont::Field::Association::POINTS, fieldB, nVerts, vtkm::CopyFlag::On));
  dataSet.AddField(vtkm::cont::make_Field(
    "fieldC", vtkm::cont::Field::Association::POINTS, fieldC, nVerts, vtkm::CopyFlag::On));

  return dataSet;
}

// The Condition function for non-marginal variable
// var is index of variable (data array)
// binId is bin index
// return true indicates considering this bin into final
// marginal histogram
// In this example, we have three variable var0, var1, var2
// the condition is P(Var0, Var2 | 1<Var1<4)
// because var0 and var2 are the marginal variable, we do not care the case var==0 or var==2
// it supposes that we do not have input as var ==0 or 2
struct VariableCondition
{
  VTKM_EXEC
  bool operator()(vtkm::Id var, vtkm::Id binId) const
  {
    if (var == 1)
    {
      if (1 < binId && binId < 4)
        return true;
    }

    return false;
  }
};

//
// Create a dataset with known point data and cell data (statistical distributions)
// Extract arrays of point and cell fields
// Create output structure to hold histogram bins
// Run FieldHistogram filter
//
void TestNDimsHistMarginalization()
{
  // Setup and calculate a ND-histogram without marginalization first
  // Create the output bin array
  vtkm::cont::ArrayHandle<vtkm::Range> range;
  vtkm::cont::ArrayHandle<vtkm::Float64> delta;

  // Data attached is the poisson distribution
  vtkm::cont::DataSet ds = MakeTestDataSet();

  // Compute Nd histogram first
  vtkm::worklet::NDimsHistogram ndHistogram;

  // Set the number of data points
  ndHistogram.SetNumOfDataPoints(ds.GetField(0).GetNumberOfValues());

  // Add field one by one
  vtkm::Range rangeFieldA;
  vtkm::Float64 deltaFieldA;
  ndHistogram.AddField(ds.GetField("fieldA").GetData(), 10, rangeFieldA, deltaFieldA);

  vtkm::Range rangeFieldB;
  vtkm::Float64 deltaFieldB;
  ndHistogram.AddField(ds.GetField("fieldB").GetData(), 10, rangeFieldB, deltaFieldB);

  vtkm::Range rangeFieldC;
  vtkm::Float64 deltaFieldC;
  ndHistogram.AddField(ds.GetField("fieldC").GetData(), 10, rangeFieldC, deltaFieldC);

  // the return binIds and freqs is sparse distribution representation
  // (we do not keep the 0 frequency entities)
  // e.g. we have three variable(data arrays) in this example
  // binIds[0, 1, 2][j] is a combination of bin ID of three variable,
  // freqs[j] is the freqncy of this bin IDs combination
  std::vector<vtkm::cont::ArrayHandle<vtkm::Id>> binIds;
  vtkm::cont::ArrayHandle<vtkm::Id> freqs;
  ndHistogram.Run(binIds, freqs);

  // setup for histogram marginalization
  // use a bool array to indicate the marginal variable (true -> marginal variable)
  // the length of this array has to be equal to number of input variables
  bool marginalVariableAry[3] = { true, false, true };
  vtkm::cont::ArrayHandle<bool> marginalVariable =
    vtkm::cont::make_ArrayHandle(marginalVariableAry, 3);

  std::vector<vtkm::Id> numberOfBinsVec(3, 10);
  vtkm::cont::ArrayHandle<vtkm::Id> numberOfBins = vtkm::cont::make_ArrayHandle(numberOfBinsVec);

  // calculate marginal histogram by the setup (return sparse represetnation)
  std::vector<vtkm::cont::ArrayHandle<vtkm::Id>> marginalBinIds;
  vtkm::cont::ArrayHandle<vtkm::Id> marginalFreqs;
  vtkm::worklet::NDimsHistMarginalization ndHistMarginalization;
  ndHistMarginalization.Run(binIds,
                            freqs,
                            numberOfBins,
                            marginalVariable,
                            VariableCondition(),
                            marginalBinIds,
                            marginalFreqs);

  // Ground truth ND histogram
  vtkm::Id gtNonSparseBins = 40;
  vtkm::Id gtIdx0[40] = { 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4,
                          4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 9 };
  vtkm::Id gtIdx1[40] = { 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 7, 0, 1, 2, 3, 4, 5, 0, 1,
                          2, 3, 4, 5, 7, 8, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 1, 2, 0, 1 };
  vtkm::Id gtFreq[40] = { 1,  2, 1, 2, 1, 4, 7, 6, 3, 2, 2, 1, 6, 6, 8, 6, 2, 2, 6, 9,
                          10, 2, 5, 1, 1, 1, 6, 7, 9, 6, 3, 3, 2, 3, 2, 2, 3, 2, 1, 1 };

  // Check result
  vtkm::Id nonSparseBins = marginalBinIds[0].GetPortalControl().GetNumberOfValues();
  VTKM_TEST_ASSERT(nonSparseBins == gtNonSparseBins,
                   "Incorrect ND-histogram marginalization results");

  for (int i = 0; i < nonSparseBins; i++)
  {
    vtkm::Id idx0 = marginalBinIds[0].GetPortalControl().Get(i);
    vtkm::Id idx1 = marginalBinIds[1].GetPortalControl().Get(i);
    vtkm::Id f = marginalFreqs.GetPortalControl().Get(i);
    VTKM_TEST_ASSERT(idx0 == gtIdx0[i] && idx1 == gtIdx1[i] && f == gtFreq[i],
                     "Incorrect ND-histogram marginalization results");
  }
} // TestNDimsHistMarginalization
}

int UnitTestNDimsHistMarginalization(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestNDimsHistMarginalization, argc, argv);
}
