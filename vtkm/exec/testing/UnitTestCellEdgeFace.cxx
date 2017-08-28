//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2016 Sandia Corporation.
//  Copyright 2016 UT-Battelle, LLC.
//  Copyright 2016 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#include <vtkm/exec/CellEdge.h>
#include <vtkm/exec/CellFace.h>

#include <vtkm/CellShape.h>
#include <vtkm/CellTraits.h>

#include <vtkm/exec/FunctorBase.h>

#include <vtkm/testing/Testing.h>

#include <set>
#include <vector>

namespace
{

using EdgeType = vtkm::Vec<vtkm::IdComponent, 2>;

void MakeEdgeCononical(EdgeType& edge)
{
  if (edge[1] < edge[0])
  {
    std::swap(edge[0], edge[1]);
  }
}

struct TestCellFacesFunctor
{
  template <typename CellShapeTag>
  void DoTest(vtkm::IdComponent numPoints,
              CellShapeTag shape,
              vtkm::CellTopologicalDimensionsTag<3>) const
  {
    // Stuff to fake running in the execution environment.
    char messageBuffer[256];
    messageBuffer[0] = '\0';
    vtkm::exec::internal::ErrorMessageBuffer errorMessage(messageBuffer, 256);
    vtkm::exec::FunctorBase workletProxy;
    workletProxy.SetErrorMessageBuffer(errorMessage);

    std::vector<vtkm::Id> pointIndexProxyBuffer(static_cast<std::size_t>(numPoints));
    for (std::size_t index = 0; index < pointIndexProxyBuffer.size(); ++index)
    {
      pointIndexProxyBuffer[index] = static_cast<vtkm::Id>(1000000 - index);
    }
    vtkm::VecCConst<vtkm::Id> pointIndexProxy(&pointIndexProxyBuffer.at(0), numPoints);

    vtkm::IdComponent numEdges = vtkm::exec::CellEdgeNumberOfEdges(numPoints, shape, workletProxy);
    VTKM_TEST_ASSERT(numEdges > 0, "No edges?");

    std::set<EdgeType> edgeSet;
    for (vtkm::IdComponent edgeIndex = 0; edgeIndex < numEdges; edgeIndex++)
    {
      EdgeType edge = vtkm::exec::CellEdgeLocalIndices(numPoints, edgeIndex, shape, workletProxy);
      VTKM_TEST_ASSERT(edge[0] >= 0, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[0] < numPoints, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[1] >= 0, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[1] < numPoints, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[0] != edge[1], "Degenerate edge.");
      MakeEdgeCononical(edge);
      VTKM_TEST_ASSERT(edge[0] < edge[1], "Internal test error: MakeEdgeCononical failed");
      VTKM_TEST_ASSERT(edgeSet.find(edge) == edgeSet.end(), "Found duplicate edge");
      edgeSet.insert(edge);

      vtkm::Id2 cononicalEdgeId =
        vtkm::exec::CellEdgeCononicalId(numPoints, edgeIndex, shape, pointIndexProxy, workletProxy);
      VTKM_TEST_ASSERT(cononicalEdgeId[0] > 0, "Not using global ids?");
      VTKM_TEST_ASSERT(cononicalEdgeId[0] < cononicalEdgeId[1], "Bad order.");
    }

    vtkm::IdComponent numFaces = vtkm::exec::CellFaceNumberOfFaces(shape, workletProxy);
    VTKM_TEST_ASSERT(numFaces > 0, "No faces?");

    std::set<EdgeType> edgesFoundInFaces;
    for (vtkm::IdComponent faceIndex = 0; faceIndex < numFaces; faceIndex++)
    {
      vtkm::VecCConst<vtkm::IdComponent> facePoints =
        vtkm::exec::CellFaceLocalIndices(faceIndex, shape, workletProxy);
      vtkm::IdComponent numPointsInFace = facePoints.GetNumberOfComponents();
      VTKM_TEST_ASSERT(numPointsInFace >= 3, "Face has fewer points than a triangle.");

      for (vtkm::IdComponent pointIndex = 0; pointIndex < numPointsInFace; pointIndex++)
      {
        VTKM_TEST_ASSERT(facePoints[pointIndex] >= 0, "Invalid point index for face.");
        VTKM_TEST_ASSERT(facePoints[pointIndex] <= numPoints, "Invalid point index for face.");
        EdgeType edge;
        if (pointIndex < numPointsInFace - 1)
        {
          edge = EdgeType(facePoints[pointIndex], facePoints[pointIndex + 1]);
        }
        else
        {
          edge = EdgeType(facePoints[0], facePoints[pointIndex]);
        }
        MakeEdgeCononical(edge);
        VTKM_TEST_ASSERT(edgeSet.find(edge) != edgeSet.end(), "Edge in face not in cell's edges");
        edgesFoundInFaces.insert(edge);
      }

      vtkm::Id3 cononicalFaceId =
        vtkm::exec::CellFaceCononicalId(faceIndex, shape, pointIndexProxy, workletProxy);
      VTKM_TEST_ASSERT(cononicalFaceId[0] > 0, "Not using global ids?");
      VTKM_TEST_ASSERT(cononicalFaceId[0] < cononicalFaceId[1], "Bad order.");
      VTKM_TEST_ASSERT(cononicalFaceId[1] < cononicalFaceId[2], "Bad order.");
    }
    VTKM_TEST_ASSERT(edgesFoundInFaces.size() == edgeSet.size(),
                     "Faces did not contain all edges in cell");
  }

  // Case of cells that have 2 dimensions (no faces)
  template <typename CellShapeTag>
  void DoTest(vtkm::IdComponent numPoints,
              CellShapeTag shape,
              vtkm::CellTopologicalDimensionsTag<2>) const
  {
    // Stuff to fake running in the execution environment.
    char messageBuffer[256];
    messageBuffer[0] = '\0';
    vtkm::exec::internal::ErrorMessageBuffer errorMessage(messageBuffer, 256);
    vtkm::exec::FunctorBase workletProxy;
    workletProxy.SetErrorMessageBuffer(errorMessage);

    std::vector<vtkm::Id> pointIndexProxyBuffer(static_cast<std::size_t>(numPoints));
    for (std::size_t index = 0; index < pointIndexProxyBuffer.size(); ++index)
    {
      pointIndexProxyBuffer[index] = static_cast<vtkm::Id>(1000000 - index);
    }
    vtkm::VecCConst<vtkm::Id> pointIndexProxy(&pointIndexProxyBuffer.at(0), numPoints);

    vtkm::IdComponent numEdges = vtkm::exec::CellEdgeNumberOfEdges(numPoints, shape, workletProxy);
    VTKM_TEST_ASSERT(numEdges == numPoints, "Polygons should have same number of points and edges");

    std::set<EdgeType> edgeSet;
    for (vtkm::IdComponent edgeIndex = 0; edgeIndex < numEdges; edgeIndex++)
    {
      EdgeType edge = vtkm::exec::CellEdgeLocalIndices(numPoints, edgeIndex, shape, workletProxy);
      VTKM_TEST_ASSERT(edge[0] >= 0, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[0] < numPoints, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[1] >= 0, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[1] < numPoints, "Bad index in edge.");
      VTKM_TEST_ASSERT(edge[0] != edge[1], "Degenerate edge.");
      MakeEdgeCononical(edge);
      VTKM_TEST_ASSERT(edge[0] < edge[1], "Internal test error: MakeEdgeCononical failed");
      VTKM_TEST_ASSERT(edgeSet.find(edge) == edgeSet.end(), "Found duplicate edge");
      edgeSet.insert(edge);

      vtkm::Id2 cononicalEdgeId =
        vtkm::exec::CellEdgeCononicalId(numPoints, edgeIndex, shape, pointIndexProxy, workletProxy);
      VTKM_TEST_ASSERT(cononicalEdgeId[0] > 0, "Not using global ids?");
      VTKM_TEST_ASSERT(cononicalEdgeId[0] < cononicalEdgeId[1], "Bad order.");
    }

    vtkm::IdComponent numFaces = vtkm::exec::CellFaceNumberOfFaces(shape, workletProxy);
    VTKM_TEST_ASSERT(numFaces == 0, "Non 3D shape should have no faces");
  }

  // Less important case of cells that have less than 2 dimensions
  // (no faces or edges)
  template <typename CellShapeTag, vtkm::IdComponent NumDimensions>
  void DoTest(vtkm::IdComponent numPoints,
              CellShapeTag shape,
              vtkm::CellTopologicalDimensionsTag<NumDimensions>) const
  {
    // Stuff to fake running in the execution environment.
    char messageBuffer[256];
    messageBuffer[0] = '\0';
    vtkm::exec::internal::ErrorMessageBuffer errorMessage(messageBuffer, 256);
    vtkm::exec::FunctorBase workletProxy;
    workletProxy.SetErrorMessageBuffer(errorMessage);

    vtkm::IdComponent numEdges = vtkm::exec::CellEdgeNumberOfEdges(numPoints, shape, workletProxy);
    VTKM_TEST_ASSERT(numEdges == 0, "0D or 1D shape should have no edges");

    vtkm::IdComponent numFaces = vtkm::exec::CellFaceNumberOfFaces(shape, workletProxy);
    VTKM_TEST_ASSERT(numFaces == 0, "Non 3D shape should have no faces");
  }

  template <typename CellShapeTag>
  void TryShapeWithNumPoints(vtkm::IdComponent numPoints, CellShapeTag) const
  {
    std::cout << "--- Test shape tag directly"
              << " (" << numPoints << " points)" << std::endl;
    this->DoTest(numPoints,
                 CellShapeTag(),
                 typename vtkm::CellTraits<CellShapeTag>::TopologicalDimensionsTag());

    std::cout << "--- Test generic shape tag"
              << " (" << numPoints << " points)" << std::endl;
    this->DoTest(numPoints,
                 vtkm::CellShapeTagGeneric(CellShapeTag::Id),
                 typename vtkm::CellTraits<CellShapeTag>::TopologicalDimensionsTag());
  }

  template <typename CellShapeTag>
  void operator()(CellShapeTag) const
  {
    this->TryShapeWithNumPoints(vtkm::CellTraits<CellShapeTag>::NUM_POINTS, CellShapeTag());
  }

  void operator()(vtkm::CellShapeTagPolygon) const
  {
    for (vtkm::IdComponent numPoints = 3; numPoints < 7; numPoints++)
    {
      this->TryShapeWithNumPoints(numPoints, vtkm::CellShapeTagPolygon());
    }
  }
};

void TestAllShapes()
{
  vtkm::testing::Testing::TryAllCellShapes(TestCellFacesFunctor());
}

} // anonymous namespace

int UnitTestCellEdgeFace(int, char* [])
{
  return vtkm::testing::Testing::Run(TestAllShapes);
}
