//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_exec_CellLocatorBoundingIntervalHierarchy_h
#define vtk_m_exec_CellLocatorBoundingIntervalHierarchy_h

#include <vtkm/TopologyElementTag.h>
#include <vtkm/VecFromPortalPermute.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/exec/CellInside.h>
#include <vtkm/exec/ParametricCoordinates.h>

namespace vtkm
{
namespace exec
{




struct CellLocatorBoundingIntervalHierarchyNode
{
#if defined(VTKM_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnested-anon-types"
#endif // gcc || clang
  vtkm::IdComponent Dimension;
  vtkm::Id ParentIndex;
  vtkm::Id ChildIndex;
  union {
    struct
    {
      vtkm::FloatDefault LMax;
      vtkm::FloatDefault RMin;
    } Node;
    struct
    {
      vtkm::Id Start;
      vtkm::Id Size;
    } Leaf;
  };
#if defined(VTKM_CLANG)
#pragma GCC diagnostic pop
#endif // gcc || clang

  VTKM_EXEC_CONT
  CellLocatorBoundingIntervalHierarchyNode()
    : Dimension()
    , ParentIndex()
    , ChildIndex()
    , Node{ 0, 0 }
  {
  }
}; // struct CellLocatorBoundingIntervalHierarchyNode

template <typename CellSetType>
class VTKM_ALWAYS_EXPORT CellLocatorBoundingIntervalHierarchy
{
  using NodeArrayHandle =
    vtkm::cont::ArrayHandle<vtkm::exec::CellLocatorBoundingIntervalHierarchyNode>;
  using CellIdArrayHandle = vtkm::cont::ArrayHandle<vtkm::Id>;

public:
  VTKM_CONT
  CellLocatorBoundingIntervalHierarchy(
    const NodeArrayHandle& nodes,
    const CellIdArrayHandle& cellIds,
    const CellSetType& cellSet,
    const vtkm::cont::CoordinateSystem::MultiplexerArrayType& coords,
    vtkm::cont::DeviceAdapterId device,
    vtkm::cont::Token& token)
    : Nodes(nodes.PrepareForInput(device, token))
    , CellIds(cellIds.PrepareForInput(device, token))
    , CellSet(cellSet.PrepareForInput(device, VisitType(), IncidentType(), token))
    , Coords(coords.PrepareForInput(device, token))
  {
  }


  VTKM_EXEC
  vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                           vtkm::Id& cellId,
                           vtkm::Vec3f& parametric) const
  {
    cellId = -1;
    vtkm::Id nodeIndex = 0;
    FindCellState state = FindCellState::EnterNode;

    while ((cellId < 0) && !((nodeIndex == 0) && (state == FindCellState::AscendFromNode)))
    {
      switch (state)
      {
        case FindCellState::EnterNode:
          VTKM_RETURN_ON_ERROR(this->EnterNode(state, point, cellId, nodeIndex, parametric));
          break;
        case FindCellState::AscendFromNode:
          this->AscendFromNode(state, nodeIndex);
          break;
        case FindCellState::DescendLeftChild:
          this->DescendLeftChild(state, point, nodeIndex);
          break;
        case FindCellState::DescendRightChild:
          this->DescendRightChild(state, point, nodeIndex);
          break;
      }
    }

    if (cellId >= 0)
    {
      return vtkm::ErrorCode::Success;
    }
    else
    {
      printf("BIH CellNotFound:: %d\n", __LINE__);
      return vtkm::ErrorCode::CellNotFound;
    }
  }

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorBoundingIntervalHierarchy* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorBoundingIntervalHierarchy* operator->() const { return this; }

private:
  enum struct FindCellState
  {
    EnterNode,
    AscendFromNode,
    DescendLeftChild,
    DescendRightChild
  };

  VTKM_EXEC
  vtkm::ErrorCode EnterNode(FindCellState& state,
                            const vtkm::Vec3f& point,
                            vtkm::Id& cellId,
                            vtkm::Id nodeIndex,
                            vtkm::Vec3f& parametric) const
  {
    VTKM_ASSERT(state == FindCellState::EnterNode);

    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& node = this->Nodes.Get(nodeIndex);

    if (node.ChildIndex < 0)
    {
      // In a leaf node. Look for a containing cell.
      VTKM_RETURN_ON_ERROR(this->FindInLeaf(point, parametric, node, cellId));
      state = FindCellState::AscendFromNode;
    }
    else
    {
      state = FindCellState::DescendLeftChild;
    }
    return vtkm::ErrorCode::Success;
  }

  VTKM_EXEC
  void AscendFromNode(FindCellState& state, vtkm::Id& nodeIndex) const
  {
    VTKM_ASSERT(state == FindCellState::AscendFromNode);

    vtkm::Id childNodeIndex = nodeIndex;
    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& childNode =
      this->Nodes.Get(childNodeIndex);
    nodeIndex = childNode.ParentIndex;
    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& parentNode =
      this->Nodes.Get(nodeIndex);

    if (parentNode.ChildIndex == childNodeIndex)
    {
      // Ascending from left child. Descend into the right child.
      state = FindCellState::DescendRightChild;
    }
    else
    {
      VTKM_ASSERT(parentNode.ChildIndex + 1 == childNodeIndex);
      // Ascending from right child. Ascend again. (Don't need to change state.)
    }
  }

  VTKM_EXEC
  void DescendLeftChild(FindCellState& state, const vtkm::Vec3f& point, vtkm::Id& nodeIndex) const
  {
    VTKM_ASSERT(state == FindCellState::DescendLeftChild);

    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& node = this->Nodes.Get(nodeIndex);
    const vtkm::FloatDefault& coordinate = point[node.Dimension];
    if (coordinate <= node.Node.LMax)
    {
      // Left child does contain the point. Do the actual descent.
      nodeIndex = node.ChildIndex;
      state = FindCellState::EnterNode;
    }
    else
    {
      // Left child does not contain the point. Skip to the right child.
      state = FindCellState::DescendRightChild;
    }
  }

  VTKM_EXEC
  void DescendRightChild(FindCellState& state, const vtkm::Vec3f& point, vtkm::Id& nodeIndex) const
  {
    VTKM_ASSERT(state == FindCellState::DescendRightChild);

    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& node = this->Nodes.Get(nodeIndex);
    const vtkm::FloatDefault& coordinate = point[node.Dimension];
    if (coordinate >= node.Node.RMin)
    {
      // Right child does contain the point. Do the actual descent.
      nodeIndex = node.ChildIndex + 1;
      state = FindCellState::EnterNode;
    }
    else
    {
      // Right child does not contain the point. Skip to ascent
      state = FindCellState::AscendFromNode;
    }
  }

  VTKM_EXEC vtkm::ErrorCode FindInLeaf(
    const vtkm::Vec3f& point,
    vtkm::Vec3f& parametric,
    const vtkm::exec::CellLocatorBoundingIntervalHierarchyNode& node,
    vtkm::Id& containingCellId) const
  {
    using IndicesType = typename CellSetPortal::IndicesType;
    for (vtkm::Id i = node.Leaf.Start; i < node.Leaf.Start + node.Leaf.Size; ++i)
    {
      vtkm::Id cellId = this->CellIds.Get(i);
      IndicesType cellPointIndices = this->CellSet.GetIndices(cellId);
      vtkm::VecFromPortalPermute<IndicesType, CoordsPortal> cellPoints(&cellPointIndices,
                                                                       this->Coords);
      bool found;
      VTKM_RETURN_ON_ERROR(this->IsPointInCell(
        point, parametric, this->CellSet.GetCellShape(cellId), cellPoints, found));
      if (found)
      {
        containingCellId = cellId;
        return vtkm::ErrorCode::Success;
      }
    }
    containingCellId = -1;
    return vtkm::ErrorCode::Success;
  }

  template <typename CoordsType, typename CellShapeTag>
  VTKM_EXEC static vtkm::ErrorCode IsPointInCell(const vtkm::Vec3f& point,
                                                 vtkm::Vec3f& parametric,
                                                 CellShapeTag cellShape,
                                                 const CoordsType& cellPoints,
                                                 bool& isInside)
  {
    isInside = false;
    VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
      cellPoints, point, cellShape, parametric));
    isInside = vtkm::exec::CellInside(parametric, cellShape);
    return vtkm::ErrorCode::Success;
  }

  using VisitType = vtkm::TopologyElementTagCell;
  using IncidentType = vtkm::TopologyElementTagPoint;
  using NodePortal = typename NodeArrayHandle::ReadPortalType;
  using CellIdPortal = typename CellIdArrayHandle::ReadPortalType;
  using CellSetPortal =
    typename CellSetType::template ExecConnectivityType<VisitType, IncidentType>;
  using CoordsPortal = typename vtkm::cont::CoordinateSystem::MultiplexerArrayType::ReadPortalType;

  NodePortal Nodes;
  CellIdPortal CellIds;
  CellSetPortal CellSet;
  CoordsPortal Coords;
}; // class CellLocatorBoundingIntervalHierarchy

} // namespace exec

} // namespace vtkm

#endif //vtk_m_exec_CellLocatorBoundingIntervalHierarchy_h
