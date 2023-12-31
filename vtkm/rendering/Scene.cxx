//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering/Scene.h>

#include <vector>

namespace vtkm
{
namespace rendering
{

struct Scene::InternalsType
{
  std::vector<vtkm::rendering::Actor> Actors;
};

Scene::Scene()
  : Internals(new InternalsType)
{
}

void Scene::AddActor(vtkm::rendering::Actor actor)
{
  this->Internals->Actors.push_back(std::move(actor));
}

const vtkm::rendering::Actor& Scene::GetActor(vtkm::IdComponent index) const
{
  return this->Internals->Actors[static_cast<std::size_t>(index)];
}

vtkm::IdComponent Scene::GetNumberOfActors() const
{
  return static_cast<vtkm::IdComponent>(this->Internals->Actors.size());
}

void Scene::Render(vtkm::rendering::Mapper& mapper,
                   vtkm::rendering::Canvas& canvas,
                   const vtkm::rendering::Camera& camera) const
{
  for (const auto& actor : this->Internals->Actors)
  {
    actor.Render(mapper, canvas, camera);
  }
}

vtkm::Bounds Scene::GetSpatialBounds() const
{
  vtkm::Bounds bounds;
  for (const auto& actor : this->Internals->Actors)
  {
    bounds.Include(actor.GetSpatialBounds());
  }

  return bounds;
}
}
} // namespace vtkm::rendering
