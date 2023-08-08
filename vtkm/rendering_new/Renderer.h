//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_vtkh_rendering_Renderer_h
#define vtkm_rendering_vtkh_rendering_Renderer_h

#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering_new/Plot.h>
#include <vtkm/rendering_new/compositing/Image.h>
#include <vtkm/rendering_new/vtkm_rendering_new_export.h>

#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

#include <vector>


/*
TODO:
Put Actor in.
remove Range (it's held in the actor)
remove Bounds (it's in Actor)
Finalize the global/local bounds in Actor and renderer.

*/

namespace vtkm
{
namespace rendering_new
{

class Compositor;

class VTKM_RENDERING_NEW_EXPORT Renderer
{
public:
  typedef std::shared_ptr<vtkm::rendering::Mapper> vtkmMapperPtr;

  Renderer();
  virtual ~Renderer();
  virtual std::string GetName() const = 0;

  virtual void SetInput(const vtkm::rendering::Actor& actor) { this->Actor = actor; }
  virtual void SetColorTable(const vtkm::cont::ColorTable& ct) { this->Actor.SetColorTable(ct); }
  vtkm::cont::ColorTable GetColorTable() const { return this->Actor.GetColorTable(); }
  std::string GetFieldName() const { return this->Actor.GetScalarFieldName(); }
  vtkm::Range GetScalarRange() const { return this->Actor.GetScalarRange(); }
  bool GetHasColorTable() const { return this->HasColorTable; }

  virtual void SetShadingOn(bool vtkmNotUsed(on)) {} // do nothing by default;
  virtual void Update(vtkm::rendering_new::Plot& plot);

  void SetDoComposite(bool do_composite) { this->DoComposite = do_composite; }

protected:
  // image related data with cinema support
  vtkm::rendering_new::Compositor* Compositor = nullptr;
  bool DoComposite = true;
  bool HasColorTable = true;
  vtkm::rendering::Actor Actor;
  vtkmMapperPtr Mapper;

  // methods
  virtual void PreExecute(vtkm::rendering_new::Plot& plot);
  virtual void PostExecute(vtkm::rendering_new::Plot& plot);
  virtual void DoExecute(vtkm::rendering_new::Plot& plot);

  void CheckForRequiredField(const std::string& field_name);

  virtual void Composite(vtkm::rendering_new::Plot& plot);
  void ImageToCanvas(Image& image, vtkm::rendering::Canvas& canvas, bool get_depth);
};

}
} // namespace vtkm::rendering_new

#endif //vtkm_rendering_vtkh_rendering_Renderer_h
