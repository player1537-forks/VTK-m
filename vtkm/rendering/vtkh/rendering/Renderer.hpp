//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_rendering_Renderer_h
#define vtkm_rendering_rendering_Renderer_h

#include <vtkm/rendering/vtkm_rendering_export.h>

#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/rendering/vtkh/rendering/Plot.h>
#include <vtkm/rendering/vtkh/compositing/Image.hpp>

#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>


/*
TODO:
Put Actor in.
remove Range (it's held in the actor)
remove Bounds (it's in Actor)
Finalize the global/local bounds in Actor and renderer.

*/

namespace vtkh
{

class Compositor;

class VTKM_RENDERING_EXPORT Renderer
{
public:
  typedef std::shared_ptr<vtkm::rendering::Mapper> vtkmMapperPtr;

  Renderer();
  virtual ~Renderer();
  virtual std::string GetName() const = 0;

//Actor stuff beg.
  virtual void SetInput(const vtkm::rendering::Actor& actor) {this->Actor = actor;}
  virtual void SetColorTable(const vtkm::cont::ColorTable &ct) {this->Actor.SetColorTable(ct);}
  vtkm::cont::ColorTable      GetColorTable() const {   return this->Actor.GetColorTable(); }
  std::string                 GetFieldName() const { return this->Actor.GetScalarFieldName(); }
  vtkm::Range                 GetScalarRange() const { return this->Actor.GetScalarRange(); }
  bool                        GetHasColorTable() const {   return this->HasColorTable; }
//Actor stuff end.

  virtual void SetShadingOn(bool vtkmNotUsed(on)) {}   // do nothing by default;
  virtual void Update();

  void AddPlot(vtkh::Plot &plot) {   this->Plots.push_back(plot);}
  void ClearPlots() { this->Plots.clear(); }

  void SetDoComposite(bool do_composite) {   this->DoComposite = do_composite; }
  void SetPlots(const std::vector<Plot> &plots) {   this->Plots = plots; }

  int                         GetNumberOfPlots() const { return static_cast<int>(this->Plots.size());}
  std::vector<Plot>         GetPlots() const {   return this->Plots;}


protected:
  // image related data with cinema support
  vtkh::Compositor                        *Compositor = nullptr;
  bool                                     DoComposite = true;
  bool                                     HasColorTable = true;
  vtkm::rendering::Actor Actor;
  vtkmMapperPtr Mapper;
  std::vector<vtkh::Plot> Plots; //used in Composite and DoExecute

  // methods
  virtual void PreExecute(); // override;
  virtual void PostExecute(); // override;
  virtual void DoExecute(); // override;

  void CheckForRequiredField(const std::string &field_name);

  virtual void Composite();
  void ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
};

} // namespace vtkh

#endif //vtkm_rendering_rendering_Renderer_h
