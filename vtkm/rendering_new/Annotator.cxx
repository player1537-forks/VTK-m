//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/rendering_new/Annotator.h>

namespace vtkm
{
namespace rendering_new
{

Annotator::Annotator(vtkm::rendering::Canvas& canvas,
                     vtkm::rendering::Camera& camera,
                     vtkm::Bounds bounds)
  : Canvas(canvas)
  , Camera(camera)
  , Bounds(bounds)
{
  this->Is3d = this->Camera.GetMode() == vtkm::rendering::Camera::Mode::ThreeD;
  this->WorldAnnotator = this->Canvas.CreateWorldAnnotator();
  // add defualt color bar positions
  vtkm::Bounds p1(vtkm::Range(0.84, 0.92), vtkm::Range(+0.1, +0.8), vtkm::Range(0, 0));
  vtkm::Bounds p2(vtkm::Range(0.84, 0.92), vtkm::Range(-0.8, -0.1), vtkm::Range(0, 0));
  vtkm::Bounds p3(vtkm::Range(-0.8, -0.72), vtkm::Range(+0.1, +0.8), vtkm::Range(0, 0));
  vtkm::Bounds p4(vtkm::Range(-0.8, -0.72), vtkm::Range(-0.8, -0.1), vtkm::Range(0, 0));

  this->ColorBarPos.push_back(p1);
  this->ColorBarPos.push_back(p2);
  this->ColorBarPos.push_back(p3);
  this->ColorBarPos.push_back(p4);
}

Annotator::~Annotator()
{
  delete this->WorldAnnotator;
}

void Annotator::RenderScreenAnnotations(const std::vector<std::string>& field_names,
                                        const std::vector<vtkm::Range>& ranges,
                                        const std::vector<vtkm::cont::ColorTable>& color_tables)
{
  this->Canvas.SetViewToScreenSpace(this->Camera, true);
  // currently we only support 4 color bars, so grab the first 4
  int num_bars = std::min(int(field_names.size()), 4);
  this->Canvas.BeginTextRenderingBatch();
  this->WorldAnnotator->BeginLineRenderingBatch();
  for (int i = 0; i < num_bars; ++i)
  {
    this->ColorBarAnnotation.SetRange(ranges[i], 5);
    this->ColorBarAnnotation.SetFieldName(field_names[i]);
    this->ColorBarAnnotation.SetPosition(this->ColorBarPos[i]);
    this->ColorBarAnnotation.SetColorTable(color_tables[i]);
    this->ColorBarAnnotation.Render(this->Camera, *this->WorldAnnotator, this->Canvas);
  }
  this->WorldAnnotator->EndLineRenderingBatch();
  this->Canvas.EndTextRenderingBatch();
}

void Annotator::RenderWorldAnnotations(vtkm::Vec<float, 3> axis_scale)
{
  if (!this->Is3d)
    return;
  this->Canvas.SetViewToWorldSpace(this->Camera, false);

  this->Canvas.BeginTextRenderingBatch();
  vtkm::Float64 xmin = this->Bounds.X.Min, xmax = this->Bounds.X.Max;
  vtkm::Float64 ymin = this->Bounds.Y.Min, ymax = this->Bounds.Y.Max;
  vtkm::Float64 zmin = this->Bounds.Z.Min, zmax = this->Bounds.Z.Max;
  vtkm::Float64 dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
  vtkm::Float64 size = vtkm::Sqrt(dx * dx + dy * dy + dz * dz);

  //TODO: get forground color
  this->WorldAnnotator->BeginLineRenderingBatch();
  this->BoxAnnotation.SetColor(this->Canvas.GetForegroundColor());
  this->BoxAnnotation.SetExtents(this->Bounds);
  this->BoxAnnotation.Render(this->Camera, *this->WorldAnnotator);
  vtkm::Vec<vtkm::Float32, 3> lookAt = this->Camera.GetLookAt();
  vtkm::Vec<vtkm::Float32, 3> position = this->Camera.GetPosition();
  bool xtest = lookAt[0] > position[0];
  bool ytest = lookAt[1] > position[1];
  bool ztest = lookAt[2] > position[2];
  this->WorldAnnotator->EndLineRenderingBatch();

  const bool outsideedges = true; // if false, do closesttriad
  if (outsideedges)
  {
    xtest = !xtest;
    //ytest = !ytest;
  }

  //vtkm::Float64 xrel = vtkm::Abs(dx) / size;
  //vtkm::Float64 yrel = vtkm::Abs(dy) / size;
  //vtkm::Float64 zrel = vtkm::Abs(dz) / size;
  float major_tick_size = size / 40.f;
  float minor_tick_size = size / 80.f;

  this->WorldAnnotator->BeginLineRenderingBatch();
  this->XAxisAnnotation.SetAxis(0);
  this->XAxisAnnotation.SetColor(this->Canvas.GetForegroundColor());
  this->XAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->XAxisAnnotation.SetWorldPosition(
    xmin, ytest ? ymin : ymax, ztest ? zmin : zmax, xmax, ytest ? ymin : ymax, ztest ? zmin : zmax);
  this->XAxisAnnotation.SetRange(xmin * axis_scale[0], xmax * axis_scale[0]);
  this->XAxisAnnotation.SetMajorTickSize(major_tick_size, 0);
  this->XAxisAnnotation.SetMinorTickSize(minor_tick_size, 0);
  this->XAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  this->XAxisAnnotation.SetMoreOrLessTickAdjustment(-1);
  //this->XAxisAnnotation.SetMoreOrLessTickAdjustment(xrel < .3 ? -1 : 0);
  this->XAxisAnnotation.Render(this->Camera, *this->WorldAnnotator, this->Canvas);

  this->YAxisAnnotation.SetAxis(1);
  this->YAxisAnnotation.SetColor(this->Canvas.GetForegroundColor());
  this->YAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->YAxisAnnotation.SetWorldPosition(
    xtest ? xmin : xmax, ymin, ztest ? zmin : zmax, xtest ? xmin : xmax, ymax, ztest ? zmin : zmax);
  this->YAxisAnnotation.SetRange(ymin * axis_scale[1], ymax * axis_scale[0]);
  this->YAxisAnnotation.SetMajorTickSize(major_tick_size, 0);
  this->YAxisAnnotation.SetMinorTickSize(minor_tick_size, 0);
  this->YAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  this->YAxisAnnotation.SetMoreOrLessTickAdjustment(-1);
  //this->YAxisAnnotation.SetMoreOrLessTickAdjustment(yrel < .3 ? -1 : 0);
  this->YAxisAnnotation.Render(this->Camera, *this->WorldAnnotator, this->Canvas);

  this->ZAxisAnnotation.SetAxis(2);
  this->ZAxisAnnotation.SetColor(this->Canvas.GetForegroundColor());
  this->ZAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->ZAxisAnnotation.SetWorldPosition(
    xtest ? xmin : xmax, ytest ? ymin : ymax, zmin, xtest ? xmin : xmax, ytest ? ymin : ymax, zmax);
  this->ZAxisAnnotation.SetRange(zmin * axis_scale[2], zmax * axis_scale[2]);
  this->ZAxisAnnotation.SetMajorTickSize(major_tick_size, 0);
  this->ZAxisAnnotation.SetMinorTickSize(minor_tick_size, 0);
  this->ZAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  //this->ZAxisAnnotation.SetMoreOrLessTickAdjustment(zrel < .3 ? -1 : 0);
  this->ZAxisAnnotation.SetMoreOrLessTickAdjustment(-1);
  this->ZAxisAnnotation.Render(this->Camera, *this->WorldAnnotator, this->Canvas);
  this->WorldAnnotator->EndLineRenderingBatch();

  this->Canvas.EndTextRenderingBatch();
}

}
} //namespace vtkm::rendering_new
