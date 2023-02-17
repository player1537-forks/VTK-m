//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/rendering/View3D.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/thirdparty/diy/diy.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace rendering
{

View3D::View3D(const vtkm::rendering::Scene& scene,
               const vtkm::rendering::Mapper& mapper,
               const vtkm::rendering::Canvas& canvas,
               const vtkm::rendering::Color& backgroundColor,
               const vtkm::rendering::Color& foregroundColor)
  : View(scene, mapper, canvas, backgroundColor, foregroundColor)
{
}

View3D::View3D(const vtkm::rendering::Scene& scene,
               const vtkm::rendering::Mapper& mapper,
               const vtkm::rendering::Canvas& canvas,
               const vtkm::rendering::Camera& camera,
               const vtkm::rendering::Color& backgroundColor,
               const vtkm::rendering::Color& foregroundColor)
  : View(scene, mapper, canvas, camera, backgroundColor, foregroundColor)
{
}

View3D::~View3D() {}

void View3D::Paint()
{
  this->GetCanvas().Clear();
  this->RenderAnnotations();
  this->GetScene().Render(this->GetMapper(), this->GetCanvas(), this->GetCamera());

#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.size() == 1)
    return;

  this->Compositor.SetCompositeMode(vtkm::rendering::compositing::Compositor::Z_BUFFER_SURFACE);
  //volume render
  this->Compositor.SetCompositeMode(vtkm::rendering::compositing::Compositor::VIS_ORDER_BLEND);
  /*
  auto colors = (this->GetCanvas().GetColorBuffer().WritePortal().GetArray())[0][0];
  auto depths = (this->GetCanvas().GetDepthBuffer().WritePortal().GetArray());
  //auto colors = &GetVTKMPointer(this->GetCanvas().GetColorBuffer())[0][0];
  //auto depths = GetVTKMPointer(this->GetCanvas().GetDepthBuffer());
  */
  this->Compositor.AddImage(this->GetCanvas());
  auto result = this->Compositor.Composite();

  //Rank 0 has the composited result, so put it into the Canvas.
  if (comm.rank() == 0)
  {
    this->GetCanvas().Clear();
    auto colors = this->GetCanvas().GetColorBuffer();
    auto depths = this->GetCanvas().GetDepthBuffer();

    int size = this->GetCanvas().GetWidth() * this->GetCanvas().GetHeight();
    for (int i = 0; i < size; i++)
    {
      const int offset = i * 4;
      vtkm::Vec4f_32 rgba;
      for (int j = 0; j < 4; j++)
        rgba[j] = static_cast<vtkm::FloatDefault>(result.Pixels[offset + j] / 255.f);

      colors.WritePortal().Set(i, rgba);
      depths.WritePortal().Set(i, result.Depths[i]);
    }
  }
#endif
}

void View3D::RenderScreenAnnotations()
{
  vtkm::Range scalarRange;

  int numActors = this->GetScene().GetNumberOfActors();
  if (numActors > 0)
    scalarRange = this->GetScene().GetActor(0).GetScalarRange();

  int totNumActors = numActors;

  /*
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  vtkm::Float64 minVal = scalarRange.Min, maxVal = scalarRange.Max;

  MPI_Comm mpiComm = vtkmdiy::mpi::mpi_cast(comm.handle());
  int totNumActors = 0;
  vtkm::Float64 minVal_res = 0, maxVal_res = 0;
  MPI_Reduce(&numActors, &totNumActors, 1, MPI_INT, MPI_SUM, 0, mpiComm);
  MPI_Reduce(&minVal, &minVal_res, 1, MPI_DOUBLE, MPI_MIN, 0, mpiComm);
  MPI_Reduce(&maxVal, &maxVal_res, 1, MPI_DOUBLE, MPI_MAX, 0, mpiComm);
  if (comm.rank() != 0)
    return;

  scalarRange.Min = minVal_res;
  scalarRange.Max = maxVal_res;
#endif

  std::cout<<"totNumActors= "<<totNumActors<<" range= "<<scalarRange<<std::endl;

  //DRP
  //This assumes that rank 0 has an actor!!
  */

  if (totNumActors > 0)
  {
    this->GetCanvas().BeginTextRenderingBatch();
    this->GetWorldAnnotator().BeginLineRenderingBatch();
    //this->ColorBarAnnotation.SetAxisColor(vtkm::rendering::Color(1,1,1));
    this->ColorBarAnnotation.SetFieldName(this->GetScene().GetActor(0).GetScalarField().GetName());
    this->ColorBarAnnotation.SetRange(scalarRange, 5);
    this->ColorBarAnnotation.SetColorTable(this->GetScene().GetActor(0).GetColorTable());
    this->ColorBarAnnotation.Render(
      this->GetCamera(), this->GetWorldAnnotator(), this->GetCanvas());
    this->GetWorldAnnotator().EndLineRenderingBatch();
    this->GetCanvas().EndTextRenderingBatch();
  }
}

void View3D::RenderWorldAnnotations()
{
  vtkm::Bounds bounds = this->GetScene().GetSpatialBounds();

#ifdef VTKM_ENABLE_MPI
  //For parallel, get the collective bounds.
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  vtkm::Float64 mins[3], maxs[3], mins_res[3], maxs_res[3];
  mins[0] = bounds.X.Min;
  mins[1] = bounds.Y.Min;
  mins[2] = bounds.Z.Min;
  maxs[0] = bounds.X.Max;
  maxs[1] = bounds.Y.Max;
  maxs[2] = bounds.Z.Max;

  //DRP
  //what if a scene has NO actors??

  MPI_Comm mpiComm = vtkmdiy::mpi::mpi_cast(comm.handle());
  MPI_Reduce(mins, mins_res, 3, MPI_DOUBLE, MPI_MIN, 0, mpiComm);
  MPI_Reduce(maxs, maxs_res, 3, MPI_DOUBLE, MPI_MAX, 0, mpiComm);
  if (comm.rank() != 0)
    return;

  bounds.X.Min = mins_res[0];
  bounds.Y.Min = mins_res[1];
  bounds.Z.Min = mins_res[2];
  bounds.X.Max = maxs_res[0];
  bounds.Y.Max = maxs_res[1];
  bounds.Z.Max = maxs_res[2];

#endif

  this->GetCanvas().BeginTextRenderingBatch();
  vtkm::Float64 xmin = bounds.X.Min, xmax = bounds.X.Max;
  vtkm::Float64 ymin = bounds.Y.Min, ymax = bounds.Y.Max;
  vtkm::Float64 zmin = bounds.Z.Min, zmax = bounds.Z.Max;
  vtkm::Float64 dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
  vtkm::Float64 size = vtkm::Sqrt(dx * dx + dy * dy + dz * dz);

  this->GetWorldAnnotator().BeginLineRenderingBatch();
  this->BoxAnnotation.SetColor(Color(.5f, .5f, .5f));
  this->BoxAnnotation.SetExtents(this->GetScene().GetSpatialBounds());
  this->BoxAnnotation.Render(this->GetCamera(), this->GetWorldAnnotator());
  this->GetWorldAnnotator().EndLineRenderingBatch();

  vtkm::Vec3f_32 lookAt = this->GetCamera().GetLookAt();
  vtkm::Vec3f_32 position = this->GetCamera().GetPosition();
  bool xtest = lookAt[0] > position[0];
  bool ytest = lookAt[1] > position[1];
  bool ztest = lookAt[2] > position[2];

  const bool outsideedges = true; // if false, do closesttriad
  if (outsideedges)
  {
    xtest = !xtest;
    //ytest = !ytest;
  }

  vtkm::Float64 xrel = vtkm::Abs(dx) / size;
  vtkm::Float64 yrel = vtkm::Abs(dy) / size;
  vtkm::Float64 zrel = vtkm::Abs(dz) / size;

  this->GetWorldAnnotator().BeginLineRenderingBatch();
  this->XAxisAnnotation.SetAxis(0);
  this->XAxisAnnotation.SetColor(AxisColor);
  this->XAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->XAxisAnnotation.SetWorldPosition(
    xmin, ytest ? ymin : ymax, ztest ? zmin : zmax, xmax, ytest ? ymin : ymax, ztest ? zmin : zmax);
  this->XAxisAnnotation.SetRange(xmin, xmax);
  this->XAxisAnnotation.SetMajorTickSize(size / 40.f, 0);
  this->XAxisAnnotation.SetMinorTickSize(size / 80.f, 0);
  this->XAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  this->XAxisAnnotation.SetMoreOrLessTickAdjustment(xrel < .3 ? -1 : 0);
  this->XAxisAnnotation.Render(this->GetCamera(), this->GetWorldAnnotator(), this->GetCanvas());

  this->YAxisAnnotation.SetAxis(1);
  this->YAxisAnnotation.SetColor(AxisColor);
  this->YAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->YAxisAnnotation.SetWorldPosition(
    xtest ? xmin : xmax, ymin, ztest ? zmin : zmax, xtest ? xmin : xmax, ymax, ztest ? zmin : zmax);
  this->YAxisAnnotation.SetRange(ymin, ymax);
  this->YAxisAnnotation.SetMajorTickSize(size / 40.f, 0);
  this->YAxisAnnotation.SetMinorTickSize(size / 80.f, 0);
  this->YAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  this->YAxisAnnotation.SetMoreOrLessTickAdjustment(yrel < .3 ? -1 : 0);
  this->YAxisAnnotation.Render(this->GetCamera(), this->GetWorldAnnotator(), this->GetCanvas());

  this->ZAxisAnnotation.SetAxis(2);
  this->ZAxisAnnotation.SetColor(AxisColor);
  this->ZAxisAnnotation.SetTickInvert(xtest, ytest, ztest);
  this->ZAxisAnnotation.SetWorldPosition(
    xtest ? xmin : xmax, ytest ? ymin : ymax, zmin, xtest ? xmin : xmax, ytest ? ymin : ymax, zmax);
  this->ZAxisAnnotation.SetRange(zmin, zmax);
  this->ZAxisAnnotation.SetMajorTickSize(size / 40.f, 0);
  this->ZAxisAnnotation.SetMinorTickSize(size / 80.f, 0);
  this->ZAxisAnnotation.SetLabelFontOffset(vtkm::Float32(size / 15.f));
  this->ZAxisAnnotation.SetMoreOrLessTickAdjustment(zrel < .3 ? -1 : 0);
  this->ZAxisAnnotation.Render(this->GetCamera(), this->GetWorldAnnotator(), this->GetCanvas());
  this->GetWorldAnnotator().EndLineRenderingBatch();

  this->GetCanvas().EndTextRenderingBatch();
}
}
} // namespace vtkm::rendering
