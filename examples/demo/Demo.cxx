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
#include <vtkm/cont/Initialize.h>
#include <vtkm/source/Tangle.h>

#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/MapperVolume.h>
#include <vtkm/rendering/MapperWireframer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>

#include <vtkm/rendering/vtkh/rendering/RayTracer.hpp>
#include <vtkm/rendering/vtkh/rendering/ScalarRenderer.hpp>
#include <vtkm/rendering/vtkh/rendering/Scene.hpp>
#include <vtkm/rendering/vtkh/rendering/VolumeRenderer.hpp>

#include <vtkm/filter/contour/Contour.h>

// This example creates a simple data set and uses VTK-m's rendering engine to render an image and
// write that image to a file. It then computes an isosurface on the input data set and renders
// this output data set in a separate image file

using vtkm::rendering::CanvasRayTracer;
using vtkm::rendering::MapperRayTracer;
using vtkm::rendering::MapperVolume;
using vtkm::rendering::MapperWireframer;

int main(int argc, char* argv[])
{
  vtkm::cont::Initialize(argc, argv, vtkm::cont::InitializeOptions::Strict);

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  vtkm::source::Tangle tangle;
  tangle.SetPointDimensions({ 50, 50, 50 });
  if (comm.rank() == 0)
    tangle.SetOrigin({ 0.f, 0.f, 0.f });
  else
    tangle.SetOrigin({ 1.f, 0.f, 0.f });
  vtkm::cont::DataSet tangleData = tangle.Execute();
  std::string fieldName = "tangle";
  vtkm::cont::PartitionedDataSet pds(tangleData);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  // Set up a camera for rendering the input data
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.Azimuth(-30.f);
  camera.Elevation(-30.f);
  /*
  camera.SetLookAt(vtkm::Vec3f_32(0.5, 0.5, 0.5));
  camera.SetViewUp(vtkm::make_Vec(0.f, 1.f, 0.f));
  camera.SetClippingRange(1.f, 10.f);
  camera.SetFieldOfView(60.f);
  camera.SetPosition(vtkm::Vec3f_32(1.5, 1.5, 1.5));
  */


  vtkm::cont::ColorTable colorTable("inferno");

  vtkh::VolumeRenderer tracer;
  tracer.SetColorTable(colorTable);
  tracer.SetInput(&pds);
  tracer.SetField(fieldName);

  std::string imgName = "volume";
  vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, imgName);
  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();

  if (1)
  {
    // Compute an isosurface:
    vtkm::filter::contour::Contour filter;
    // [min, max] of the tangle field is [-0.887, 24.46]:
    filter.SetIsoValue(3.0);
    filter.SetActiveField(fieldName);
    vtkm::cont::DataSet isoData = filter.Execute(tangleData);

    pds = vtkm::cont::PartitionedDataSet(isoData);
    vtkm::Bounds bounds = pds.GetGlobalBounds();
    camera.ResetToBounds(bounds);

    vtkh::RayTracer rayTracer;
    rayTracer.SetInput(&pds);
    rayTracer.SetField(fieldName);
    imgName = "rayTracer";
    vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, imgName);

    vtkh::Scene scene;
    scene.AddRender(render);
    scene.AddRenderer(&rayTracer);
    scene.Render();
  }

  /*
  // Background color:
  vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);
  vtkm::rendering::Actor actor(tangleData.GetCellSet(),
                               tangleData.GetCoordinateSystem(),
                               tangleData.GetField(fieldName),
                               colorTable);
  vtkm::rendering::Scene scene;
  scene.AddActor(actor);
  // 2048x2048 pixels in the canvas:
  CanvasRayTracer canvas(2048, 2048);
  // Create a view and use it to render the input data using OS Mesa

  vtkm::rendering::View3D view(scene, MapperVolume(), canvas, camera, bg);
  view.Paint();
  view.SaveAs("volume.png");

  // Compute an isosurface:
  vtkm::filter::contour::Contour filter;
  // [min, max] of the tangle field is [-0.887, 24.46]:
  filter.SetIsoValue(3.0);
  filter.SetActiveField(fieldName);
  vtkm::cont::DataSet isoData = filter.Execute(tangleData);

  // Render a separate image with the output isosurface
  vtkm::rendering::Actor isoActor(
    isoData.GetCellSet(), isoData.GetCoordinateSystem(), isoData.GetField(fieldName), colorTable);
  // By default, the actor will automatically scale the scalar range of the color table to match
  // that of the data. However, we are coloring by the scalar that we just extracted a contour
  // from, so we want the scalar range to match that of the previous image.
  isoActor.SetScalarRange(actor.GetScalarRange());
  vtkm::rendering::Scene isoScene;
  isoScene.AddActor(isoActor);

  // Wireframe surface:
  vtkm::rendering::View3D isoView(isoScene, MapperWireframer(), canvas, camera, bg);
  isoView.Paint();
  isoView.SaveAs("isosurface_wireframer.png");

  // Smooth surface:
  vtkm::rendering::View3D solidView(isoScene, MapperRayTracer(), canvas, camera, bg);
  solidView.Paint();
  solidView.SaveAs("isosurface_raytracer.png");
  */

  return 0;
}
