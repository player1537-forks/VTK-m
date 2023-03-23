//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/testing/Testing.h>
#include <vtkm/rendering/compositing/Compositor.h>
#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/compositing/ImageCompositor.h>
#include <vtkm/rendering/testing/t_vtkm_test_utils.h>
#include <vtkm/source/Tangle.h>

#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperConnectivity.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/MapperVolume.h>
#include <vtkm/rendering/MapperWireframer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/testing/RenderTest.h>

#include <vtkm/io/VTKDataSetReader.h>

#include <vtkm/filter/clean_grid/CleanGrid.h>

namespace
{

template <typename T>
T* GetVTKMPointer(vtkm::cont::ArrayHandle<T>& handle)
{
  return handle.WritePortal().GetArray();
}


vtkm::cont::DataSet ReadDS(int rank)
{

  std::string vtkFile;

  vtkm::io::VTKDataSetReader reader(vtkFile);
}

#if 0
vtkm::rendering::compositing::Image ConstImage(const std::size_t& width,
                                               const std::size_t& height,
                                               const vtkm::Vec4f& rgba,
                                               const vtkm::FloatDefault& depth)
{
  auto numPix = width * height;
  std::vector<vtkm::FloatDefault> rgbaVals(numPix * 4);
  std::vector<vtkm::FloatDefault> depthVals(numPix, depth);

  for (std::size_t i = 0; i < numPix; i++)
  {
    rgbaVals[i * 4 + 0] = rgba[0];
    rgbaVals[i * 4 + 1] = rgba[1];
    rgbaVals[i * 4 + 2] = rgba[2];
    rgbaVals[i * 4 + 3] = rgba[3];
  }

  vtkm::rendering::compositing::Image img(vtkm::Bounds(0, width, 0, height, 0, 1));
  img.Init(rgbaVals.data(), depthVals.data(), width, height);

  return img;
}
#endif

void TestImageComposite()
{
#if 0
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  std::size_t width = 4, height = 4;

  //res is the background, initially black.
  auto img0 = ConstImage(width, height, { 1, 0, 0, 1 }, 1.0);

  auto img1 = ConstImage(width, height, { 0, 1, 1, .5 }, 0.5);

  vtkm::rendering::compositing::Compositor compositor;

  compositor.SetCompositeMode(vtkm::rendering::compositing::Compositor::Z_BUFFER_SURFACE);
  vtkm::rendering::compositing::Image img;

  if (comm.rank() == 0)
    img = ConstImage(width, height, { 1, 0, 0, 1 }, 1.0);
  else
    img = ConstImage(width, height, { 0, 1, 1, .5 }, 0.5);

  compositor.AddImage(img.m_pixels.data(), img.m_depths.data(), width, height);

  auto res = compositor.Composite();

  //vtkm::rendering::compositing::ImageCompositor imgCompositor;
  //compositor.ZBufferComposite(res, img);
  //compositor.Blend(res, img);

  if (comm.rank() == 0)
  {
    for (int i = 0; i < width * height; i++)
    {
      std::cout << i << ": ";
      std::cout << (int)res.m_pixels[i * 4 + 0] << " ";
      std::cout << (int)res.m_pixels[i * 4 + 1] << " ";
      std::cout << (int)res.m_pixels[i * 4 + 2] << " ";
      std::cout << (int)res.m_pixels[i * 4 + 3] << " ";
      std::cout << res.m_depths[i] << std::endl;
    }
  }
#endif
}

void TestRenderComposite()
{
  using vtkm::rendering::CanvasRayTracer;
  using vtkm::rendering::MapperRayTracer;
  using vtkm::rendering::MapperVolume;
  using vtkm::rendering::MapperWireframer;

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  int numBlocks = comm.size() * 1;
  int rank = comm.rank();
  int dsPerRank = 2;

  vtkm::rendering::Camera camera;
  camera.SetLookAt(vtkm::Vec3f_32(1.0, 0.5, 0.5));
  camera.SetViewUp(vtkm::make_Vec(0.f, 1.f, 0.f));
  camera.SetClippingRange(1.f, 10.f);
  camera.SetFieldOfView(60.f);
  camera.SetPosition(vtkm::Vec3f_32(-2, 1.75, 1.75));
  vtkm::cont::ColorTable colorTable("inferno");

  // Background color:
  vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);
  vtkm::rendering::Scene scene;
  int width = 512, height = 512;
  CanvasRayTracer canvas(width, height);

  for (int i = 0; i < dsPerRank; i++)
  {
    //Create a sequence of datasets along the X direction.
    std::string fieldName = "tangle";
    vtkm::source::Tangle tangle;
    vtkm::Vec3f pt(rank * dsPerRank + i, 0, 0);
    if (rank == 1)
      std::cout << "PT= " << pt << std::endl;
    tangle.SetPointDimensions({ 50, 50, 50 });
    tangle.SetOrigin(pt);
    vtkm::cont::DataSet ds = tangle.Execute();

    vtkm::rendering::Actor actor(
      ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldName), colorTable);
    scene.AddActor(actor);
  }

  vtkm::rendering::View3D view(scene, MapperRayTracer(), canvas, camera, bg);
  view.Paint();

  canvas.SaveAs("result.png");

  /*
  auto colors = &GetVTKMPointer(canvas.GetColorBuffer())[0][0];
  auto depths = GetVTKMPointer(canvas.GetDepthBuffer());

  vtkm::rendering::compositing::Compositor compositor;
  compositor.AddImage(colors, depths, width, height);
  auto res = compositor.Composite();
  if (comm.rank() == 0)
  {
    res.Save("RESULT.png", { "" });
  }
  */
}

void TestVolumeRenderComposite(bool unstructured)
{
  using vtkm::rendering::CanvasRayTracer;
  using vtkm::rendering::MapperConnectivity;
  using vtkm::rendering::MapperVolume;

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  int numBlocks = comm.size() * 1;
  int rank = comm.rank();
  int dsPerRank = 1;

  vtkm::rendering::Camera camera;
  camera.SetLookAt(vtkm::Vec3f_32(1.0, 0.5, 0.5));
  camera.SetViewUp(vtkm::make_Vec(0.f, 1.f, 0.f));
  camera.SetClippingRange(1.f, 10.f);
  camera.SetFieldOfView(60.f);
  camera.SetPosition(vtkm::Vec3f_32(-2, 1.75, 1.75));
  vtkm::cont::ColorTable colorTable("inferno");

  colorTable.AddPointAlpha(0.0, .01f);
  colorTable.AddPointAlpha(1.0, .01f);

  // Background color:
  vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);
  vtkm::rendering::Scene scene;
  int width = 512, height = 512;
  CanvasRayTracer canvas(width, height);

  for (int i = 0; i < dsPerRank; i++)
  {
    //Create a sequence of datasets along the X direction.
    std::string fieldName = "tangle";
    vtkm::source::Tangle tangle;
    vtkm::Vec3f pt(rank * dsPerRank + i, 0, 0);
    if (rank == 1)
      std::cout << "PT= " << pt << std::endl;
    tangle.SetPointDimensions({ 50, 50, 50 });
    tangle.SetOrigin(pt);
    vtkm::cont::DataSet ds = tangle.Execute();

    if (unstructured)
    {
      vtkm::filter::clean_grid::CleanGrid cleanGrid;
      ds = cleanGrid.Execute(ds);
    }

    vtkm::rendering::Actor actor(
      ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldName), colorTable);
    scene.AddActor(actor);
  }



  if (unstructured)
  {
    vtkm::rendering::View3D view(scene, MapperConnectivity(), canvas, camera, bg);
    view.Paint();
    canvas.SaveAs("result-unstructured.png");
  }
  else
  {
    vtkm::rendering::View3D view(scene, MapperVolume(), canvas, camera, bg);
    view.Paint();
    canvas.SaveAs("result-structured.png");
  }
}

void RenderTests()
{
  //  TestImageComposite();
  //TestRenderComposite();
  TestVolumeRenderComposite(false);
  TestVolumeRenderComposite(true);
}

} //namespace

int UnitTestImageCompositing(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(RenderTests, argc, argv);
}
