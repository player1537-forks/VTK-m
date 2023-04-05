//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

//#include <vtkm/cont/testing/Testing.h>
#include <vtkm/filter/contour/ClipWithField.h>
#include <vtkm/filter/contour/Contour.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/rendering/testing/RenderTest.h>
#include <vtkm/rendering/vtkh/rendering/PointRenderer.hpp>
#include <vtkm/rendering/vtkh/rendering/RayTracer.hpp>
#include <vtkm/rendering/vtkh/rendering/ScalarRenderer.hpp>
#include <vtkm/rendering/vtkh/rendering/Scene.hpp>
#include <vtkm/rendering/vtkh/rendering/VolumeRenderer.hpp>

#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/rendering/testing/t_vtkm_test_utils.hpp>

void ScalarRenderer()
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
    std::cout << "ScalarRenderer" << std::endl;

  int rank = comm.rank();
  const int base_size = 32;
  const int blocks_per_rank = 2;
  const int num_blocks = comm.size() * blocks_per_rank;

  vtkm::cont::PartitionedDataSet pds;
  for (int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    auto ds = CreateTestData(domain_id, num_blocks, base_size);
    pds.AppendPartition(ds);
  }

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.Azimuth(-30.f);
  camera.Elevation(-30.f);

  vtkh::ScalarRenderer tracer;

  tracer.SetInput(&pds);
  tracer.SetCamera(camera);
  tracer.Update();

  vtkm::cont::PartitionedDataSet* output = tracer.GetOutput();

  if (comm.rank() == 0)
  {
    auto result = output->GetPartition(0);
    vtkm::io::VTKDataSetWriter writer("scalar_data.vtk");
    writer.WriteDataSet(result);
  }
}

void PointRenderer(bool renderVar)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
    std::cout << "PointRenderer" << std::endl;

  const int base_size = 16;
  const int num_blocks = 2;

  vtkm::cont::PartitionedDataSet pds;
  for (int i = 0; i < num_blocks; ++i)
  {
    auto ds = CreateTestData(i, num_blocks, base_size);
    pds.AppendPartition(ds);
  }

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  std::string imgFile;
  if (renderVar)
    imgFile = "render_var_points";
  else
    imgFile = "render_points";

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(16, 36, -36));
  vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, imgFile);
  vtkh::PointRenderer renderer;
  renderer.SetInput(&pds);
  renderer.SetField("point_data_Float64");
  if (renderVar)
  {
    renderer.SetUseVariableRadius(true);
    renderer.SetRadiusDelta(1.0f);
  }

  vtkh::Scene scene;
  scene.AddRenderer(&renderer);
  scene.AddRender(render);
  scene.Render();
}

void RayTrace(bool doStructured)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    std::cout << "RayTrace with: ";
    if (doStructured)
      std::cout << "structured data ";
    else
      std::cout << "unstructured data ";
    std::cout << "NumRanks= " << comm.size() << std::endl;
  }

  int rank = comm.rank();
  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int num_blocks = comm.size() * blocks_per_rank;

  vtkm::cont::PartitionedDataSet pds;
  //vtkh::DataSet data_set;
  for (int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    auto ds = CreateTestData(domain_id, num_blocks, base_size);
    if (!doStructured)
    {
      vtkm::filter::clean_grid::CleanGrid cleanGrid;
      ds = cleanGrid.Execute(ds);
    }
    //data_set.AddDomain(ds, domain_id);
    pds.AppendPartition(ds);
  }

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(-16, -16, -16));
  camera.ResetToBounds(bounds);

  std::string imgName = "ray_trace_par";
  if (doStructured)
    imgName += "_structured";
  else
    imgName += "_unstructured";

  vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, imgName);

  vtkh::RayTracer tracer;

  tracer.SetInput(&pds);
  tracer.SetField("point_data_Float64");

  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void VolumeRender(bool doStructured)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  if (comm.rank() == 0)
  {
    std::cout << "VolumeRender with: ";
    if (doStructured)
      std::cout << "structured data ";
    else
      std::cout << "unstructured data ";
    std::cout << "NumRanks= " << comm.size() << std::endl;
  }

  int rank = comm.rank();

  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int num_blocks = comm.size() * blocks_per_rank;

  //vtkh::DataSet data_set;
  vtkm::cont::PartitionedDataSet pds;
  for (int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    auto ds = CreateTestData(domain_id, num_blocks, base_size);
    if (!doStructured)
    {
      vtkm::filter::clean_grid::CleanGrid cleanGrid;
      ds = cleanGrid.Execute(ds);
    }
    pds.AppendPartition(ds);
    //data_set.AddDomain(ds, domain_id);
  }

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(-16, -16, -16));
  camera.ResetToBounds(bounds);
  std::string imgName = "volume_par";
  if (doStructured)
    imgName += "_structured";
  else
    imgName += "_unstructured";

  vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, imgName);

  vtkm::cont::ColorTable color_map("cool to warm");
  color_map.AddPointAlpha(0.0, .05);
  color_map.AddPointAlpha(1.0, .5);

  vtkh::VolumeRenderer tracer;
  tracer.SetColorTable(color_map);
  tracer.SetInput(&pds);
  tracer.SetField("point_data_Float64");

  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void VolumeRenderBlank()
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  if (comm.rank() == 0)
  {
    std::cout << "VolumeRender blank ";
    std::cout << "NumRanks= " << comm.size() << std::endl;
  }

  int rank = comm.rank();

  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int num_blocks = comm.size() * blocks_per_rank;

  vtkm::cont::PartitionedDataSet pds;
  for (int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    auto ds = CreateTestData(domain_id, num_blocks, base_size);
    pds.AppendPartition(ds);
  }



  vtkm::filter::contour::ClipWithField min_clipper, max_clipper;
  max_clipper.SetActiveField("point_data_Float64");
  max_clipper.SetInvertClip(true);
  max_clipper.SetClipValue(40000.);
  auto maxResult = max_clipper.Execute(pds);

  min_clipper.SetActiveField("point_data_Float64");
  min_clipper.SetInvertClip(false);
  max_clipper.SetClipValue(10000.);
  auto isoOutput = min_clipper.Execute(maxResult);

  vtkm::rendering::Camera camera;
  vtkm::Vec<vtkm::Float32, 3> pos = camera.GetPosition();
  pos[0] += 10000.1;
  pos[1] += 10000.1;
  camera.SetPosition(pos);
  vtkm::Vec<vtkm::Float32, 3> look;
  look[0] = 100000.f;
  look[1] = 100000.f;
  look[2] = 100000.f;
  camera.SetLookAt(look);
  vtkh::Render render =
    vtkh::MakeRender(512, 512, camera, isoOutput, "volume_unstructured_blank_par");


  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, 0.01);
  //color_map.AddPointAlpha(1.0, 0.2);
  color_map.AddPointAlpha(1.0, 0.6);

  vtkh::VolumeRenderer tracer;
  tracer.SetColorTable(color_map);
  tracer.SetInput(&isoOutput);
  tracer.SetField("point_data_Float64");

  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void MultiRender(bool doBatch)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  if (comm.rank() == 0)
  {
    std::cout << "MultiRender ";
    if (doBatch)
      std::cout << " Batch";
    std::cout << std::endl;
  }

  int rank = comm.rank();

  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int num_blocks = comm.size() * blocks_per_rank;

  vtkm::cont::PartitionedDataSet pds;
  for (int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    auto ds = CreateTestData(domain_id, num_blocks, base_size);
    pds.AppendPartition(ds);
  }

  vtkm::filter::contour::Contour contour;
  contour.SetActiveField("point_data_Float64");
  contour.SetFieldsToPass("cell_data_Float64");
  contour.SetIsoValue(0, -1); //ask for something that doesn't exist.
  contour.SetIsoValue(1, (float)base_size * (float)num_blocks * 0.5f);

  auto iso_output = contour.Execute(pds);
  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  vtkh::Render render = vtkh::MakeRender(512, 512, camera, pds, "multi_par");
  vtkh::RayTracer tracer;
  tracer.SetInput(&iso_output);
  tracer.SetField("cell_data_Float64");

  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, .1);
  color_map.AddPointAlpha(1.0, .3);

  vtkh::VolumeRenderer v_tracer;
  v_tracer.SetColorTable(color_map);
  v_tracer.SetInput(&pds);
  v_tracer.SetField("point_data_Float64");

  vtkh::Scene scene;

  if (doBatch)
  {
    scene.SetRenderBatchSize(5);

    const int num_images = 11;
    std::vector<vtkh::Render> renders;
    for (int i = 0; i < num_images; ++i)
    {
      vtkh::Render tmp = render.Copy();
      camera.Azimuth(float(i));
      tmp.SetCamera(camera);
      std::stringstream name;
      name << "par_batch_" << i;
      tmp.SetImageName(name.str());
      renders.push_back(tmp);
    }
    scene.SetRenders(renders);
  }
  else
  {
    scene.AddRender(render);
  }

  scene.AddRenderer(&v_tracer);
  scene.AddRenderer(&tracer);
  scene.Render();
}

//TODO
// add serial versions
// Is VR_blank correct?

void RenderTests()
{
  ScalarRenderer();
  PointRenderer(true);
  PointRenderer(false);

  //Add PointRenderer no data
  MultiRender(true);
  MultiRender(false);

  RayTrace(true);
  RayTrace(false);

  VolumeRender(true);
  VolumeRender(false);

  VolumeRenderBlank();
}

int UnitTestRendering_par(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(RenderTests, argc, argv);
}
