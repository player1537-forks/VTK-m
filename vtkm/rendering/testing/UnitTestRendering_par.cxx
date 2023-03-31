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
#include <vtkm/rendering/testing/RenderTest.h>
#include <vtkm/rendering/vtkh/rendering/RayTracer.hpp>
#include <vtkm/rendering/vtkh/rendering/Scene.hpp>
#include <vtkm/rendering/vtkh/rendering/VolumeRenderer.hpp>

#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/rendering/testing/t_vtkm_test_utils.hpp>

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


void RenderTests()
{
  RayTrace(true);
  RayTrace(false);

  VolumeRender(true);
  VolumeRender(false);
}

int UnitTestRendering_par(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(RenderTests, argc, argv);
}
