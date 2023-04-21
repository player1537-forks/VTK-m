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

void PrintDescription(const std::string& testName, int blocksPerRank)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    std::cout << testName << " #ranks= " << comm.size() << " blocksPerRank= " << blocksPerRank
              << " ";
    std::cout << std::endl;
  }
}

void PrintDescription(const std::string& testName, int blocksPerRank, bool structured)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    std::cout << testName << " #ranks= " << comm.size() << " blocksPerRank= " << blocksPerRank
              << " ";
    if (structured)
      std::cout << "Structured";
    else
      std::cout << "Unstructured";

    std::cout << std::endl;
  }
}

vtkm::cont::PartitionedDataSet GenerateData(int baseSize,
                                            int blocksPerRank,
                                            bool doStructured = true)
{
  vtkm::cont::PartitionedDataSet result;

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  int rank = comm.rank();
  int totNumBlocks = blocksPerRank * comm.size();
  for (int i = 0; i < blocksPerRank; ++i)
  {
    int ID = rank * blocksPerRank + i;
    auto ds = CreateTestData(ID, totNumBlocks, baseSize);
    if (!doStructured)
    {
      vtkm::filter::clean_grid::CleanGrid cleanGrid;
      ds = cleanGrid.Execute(ds);
    }
    result.AppendPartition(ds);
  }

  return result;
}

void ScalarRenderer(int blocksPerRank)
{
  PrintDescription("ScalarRenderer", blocksPerRank);

  const int baseSize = 32;
  auto pds = GenerateData(baseSize, blocksPerRank);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.Azimuth(-30.f);
  camera.Elevation(-30.f);

  vtkh::ScalarRenderer tracer;

  std::cout << __FILE__ << " " << __LINE__ << " fix me" << std::endl;
  vtkm::rendering::Actor actor(pds, "", vtkm::cont::ColorTable("Cool to Warm"));

  tracer.SetInput(actor);
  tracer.SetCamera(camera);
  tracer.Update();

  vtkm::cont::PartitionedDataSet* output = tracer.GetOutput();
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    auto result = output->GetPartition(0);
    vtkm::io::VTKDataSetWriter writer("scalar_data.vtk");
    writer.WriteDataSet(result);
  }
}

void PointRenderer(bool renderVar, int blocksPerRank)
{
  PrintDescription("PointRenderer", blocksPerRank);

  const int baseSize = 16;
  auto pds = GenerateData(baseSize, blocksPerRank);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  std::string imgFile;
  if (renderVar)
    imgFile = "render_var_points";
  else
    imgFile = "render_points";

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(16, 36, -36));
  vtkh::Plot plot = vtkh::MakePlot(512, 512, camera, pds, imgFile);
  vtkh::PointRenderer renderer;

  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  renderer.SetInput(actor);
  if (renderVar)
  {
    renderer.SetUseVariableRadius(true);
    renderer.SetRadiusDelta(1.0f);
  }

  vtkh::Scene scene;
  scene.AddRenderer(&renderer);
  scene.AddPlot(plot);
  scene.Render();
}

void RayTrace(bool doStructured, int blocksPerRank)
{
  PrintDescription("RayTrace", blocksPerRank, doStructured);
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    std::cout << "RayTrace "
              << " blocksPerRank= " << blocksPerRank << " with: ";
    if (doStructured)
      std::cout << "structured data ";
    else
      std::cout << "unstructured data ";
    std::cout << "NumRanks= " << comm.size() << std::endl;
  }

  const int baseSize = 32;
  auto pds = GenerateData(baseSize, blocksPerRank, doStructured);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(-16, -16, -16));
  camera.ResetToBounds(bounds);

  std::string imgName = "ray_trace_par";
  if (doStructured)
    imgName += "_structured";
  else
    imgName += "_unstructured";

  vtkh::Plot plot = vtkh::MakePlot(512, 512, camera, pds, imgName);

  vtkh::RayTracer tracer;

  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  tracer.SetInput(actor);

  vtkh::Scene scene;
  scene.AddPlot(plot);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void VolumeRender(bool doStructured, int blocksPerRank)
{
  PrintDescription("VolumeRender", blocksPerRank, doStructured);

  const int baseSize = 32;
  auto pds = GenerateData(baseSize, blocksPerRank, doStructured);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(-16, -16, -16));
  camera.ResetToBounds(bounds);
  std::string imgName = "volume_par";
  if (doStructured)
    imgName += "_structured";
  else
    imgName += "_unstructured";

  vtkh::Plot plot = vtkh::MakePlot(512, 512, camera, pds, imgName);

  vtkm::cont::ColorTable color_map("cool to warm");
  color_map.AddPointAlpha(0.0, .05);
  color_map.AddPointAlpha(1.0, .5);

  //DRP: renderer needs to take an actor in ctor.
  vtkh::VolumeRenderer tracer;
  vtkm::rendering::Actor actor(pds, "point_data_Float64", color_map);
  tracer.SetInput(actor);
  tracer.SetColorTable(color_map);

  vtkh::Scene scene;
  scene.AddPlot(plot);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void VolumeRenderBlank(int blocksPerRank)
{
  PrintDescription("VolumeRenderBlank", blocksPerRank);

  const int baseSize = 32;
  auto pds = GenerateData(baseSize, blocksPerRank);

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
  vtkh::Plot plot = vtkh::MakePlot(512, 512, camera, isoOutput, "volume_unstructured_blank_par");


  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, 0.01);
  //color_map.AddPointAlpha(1.0, 0.2);
  color_map.AddPointAlpha(1.0, 0.6);

  vtkh::VolumeRenderer tracer;
  vtkm::rendering::Actor actor(pds, "point_data_Float64", color_map);
  tracer.SetInput(actor);
  tracer.SetColorTable(color_map);


  vtkh::Scene scene;
  scene.AddPlot(plot);
  scene.AddRenderer(&tracer);
  scene.Render();
}

void MultiRender(bool doBatch, int blocksPerRank)
{
  PrintDescription("MultiRender", blocksPerRank);

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  const int baseSize = 32;
  const int numBlocks = blocksPerRank * comm.size();
  auto pds = GenerateData(baseSize, blocksPerRank);

  vtkm::filter::contour::Contour contour;
  contour.SetActiveField("point_data_Float64");
  contour.SetFieldsToPass("cell_data_Float64");
  contour.SetIsoValue(0, -1); //ask for something that doesn't exist.
  contour.SetIsoValue(1, (float)baseSize * (float)numBlocks * 0.5f);

  auto iso_output = contour.Execute(pds);
  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  vtkh::Plot plot = vtkh::MakePlot(512, 512, camera, pds, "multi_par");
  vtkh::RayTracer tracer;
  vtkm::rendering::Actor actor(
    iso_output, "cell_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  tracer.SetInput(actor);

  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, .1);
  color_map.AddPointAlpha(1.0, .3);

  vtkm::rendering::Actor actor2(pds, "point_data_Float64", color_map);
  vtkh::VolumeRenderer v_tracer;
  v_tracer.SetInput(actor2);
  v_tracer.SetColorTable(color_map);


  vtkh::Scene scene;

  if (doBatch)
  {
    scene.SetRenderBatchSize(5);

    const int num_images = 11;
    std::vector<vtkh::Plot> renders;
    for (int i = 0; i < num_images; ++i)
    {
      vtkh::Plot tmp = plot.Copy();
      camera.Azimuth(float(i));
      tmp.SetCamera(camera);
      std::stringstream name;
      name << "par_batch_" << i;
      tmp.SetImageName(name.str());
      renders.push_back(tmp);
    }
    scene.SetPlots(renders);
  }
  else
  {
    scene.AddPlot(plot);
  }

  scene.AddRenderer(&v_tracer);
  scene.AddRenderer(&tracer);
  scene.Render();
}

//TODO
// serial pointrenderer fails.
// add serial versions
// Is VR_blank correct?


//Issues:
// Debug  :
//          Volume render hangs with nb=2
//
// Release: ScalarRender nb > 1
//          PointRenderer -- all
//          MultiRender nb > 1
//          VolumeRender hangs? with nb=2
//
void RenderTests()
{
  //std::vector<int> blocksPerRank = {1,2,3};
  std::vector<int> blocksPerRank = { 1 };
  std::vector<bool> flags = { true, false };

  for (auto nb : blocksPerRank)
  {
    //Actor wants a fieldName
    //ScalarRenderer(nb);

    for (auto v : flags)
    {
      //PointRenderer(v, nb);

      //Add PointRenderer no data
      MultiRender(v, nb);
      RayTrace(v, nb);
      VolumeRender(v, nb);
      VolumeRenderBlank(nb);
    }
  }
}

int UnitTestRendering_par(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(RenderTests, argc, argv);
}
