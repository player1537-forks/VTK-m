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
#include <vtkm/rendering_new/PointRenderer.h>
#include <vtkm/rendering_new/RayTracer.h>
#include <vtkm/rendering_new/ScalarRenderer.h>
#include <vtkm/rendering_new/Scene.h>
#include <vtkm/rendering_new/VolumeRenderer.h>

#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/rendering/testing/t_vtkm_test_utils.hpp>

void PrintDescription(const std::string& testName, int blocksPerRank)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    std::cout << std::endl;
    std::cout << "*******************************************************************" << std::endl;
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
    std::cout << std::endl;
    std::cout << "*******************************************************************" << std::endl;
    std::cout << testName << " #ranks= " << comm.size() << " blocksPerRank= " << blocksPerRank
              << " ";
    if (structured)
      std::cout << "Structured";
    else
      std::cout << "Unstructured";

    std::cout << std::endl;
  }
}

std::string GetImageName(const std::string& testName, bool structured, int blocksPerRank)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  std::string imgName = testName + "_rank" + std::to_string(comm.size());
  if (structured)
    imgName += "_struct";
  else
    imgName += "_unstruct";
  imgName += "_bpr" + std::to_string(blocksPerRank);

  return imgName;
}

std::string GetImageName(const std::string& testName, int blocksPerRank)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();

  std::string imgName = testName + "_rank" + std::to_string(comm.size());
  imgName += "_bpr" + std::to_string(blocksPerRank);

  return imgName;
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
  pds.PrintSummary(std::cout);

  vtkm::Bounds bounds = pds.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.Azimuth(-30.f);
  camera.Elevation(-30.f);
  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, "TMP");

  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));

  vtkm::rendering_new::ScalarRenderer tracer;
  tracer.SetInput(actor);
  plot.AddRenderer(&tracer);
  vtkm::rendering_new::Scene scene;
  scene.AddPlot(plot);
  scene.Render();

  vtkm::cont::PartitionedDataSet* output = tracer.GetOutput();
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    auto result = output->GetPartition(0);
    vtkm::io::VTKDataSetWriter writer("scalar_data.vtk");
    writer.WriteDataSet(result);
  }

  /*
  std::cout<<__FILE__<<" "<<__LINE__<<" fix me"<<std::endl;

  vtkm::rendering::Actor actor(pds, "", vtkm::cont::ColorTable("Cool to Warm"));

  tracer.SetInput(actor);
  tracer.SetCamera(camera);
  std::vector<vtkm::rendering_new::Plot> plots;
  std::cout<<"FIX ME!!!!!!!!!!"<<std::endl;
  tracer.Update(plots);

  vtkm::cont::PartitionedDataSet* output = tracer.GetOutput();
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() == 0)
  {
    auto result = output->GetPartition(0);
    vtkm::io::VTKDataSetWriter writer("scalar_data.vtk");
    writer.WriteDataSet(result);
  }
  */
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
  imgFile = GetImageName(imgFile, blocksPerRank);

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  camera.SetPosition(vtkm::Vec<vtkm::Float64, 3>(16, 36, -36));
  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, imgFile);
  vtkm::rendering_new::PointRenderer renderer;

  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  renderer.SetInput(actor);
  if (renderVar)
  {
    renderer.SetUseVariableRadius(true);
    renderer.SetRadiusDelta(1.0f);
  }

  vtkm::rendering_new::Scene scene;
  plot.AddRenderer(&renderer);
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

  std::string imgName = GetImageName("ray_trace_par", doStructured, blocksPerRank);

  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, imgName);

  vtkm::rendering_new::RayTracer tracer;

  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  tracer.SetInput(actor);
  plot.AddRenderer(&tracer);

  vtkm::rendering_new::Scene scene;
  scene.AddPlot(plot);
  scene.Render();
}

void RayTrace2(bool doStructured, int blocksPerRank)
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

  std::string imgName = GetImageName("RayTrace2", doStructured, blocksPerRank);
  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, imgName);


  vtkm::rendering_new::RayTracer rendererRT;
  vtkm::rendering::Actor actor(pds, "point_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  rendererRT.SetInput(actor);

  plot.AddRenderer(&rendererRT);

  vtkm::rendering_new::Scene scene;
  scene.AddPlot(plot);
  //scene.AddRenderer(&rendererRT);
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
  std::string imgName = GetImageName("volume_par", doStructured, blocksPerRank);

  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, imgName);

  vtkm::cont::ColorTable color_map("cool to warm");
  color_map.AddPointAlpha(0.0, .05);
  color_map.AddPointAlpha(1.0, .5);

  //DRP: renderer needs to take an actor in ctor.
  vtkm::rendering_new::VolumeRenderer tracer;
  vtkm::rendering::Actor actor(pds, "point_data_Float64", color_map);
  tracer.SetInput(actor);
  tracer.SetColorTable(color_map);

  plot.AddRenderer(&tracer);
  vtkm::rendering_new::Scene scene;
  scene.AddPlot(plot);
  //scene.AddRenderer(&tracer);
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
  std::string imgName = GetImageName("volume_blank", false, blocksPerRank);
  vtkm::rendering_new::Plot plot =
    vtkm::rendering_new::MakePlot(512, 512, camera, isoOutput, imgName);

  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, 0.01);
  //color_map.AddPointAlpha(1.0, 0.2);
  color_map.AddPointAlpha(1.0, 0.6);

  vtkm::rendering_new::VolumeRenderer tracer;
  vtkm::rendering::Actor actor(pds, "point_data_Float64", color_map);
  tracer.SetInput(actor);
  tracer.SetColorTable(color_map);
  plot.AddRenderer(&tracer);

  vtkm::rendering_new::Scene scene;
  scene.AddPlot(plot);
  //scene.AddRenderer(&tracer);
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

  std::string imgName = GetImageName("multiRender", blocksPerRank);
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  vtkm::rendering_new::Plot plot = vtkm::rendering_new::MakePlot(512, 512, camera, pds, imgName);
  vtkm::rendering_new::RayTracer tracer;
  vtkm::rendering::Actor actor(
    iso_output, "cell_data_Float64", vtkm::cont::ColorTable("Cool to Warm"));
  tracer.SetInput(actor);

  vtkm::cont::ColorTable color_map("Cool to Warm");
  color_map.AddPointAlpha(0.0, .1);
  color_map.AddPointAlpha(1.0, .3);

  vtkm::rendering::Actor actor2(pds, "point_data_Float64", color_map);
  vtkm::rendering_new::VolumeRenderer v_tracer;
  v_tracer.SetInput(actor2);
  v_tracer.SetColorTable(color_map);


  plot.AddRenderer(&v_tracer);
  plot.AddRenderer(&tracer);

  vtkm::rendering_new::Scene scene;
  if (doBatch)
  {
    scene.SetRenderBatchSize(5);

    const int num_images = 11;
    std::vector<vtkm::rendering_new::Plot> plots;
    for (int i = 0; i < num_images; ++i)
    {
      vtkm::rendering_new::Plot tmp = plot.Copy();
      camera.Azimuth(float(i));
      tmp.SetCamera(camera);
      std::stringstream name;
      name << imgName << "_"
           << "par_batch_" << i;
      tmp.SetImageName(name.str());
      plots.push_back(tmp);
    }
    scene.SetPlots(plots);
  }
  else
  {
    scene.AddPlot(plot);
  }

  //  scene.AddRenderer(&v_tracer);
  //  scene.AddRenderer(&tracer);
  //  scene.RenderORIG();
  scene.Render();
}

//TODO
// 1 rank issues with DIY. Vicente working on this
// ScalarRenderer -- not clear what this is supposed to do.
// PointRenderer: setting bbox is off with 4 ranks.



// serial pointrenderer fails.
// add serial versions
// Is VR_blank correct?
//PointRenderer: needs ParticleMerging filter.


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
  RayTrace2(true, 1);
  PointRenderer(false, 1);
  return;

  //ScalarRenderer(1);
  //return;


  /*
  MultiRender(true, 1);
  RayTrace2(true, 1);
  VolumeRender(true, 1);
  return;
  //RayTrace2(true, 2);
  */

#if 1

  std::vector<int> blocksPerRank = { 1, 2 };
  //std::vector<int> blocksPerRank = {1};
  std::vector<bool> flags = { true, false };

  for (auto nb : blocksPerRank)
  {
    //Actor wants a fieldName
    //ScalarRenderer(nb);

    for (auto v : flags)
    {
      //PointRenderer(v, nb); //doesn't work...

      //Add PointRenderer no data
      MultiRender(v, nb);
      RayTrace2(v, nb);
      VolumeRender(v, nb);
      VolumeRenderBlank(nb);
    }
  }
#endif
}

int UnitTestRendering_par(int argc, char* argv[])
{

  return vtkm::cont::testing::Testing::Run(RenderTests, argc, argv);
}
