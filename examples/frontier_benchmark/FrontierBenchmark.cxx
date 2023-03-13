//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/Bounds.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/contour/Contour.h>
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/compositing/Image.h>
#include <vtkm/rendering/testing/RenderTest.h>
#include <vtkm/source/PerlinNoise.h>

struct BenchmarkOptions
{
  BenchmarkOptions(int argc, char** argv) { this->Parse(argc, argv); }

  void Parse(int argc, char** argv)
  {
    for (int i = 1; i < argc; ++i)
    {
      auto arg = std::string(argv[i]);
      auto seperatorPos = arg.find('=');
      if (seperatorPos == std::string::npos)
      {
        continue;
      }

      auto key = arg.substr(0, seperatorPos);
      auto val = arg.substr(seperatorPos + 1);

      if (key == "--perlin-dims")
      {
        int start = 0;
        for (int d = 0; d < 3; ++d)
        {
          auto end = val.find(',', start);
          auto dim = val.substr(start, end - start);
          this->PerlinDimensions[d] = std::stoi(dim);
          start = end + 1;
        }
      }
      else if (key == "--perlin-seed")
      {
        this->PerlinSeed = std::stoi(val);
      }
      else if (key == "--width")
      {
        this->CanvasWidth = std::stoi(val);
      }
      else if (key == "--height")
      {
        this->CanvasHeight = std::stoi(val);
      }
      else if (key == "--iters")
      {
        this->NumIterations = std::stoi(val);
      }
      else if (key == "--timing-file")
      {
        this->TimingFileName = val;
      }
    }
  }
  vtkm::Id3 PerlinDimensions = vtkm::Id3(1024, 1024, 1024);
  vtkm::IdComponent PerlinSeed = 1;
  vtkm::Id CanvasWidth = 1920;
  vtkm::Id CanvasHeight = 1080;
  vtkm::Id NumIterations = 10;
  std::string TimingFileName = "timing.csv";
};

struct MpiTopology
{
  MpiTopology(int rank, int size)
    : Rank(rank)
    , Size(size)
  {
  }

  int Rank;
  int Size;

  int XRank;
  int YRank;
  int ZRank;
  int XSize;
  int YSize;
  int ZSize;

  void SetShapeToCube()
  {
    int sizeCbrt = static_cast<int>(std::cbrt(static_cast<float>(this->Size)));
    this->ZSize = sizeCbrt;
    this->YSize = sizeCbrt;
    this->XSize = this->Size / (this->YSize * this->ZSize);

    this->XRank = this->Rank / (this->YSize * this->ZSize);
    this->YRank = (this->Rank - (this->XRank * this->YSize * this->ZSize)) / this->ZSize;
    this->ZRank =
      this->Rank - (this->XRank * this->YSize * this->ZSize) - (this->YRank * this->ZSize);
  }
};

struct IterationTimes
{
  vtkm::Float64 RenderTime = -1.0f;
  vtkm::Float64 CompositeTime = -1.0f;
  vtkm::Float64 TotalTime = -1.0f;
};

std::string GetImageName(const std::string& prefix,
                         const BenchmarkOptions& options,
                         const MpiTopology& mpiTopology)
{
  std::stringstream ss;
  ss << prefix << "_" << options.CanvasWidth << "x" << options.CanvasHeight << "_"
     << mpiTopology.Size << ".png";
  return ss.str();
}

void CollectTimeStats(vtkm::Float64 time,
                      vtkm::Float64& minTime,
                      vtkm::Float64& maxTime,
                      vtkm::Float64& avgTime)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  std::vector<vtkm::Float64> allTimes(comm.size(), -1.0f);
  if (comm.rank() == 0)
    vtkmdiy::mpi::gather(comm, time, allTimes, 0);
  else
    vtkmdiy::mpi::gather(comm, time, 0);
  if (comm.rank() == 0)
  {
    auto minMaxTime = std::minmax_element(allTimes.begin(), allTimes.end());
    minTime = *(minMaxTime.first);
    maxTime = *(minMaxTime.second);
    avgTime = std::accumulate(allTimes.begin(), allTimes.end(), 0.0) /
      static_cast<vtkm::Float64>(allTimes.size());
  }
}

void SaveTimeStats(const std::vector<IterationTimes>& stats,
                   const BenchmarkOptions& options,
                   const MpiTopology& vtkmNotUsed(mpiTopology))
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  if (comm.rank() != 0)
    return;

  std::vector<vtkm::Float64> renderTimes;
  std::vector<vtkm::Float64> compositeTimes;
  std::vector<vtkm::Float64> totalTimes;

  for (const auto& stat : stats)
  {
    renderTimes.push_back(stat.RenderTime);
    compositeTimes.push_back(stat.CompositeTime);
    totalTimes.push_back(stat.TotalTime);
  }

  auto CalculateTimeStats = [](const std::vector<vtkm::Float64>& times,
                               vtkm::Float64& minTime,
                               vtkm::Float64& maxTime,
                               vtkm::Float64& avgTime) {
    auto minMaxTime = std::minmax_element(times.begin(), times.end());
    minTime = *(minMaxTime.first);
    maxTime = *(minMaxTime.second);
    avgTime =
      std::accumulate(times.begin(), times.end(), 0.0) / static_cast<vtkm::Float64>(times.size());
  };

  vtkm::Float64 minRTime, maxRTime, avgRTime;
  vtkm::Float64 minCTime, maxCTime, avgCTime;
  vtkm::Float64 minTTime, maxTTime, avgTTime;
  CalculateTimeStats(renderTimes, minRTime, maxRTime, avgRTime);
  CalculateTimeStats(compositeTimes, minCTime, maxCTime, avgCTime);
  CalculateTimeStats(totalTimes, minTTime, maxTTime, avgTTime);

  auto ToHumanReadableMS = [](vtkm::Float64 time) { return static_cast<int>(time * 1000); };
  std::ofstream file;
  file.open(options.TimingFileName, std::ios::app);
  bool hasHeader = file.tellp() != 0;
  if (!hasHeader)
  {
    file << "World Size,Canvas Size,Min Render Time,Max Render Time,Avg Render Time,"
            "Min Composite Time,Max Composite Time,Avg Composite Time,"
            "Min Total Time,Max Total Time,Avg Total Time"
         << std::endl;
  }
  file << comm.size() << "," << options.CanvasWidth << "x" << options.CanvasHeight << ",";
  file << ToHumanReadableMS(minRTime) << "," << ToHumanReadableMS(maxRTime) << ","
       << ToHumanReadableMS(avgRTime) << ",";
  file << ToHumanReadableMS(minCTime) << "," << ToHumanReadableMS(maxCTime) << ","
       << ToHumanReadableMS(avgCTime) << ",";
  file << ToHumanReadableMS(minTTime) << "," << ToHumanReadableMS(maxTTime) << ","
       << ToHumanReadableMS(avgTTime) << std::endl;
  file.close();
}

void GenerateDataSet(const BenchmarkOptions& options,
                     const MpiTopology& mpiTopology,
                     vtkm::cont::DataSet& dataSet,
                     std::string& fieldName,
                     vtkm::Range& globalFiendRange,
                     vtkm::Bounds& globalBounds)
{
  fieldName = "perlinnoise";
  vtkm::source::PerlinNoise perlin;
  perlin.SetCellDimensions(options.PerlinDimensions);
  perlin.SetSeed(options.PerlinSeed);

  // Perlin Noise does not generate "interesting" surfaces when the scale is too small, .i.e,
  // X, Y, Z values are too close to each other. Hence we scale the perlin noise by a factor
  vtkm::FloatDefault scale = 3.0f;
  vtkm::Vec3f origin = { mpiTopology.XRank * scale,
                         mpiTopology.YRank * scale,
                         mpiTopology.ZRank * scale };
  vtkm::Vec3f maxExtent = origin + vtkm::Vec3f{ scale, scale, scale };
  perlin.SetOrigin(origin);
  perlin.SetMaxExtent(maxExtent);
  dataSet = perlin.Execute();
  vtkm::filter::contour::Contour contour;
  contour.SetNumberOfIsoValues(1);
  contour.SetIsoValue(0, 0.5f);
  contour.SetActiveField(fieldName);
  dataSet = contour.Execute(dataSet);

  globalFiendRange = { 0.0f, 1.0f };
  globalBounds = { vtkm::Vec3f{ 0.0f, 0.0f, 0.0f },
                   vtkm::Vec3f{ static_cast<vtkm::FloatDefault>(mpiTopology.XSize),
                                static_cast<vtkm::FloatDefault>(mpiTopology.YSize),
                                static_cast<vtkm::FloatDefault>(mpiTopology.ZSize) } *
                     scale };
  vtkm::Float64 boundsEps = 0.1f;
  globalBounds.Include(globalBounds.MinCorner() -
                       vtkm::Vec3f_64{ boundsEps, boundsEps, boundsEps });
  globalBounds.Include(globalBounds.MaxCorner() +
                       vtkm::Vec3f_64{ boundsEps, boundsEps, boundsEps });
}

void RunBenchmark(const BenchmarkOptions& options)
{
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MpiTopology mpiTopology(comm.rank(), comm.size());
  mpiTopology.SetShapeToCube();

  // Generate DataSet
  vtkm::cont::DataSet dataSet;
  std::string fieldName;
  vtkm::Range globalFiendRange;
  vtkm::Bounds globalBounds;
  GenerateDataSet(options, mpiTopology, dataSet, fieldName, globalFiendRange, globalBounds);

  vtkm::rendering::Scene scene;
  vtkm::cont::ColorTable colorTable("inferno");
  vtkm::rendering::Actor actor(
    dataSet.GetCellSet(), dataSet.GetCoordinateSystem(), dataSet.GetField(fieldName), colorTable);
  actor.SetScalarRange(globalFiendRange);
  scene.AddActor(actor);

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(globalBounds);
  camera.Azimuth(10.0f);
  camera.Elevation(20.0f);

  std::vector<IterationTimes> benchmarkTimes;
  vtkm::rendering::CanvasRayTracer canvas(options.CanvasWidth, options.CanvasHeight);
  for (int iter = 0; iter < options.NumIterations; iter++)
  {
    vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);
    vtkm::rendering::View3D view(scene, vtkm::rendering::MapperRayTracer(), canvas, camera, bg);
    view.Paint();
    if (comm.rank() == 0)
    {
      benchmarkTimes.push_back(
        IterationTimes{ .RenderTime = view.GetTimes()[vtkm::rendering::RENDER_TIME_KEY],
                        .CompositeTime = view.GetTimes()[vtkm::rendering::COMPOSITE_TIME_KEY],
                        .TotalTime = view.GetTimes()[vtkm::rendering::TOTAL_TIME_KEY] });
    }
  }

  SaveTimeStats(benchmarkTimes, options, mpiTopology);

  canvas.SaveAs(GetImageName("perlin_surface", options, mpiTopology));
}

int main(int argc, char* argv[])
{
  // Initialize MPI
  std::unique_ptr<vtkmdiy::mpi::environment> mpiEnv = nullptr;
  if (!vtkmdiy::mpi::environment::initialized())
  {
    mpiEnv.reset(new vtkmdiy::mpi::environment(argc, argv));
  }

  // Initialize VTK-m
  vtkm::cont::Initialize(argc, argv, vtkm::cont::InitializeOptions::None);

  BenchmarkOptions options(argc, argv);
  RunBenchmark(options);
  return 0;
}