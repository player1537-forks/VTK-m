//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "Benchmarker.h"

#include <vtkm/TypeTraits.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/contour/Contour.h>

#include <vtkm/source/PerlinNoise.h>
#include <vtkm/source/Tangle.h>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/raytracing/RayTracer.h>
#include <vtkm/rendering/raytracing/SphereIntersector.h>
#include <vtkm/rendering/raytracing/TriangleExtractor.h>

#include <vtkm/rendering/vtkh/rendering/Plot.h>
#include <vtkm/rendering/vtkh/rendering/RayTracer.h>
#include <vtkm/rendering/vtkh/rendering/Scene.h>

#include <vtkm/exec/FunctorBase.h>

#include <sstream>
#include <string>
#include <vector>

namespace
{

enum CameraModeType
{
  Static = 0,
  Orbit = 1,
};

enum PerlinModeType
{
  Grow = 0,
  Subdivide = 1,
};

struct BenchmarkOptionsType
{
  vtkm::Id3 PerlinDimensions = vtkm::Id3(128, 128, 128);
  vtkm::IdComponent PerlinSeed = 1;
  vtkm::Float32 PerlinScale = 3.0f;
  vtkm::Id CanvasWidth = 128;
  vtkm::Id CanvasHeight = 128;
  vtkm::Id NumIterations = 2;
  std::string TimingFileName = "timing.csv";
  std::string ImageFormat = "png";
  CameraModeType CameraMode = CameraModeType::Static;
  PerlinModeType PerlinMode = PerlinModeType::Grow;
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

vtkm::cont::PartitionedDataSet GenerateDataSet(const BenchmarkOptionsType& options,
                                               const MpiTopology& mpiTopology,
                                               std::string& fieldName,
                                               vtkm::Range& globalFieldRange,
                                               vtkm::Bounds& globalBounds)
{
  fieldName = "perlinnoise";
  vtkm::source::PerlinNoise perlin;
  perlin.SetCellDimensions(options.PerlinDimensions);
  perlin.SetSeed(options.PerlinSeed);

  vtkm::cont::DataSet ds;

  // Perlin Noise does not generate "interesting" surfaces when the scale is too small, .i.e,
  // X, Y, Z values are too close to each other. Hence we scale the perlin noise by a factor
  if (options.PerlinMode == PerlinModeType::Grow)
  {
    vtkm::Vec3f origin = vtkm::Vec3f_32{ static_cast<vtkm::Float32>(mpiTopology.XRank),
                                         static_cast<vtkm::Float32>(mpiTopology.YRank),
                                         static_cast<vtkm::Float32>(mpiTopology.ZRank) } *
      options.PerlinScale;
    vtkm::Vec3f maxExtent = origin + vtkm::Vec3f{ options.PerlinScale };
    perlin.SetOrigin(origin);
    perlin.SetMaxExtent(maxExtent);
    ds = perlin.Execute();
  }
  else if (options.PerlinMode == PerlinModeType::Subdivide)
  {
    vtkm::Vec3f_32 blockSpacing = vtkm::Vec3f_32{ options.PerlinScale } *
      vtkm::Vec3f_32{ 1.0f / static_cast<vtkm::Float32>(mpiTopology.XSize),
                      1.0f / static_cast<vtkm::Float32>(mpiTopology.YSize),
                      1.0f / static_cast<vtkm::Float32>(mpiTopology.ZSize) };
    vtkm::Vec3f_32 origin = vtkm::Vec3f_32{ static_cast<vtkm::Float32>(mpiTopology.XRank),
                                            static_cast<vtkm::Float32>(mpiTopology.YRank),
                                            static_cast<vtkm::Float32>(mpiTopology.ZRank) } *
      blockSpacing;

    vtkm::Vec3f_32 maxExtent = origin + blockSpacing;
    perlin.SetOrigin(origin);
    perlin.SetMaxExtent(maxExtent);
    ds = perlin.Execute();
  }


  std::vector<vtkm::Float64> isoValues{ 0.4f, 0.75f };
  vtkm::filter::contour::Contour contour;
  contour.SetIsoValues(isoValues);
  contour.SetActiveField(fieldName);
  ds = contour.Execute(ds);

  auto dataSet = vtkm::cont::PartitionedDataSet(ds);

  globalBounds = dataSet.GetGlobalBounds();
  globalFieldRange = { 0.0f, 1.0f };

  // Add a small epsilon to the bounds to prevent the world annotations from being clipped
  if (options.CameraMode == CameraModeType::Orbit)
  {
    vtkm::Float64 boundsEps = 0.2f;
    globalBounds.Include(globalBounds.MinCorner() -
                         vtkm::Vec3f_64{ boundsEps, boundsEps, boundsEps });
    globalBounds.Include(globalBounds.MaxCorner() +
                         vtkm::Vec3f_64{ boundsEps, boundsEps, boundsEps });
  }

  return dataSet;
}



// Hold configuration state (e.g. active device)
vtkm::cont::InitializeResult Config;

void BenchRayTracingMPI(::benchmark::State& state)
{
  BenchmarkOptionsType options;

  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MpiTopology mpiTopology(comm.rank(), comm.size());
  mpiTopology.SetShapeToCube();

  //Generate the data
  std::string fieldName;
  vtkm::Range globalFieldRange;
  vtkm::Bounds globalBounds;
  auto dataSet = GenerateDataSet(options, mpiTopology, fieldName, globalFieldRange, globalBounds);

  vtkm::rendering::Camera camera;
  camera.ResetToBounds(globalBounds);

  vtkm::rendering::Color color(1, 1, 1);
  vtkm::rendering::Actor actor(dataSet, fieldName, color);
  vtkh::RayTracer rayTracer;
  rayTracer.SetInput(actor);

  vtkh::Plot plot =
    vtkh::MakePlot(options.CanvasWidth, options.CanvasHeight, camera, dataSet, "tmp");

  vtkm::cont::Timer timer{ Config.Device };
  for (auto _ : state)
  {
    (void)_;
    timer.Start();

    vtkh::Scene scene;
    plot.AddRenderer(&rayTracer);
    scene.AddPlot(plot);
    scene.Render();
    timer.Stop();

    if (comm.rank() == 0)
      state.SetIterationTime(timer.GetElapsedTime());
  }
}

VTKM_BENCHMARK(BenchRayTracingMPI);

} // end namespace vtkm::benchmarking

int main(int argc, char* argv[])
{
  auto opts = vtkm::cont::InitializeOptions::RequireDevice;

  std::vector<char*> args(argv, argv + argc);
  vtkm::bench::detail::InitializeArgs(&argc, args, opts);

  // Parse VTK-m options:
  Config = vtkm::cont::Initialize(argc, args.data(), opts);

  // This occurs when it is help
  if (opts == vtkm::cont::InitializeOptions::None)
  {
    std::cout << Config.Usage << std::endl;
  }
  else
  {
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(Config.Device);
  }

  // handle benchmarking related args and run benchmarks:
  VTKM_EXECUTE_BENCHMARKS(argc, args.data());
}
