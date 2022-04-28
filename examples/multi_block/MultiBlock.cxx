//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/filter/contour/Contour.h>
#include <vtkm/source/Tangle.h>
#include <vtkm/cont/cuda/internal/CudaAllocator.h>

#include <cstdlib>
#include <iostream>

int rank = 0;
int numRanks = 1;


struct Options
{
public:
  enum RunModeType {SERIAL, OPENMP, GPU};

  std::string DataPath = "";
  std::string DataFile = "";
  std::string Field = "";
  std::string MapField = "";
  int IsoLevels = 1;
  std::vector<double> IsoValues;
  int Rank = 0;
  int NumRanks = 1;
  std::string ThreadMode = "serial";
  int NumTasks = 1;
  bool Tangle = false;
  vtkm::Id NumTangle;
  vtkm::Id3 TangleDims;
  RunModeType RunMode = SERIAL;

  bool ParseOptions(int argc, char **argv)
  {
    if (argc < 3)
    {
      this->PrintUsage();
      return false;
    }

    //Put args into a list of pairs.
    std::vector<std::pair<std::string, std::vector<std::string>>> args;
    for (int i = 1; i < argc; i++)
    {
      std::string tmp = argv[i];
      if (tmp.find("--") != std::string::npos)
      {
        std::pair<std::string, std::vector<std::string>> entry;
        entry.first = tmp;
        args.push_back(entry);
      }
      else
        args.back().second.push_back(tmp);
    }

    //Fill out member data from args.
    this->DataPath = ".";
    for (const auto& a : args)
    {
      if (a.first == "--data")
      {
        if (a.second.empty()) return false;

        auto tmp = a.second[0];
        auto pos = tmp.rfind("/");
        if (pos != std::string::npos)
        {
          this->DataPath = tmp.substr(0, pos);
          this->DataFile = tmp.substr(pos+1, tmp.size()-pos-1);
        }
      }
      else if (a.first == "--tangle")
      {
        this->Tangle = true;
        this->NumTangle = std::atoi(a.second[0].c_str());
        this->TangleDims[0] = std::atoi(a.second[1].c_str());
        this->TangleDims[1] = std::atoi(a.second[2].c_str());
        this->TangleDims[2] = std::atoi(a.second[3].c_str());
      }
      else if (a.first == "--field")
      {
        if (a.second.empty()) return false;
        this->Field = a.second[0];
      }
      else if (a.first == "--mapfield")
      {
        if (a.second.empty()) return false;
        this->MapField = a.second[0];
      }
      else if (a.first == "--threading")
      {
        if (a.second.empty()) return false;
        if (a.second[0] == "serial")
          this->ThreadMode = "serial";
        else if (a.second[0] == "openmp")
          this->ThreadMode = "openmp";
        else if (a.second[0] == "task")
        {
          if (a.second.size() != 2)
            return false;
          this->NumTasks = std::stoi(a.second[1]);
          this->ThreadMode = "task";
        }
      }
      else if (a.first == "--isolevels")
      {
        if (a.second.empty()) return false;
        this->IsoLevels = std::stoi(a.second[0]);
      }
      else if (a.first == "--isovalues")
      {
        if (a.second.empty()) return false;
        this->IsoValues.clear();
        for (const auto& aa : a.second)
          this->IsoValues.push_back(std::stod(aa));
      }
      else if (a.first == "--serial")
        this->RunMode = SERIAL;
      else if (a.first == "--openmp")
        this->RunMode = OPENMP;
      else if (a.first == "--gpu")
        this->RunMode = GPU;
    }

    if (this->MapField == "")
      this->MapField = this->Field;

    if ((!this->Tangle && (this->DataFile == "" || this->ThreadMode == "")) ||
        this->Field == "")
    {
      if (this->Rank == 0) std::cerr<<"Error in options"<<std::endl;
      return false;
    }

    return true;
  }

  void PrintUsage() const
  {
    std::cerr<<"Usage: --data <dataFile> or --tangle n d0 d1 d2 --field <field> --mapfield <mapfield> --thread <serial openmp task N>"<<std::endl;
  }
};



vtkm::cont::PartitionedDataSet
ReadVisItFile(const std::string& dir, const std::string& fname)
{
  std::string tmp = dir + "/" + fname;
  std::cout<<"Opening: "<<tmp<<std::endl;
  std::ifstream visitFile;
  visitFile.open(dir + "/" + fname);

  vtkm::cont::PartitionedDataSet dataSets;

  std::string f;
  std::vector<std::string> fileNames;
  int cnt = 0;
  while (visitFile >> f)
  {
    if (cnt > 1) //skip !NBLOCKS N
      fileNames.push_back(f);
    cnt++;
  }

  int numDS = fileNames.size();
  int numPerRank = numDS / numRanks;
  int b0 = rank*numPerRank;
  int b1 = b0 + numPerRank;
  int leftOver = numDS % numRanks;

  std::cout<<"numDS= "<<numDS<<std::endl;

  std::vector<int> domainIds;
  std::vector<std::string> myFiles;
  for (int i = b0; i < b1; i++)
  {
    myFiles.push_back(fileNames[i]);
    domainIds.push_back(i);
  }
  if (leftOver > 0 && rank < leftOver)
  {
    myFiles.push_back(fileNames[numPerRank*numRanks + rank]);
    domainIds.push_back(numPerRank*numRanks + rank);
  }

  for (std::size_t i = 0; i < myFiles.size(); i++)
  {
    //std::string f = dir + "/" + myFiles[i];
    //std::cout<<"Reading: "<<f<<std::endl;
    vtkm::io::VTKDataSetReader reader(dir + "/" + myFiles[i]);
    auto ds = reader.ReadDataSet();
    dataSets.AppendPartition(ds);
  }
  return dataSets;
}

int main(int argc, char** argv)
{
  vtkm::cont::Initialize(argc, argv);

  /*
  vtkm::source::Tangle tangle({20,20,20});
  auto ds = tangle.Execute();
  vtkm::io::VTKDataSetWriter writer("out.vtk");
  writer.WriteDataSet(ds);
  return 0;
  */

  Options opts;
  if (!opts.ParseOptions(argc, argv))
  {
    opts.PrintUsage();
    return 0;
  }

  vtkm::Id numCPUThreads = 1, numGPUThreads = 1;  
  if (opts.ThreadMode == "task")
  {
    numCPUThreads = opts.NumTasks;
    numGPUThreads = opts.NumTasks;
  }
  if (opts.RunMode == Options::RunModeType::SERIAL)
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});
  else if (opts.RunMode == Options::RunModeType::OPENMP)
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP{});
  else if (opts.RunMode == Options::RunModeType::GPU)
  {
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});
#ifdef VTKM_CUDA
    if (opts.ThreadMode == "task")
    {
        vtkm::cont::cuda::internal::CudaAllocator::UsingManagedMemory();
        std::cout<<"MULTIBLOCK:                           Turn ManagedMemoryOff()"<<std::endl;
        vtkm::cont::cuda::internal::CudaAllocator::ForceManagedMemoryOff();
    }
#endif
  }

  vtkm::cont::PartitionedDataSet dataSets;
  if (opts.Tangle)
  {
    for (vtkm::Id i = 0; i < opts.NumTangle; i++)
    {
      vtkm::source::Tangle tangle(opts.TangleDims);
      dataSets.AppendPartition(tangle.Execute());
    }
  }
  //dataSets = ReadVisItFile(opts.DataPath, opts.DataFile);

  vtkm::filter::contour::Contour contour;
  contour.SetRunMultiThreadedFilter(opts.ThreadMode == "task");
  contour.SetThreadsPerCPU(numCPUThreads);
  contour.SetThreadsPerGPU(numGPUThreads);
  contour.SetGenerateNormals(true);
  contour.SetActiveField(opts.Field);

  if (!opts.IsoValues.empty())
  {
    for (std::size_t i = 0; i < opts.IsoValues.size(); i++)
      contour.SetIsoValue((vtkm::Id)i, opts.IsoValues[i]);
  }
  else
  {
    auto field = dataSets.GetPartition(0).GetField(opts.Field);
    vtkm::Range range;
    field.GetRange(&range);
    vtkm::FloatDefault dR = (range.Max-range.Min) / static_cast<vtkm::FloatDefault>(opts.IsoLevels+1);
    vtkm::FloatDefault v = dR;
    for (int i = 0; i < opts.IsoLevels; i++)
    {
      contour.SetIsoValue((vtkm::Id)i, v);
      v += dR;
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  auto result = contour.Execute(dataSets);
  auto t2 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;

  return 0;
}
