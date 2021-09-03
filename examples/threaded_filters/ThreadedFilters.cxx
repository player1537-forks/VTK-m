//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/filter/Contour.h>

#include <iostream>
#include <chrono>

#include <vtkm/filter/TaskQueue.h>

struct Options
{
public:
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
        {
          this->ThreadMode = "openmp";
          if (a.second.size() != 2) return false;
          this->NumTasks = std::stoi(a.second[1]);
        }
        else if (a.second[0] == "task")
        {
          if (a.second.size() != 2) return false;
          this->NumTasks = std::stoi(a.second[1]);
          this->ThreadMode = "task";
        }
        else if (a.second[0] == "gpu")
        {
          if (a.second.size() != 2) return false;
          this->NumTasks = std::stoi(a.second[1]);
          this->ThreadMode = "gpu";
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
    }

    if (this->MapField == "")
      this->MapField = this->Field;

    if (this->DataFile == "" ||
        this->ThreadMode == "" ||
        this->Field == "")
    {
      if (this->Rank == 0) std::cerr<<"Error in options"<<std::endl;
      return false;
    }

    return true;
  }

  void PrintUsage() const
  {
    std::cerr<<"Usage: --data <dataFile> --field <field> --mapfield <mapfield> --thread <serial openmp task N>"<<std::endl;
  }

  void SetFilter(vtkm::filter::Contour& f) const
  {
    f.SetGenerateNormals(false);
    f.SetActiveField(this->Field);
    f.SetFieldsToPass(this->MapField);
    f.SetNumberOfIsoValues(this->IsoValues.size());
    for (vtkm::Id i = 0; i < (vtkm::Id)this->IsoValues.size(); i++)
      f.SetIsoValue(i, this->IsoValues[i]);
  }
};


int rank = 0, numRanks = 1;
vtkm::cont::PartitionedDataSet
ReadVisItFile(const std::string& dir, const std::string& fname)
{
  std::string tmp = dir + "/" + fname;
  std::cout<<"Opening: "<<tmp<<std::endl;
  std::ifstream visitFile;
  visitFile.open(dir + "/" + fname);

  vtkm::cont::PartitionedDataSet pds;

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

  int n = omp_get_max_threads();
  omp_set_num_threads(20);
#pragma omp parallel for
  for (int i = 0; i < (int)myFiles.size(); i++)
  {
    //std::string f = dir + "/" + myFiles[i];
    //std::cout<<"Reading: "<<f<<std::endl;
    vtkm::io::VTKDataSetReader reader(dir + "/" + myFiles[i]);
    auto ds = reader.ReadDataSet();
#pragma omp critical
    {
      pds.AppendPartition(ds);
    }
  }
  omp_set_num_threads(n);
  std::cout<<"Setting #threads back to "<<n<<std::endl;

  return pds;
}

void
RunTask(vtkm::filter::DataSetQueue *input,
        vtkm::filter::DataSetQueue *output,
        const Options& o,
        int* domCount)
{
  std::pair<vtkm::Id, vtkm::cont::DataSet> entry;

//  std::thread::id this_id = std::this_thread::get_id();
//  std::cout<<"RunTask_"<<this_id<<std::endl;

  while (input->GetTask(entry))
  {
    if (!entry.second.HasField(o.Field))
      continue;

    vtkm::filter::Contour iso;
    o.SetFilter(iso);
    auto res = iso.Execute(entry.second);
    output->Push(std::make_pair(entry.first, std::move(res)));
    (*domCount)++;
  }
}

vtkm::cont::PartitionedDataSet
Run(const vtkm::cont::PartitionedDataSet& input, const Options& o)
{
  vtkm::Id numBlocks = input.GetNumberOfPartitions();

  vtkm::cont::PartitionedDataSet output;

  std::cout<<"Thread Mode= "<<o.ThreadMode<<" : "<<o.NumTasks<<std::endl;

  if (o.ThreadMode == "serial")
  {
    vtkm::filter::Contour iso;
    o.SetFilter(iso);

    auto t1 = std::chrono::high_resolution_clock::now();
    for (vtkm::Id i = 0; i < numBlocks; i++)
    {
      auto out = iso.Execute(input.GetPartition(i));
      output.AppendPartition(out);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout<<"    Timer1= "<<dt.count()<<std::endl;
  }
  else if (o.ThreadMode == "openmp")
  {
    omp_set_num_threads(o.NumTasks);
    std::cout<<"Num omp threads: "<<omp_get_num_threads()<<" "<<omp_get_max_threads()<<std::endl;
    int numThreads = omp_get_max_threads();
    std::vector<int> domCounts(numThreads, 0);

    auto t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(static, 1)
    for (vtkm::Id i = 0; i < numBlocks; i++)
    {
      int tid = omp_get_thread_num();
      domCounts[tid]++;

      vtkm::filter::Contour iso;
      o.SetFilter(iso);
      auto out = iso.Execute(input.GetPartition(i));
#pragma omp critical
      {
        output.AppendPartition(out);
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout<<"    Timer1= "<<dt.count()<<std::endl;

    std::cout<<"Doms per thread: ";
    for (const auto& c : domCounts) std::cout<<c<<" ";
    std::cout<<std::endl;
  }
  else if (o.ThreadMode == "task" || o.ThreadMode == "gpu")
  {
    vtkm::filter::DataSetQueue in(input), out;
    std::vector<std::thread> threads;

    if (o.ThreadMode == "gpu")
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});

    std::vector<int> domCounts(o.NumTasks, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < o.NumTasks; i++)
    {
      threads.push_back(std::thread(RunTask, &in, &out, o, &domCounts[i]));
    }
    for (auto& t : threads)
      t.join();

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout<<"    Timer1= "<<dt.count()<<std::endl;

    output = out.Get();

    std::cout<<"Doms per thread: ";
    for (const auto& c : domCounts) std::cout<<c<<" ";
    std::cout<<std::endl;
  }

  return output;
}


int main(int argc, char** argv)
{
  auto vtkm_opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, vtkm_opts);

  auto& tracker = vtkm::cont::GetRuntimeDeviceTracker();
  const bool runOnCuda = tracker.CanRunOn(vtkm::cont::DeviceAdapterTagCuda{});
  const bool runOnOpenMP = tracker.CanRunOn(vtkm::cont::DeviceAdapterTagOpenMP{});
  const bool runOnTbb = tracker.CanRunOn(vtkm::cont::DeviceAdapterTagTBB{});

  std::cout<<"CanRunOn: C.O.T= "<<runOnCuda<<" "<<runOnOpenMP<<" "<<runOnTbb<<std::endl;

  Options opts;
  if (!opts.ParseOptions(argc, argv))
  {
    opts.PrintUsage();
    return 0;
  }

  auto dataSets = ReadVisItFile(opts.DataPath, opts.DataFile);

  auto t1 = std::chrono::high_resolution_clock::now();
  Run(dataSets, opts);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;


  return 1;
}
