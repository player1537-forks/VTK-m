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

#include <vtkm/rendering_new/Scene.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

namespace vtkm
{
namespace rendering_new
{

void Scene::SetRenderBatchSize(int batch_size)
{
  if (batch_size < 1)
  {
    throw vtkm::cont::ErrorBadValue("Render batch size must be greater than 0");
  }
  this->RenderInBatches = true;
  this->BatchSize = batch_size;
}

void Scene::Render()
{
  if (this->RenderInBatches)
  {
    for (auto& plot : this->Plots)
      plot.Render();
  }
  else
  {
    for (auto& plot : this->Plots)
      plot.Render();
  }
}

//
//old stuff
//
#if 0
void
Scene::Render()
{
  std::vector<vtkm::Range> ranges;
  std::vector<std::string> fieldNames;
  std::vector<vtkm::cont::ColorTable> colorTables;
  bool firstTime = true;

  //ignore cinema right now.....
  std::size_t numPlots = this->Plots.size();
  for (std::size_t i = 0; i < numPlots; i++)
  {
    auto& plot = this->Plots[i];
    plot.GetCanvas().Clear();
    auto renderer = plot.Renderers[0];

    //opaque pass...
    //last plot
    bool lastPlot = (i == numPlots-1);
    renderer->SetDoComposite(lastPlot);
    //renderer->SetPlots({plot});

    renderer->Update(plot);
    //renderer->ClearPlots();

    if (firstTime)
    {
      const auto& r = plot.Renderers[0];
      if (r->GetHasColorTable())
      {
        ranges.push_back(r->GetScalarRange());
        fieldNames.push_back(r->GetFieldName());
        colorTables.push_back(r->GetColorTable());
      }
      firstTime = false;
    }

    //batch
    plot.RenderWorldAnnotations();
    plot.RenderScreenAnnotations(fieldNames, ranges, colorTables);
    plot.RenderBackground();
    plot.Save();
  }
}

void
Scene::RenderORIG()
{
  std::vector<vtkm::Range> ranges;
  std::vector<std::string> field_names;
  std::vector<vtkm::cont::ColorTable> color_tables;
  bool do_once = true;

  //
  // We are going to render images in batches. With databases
  // like Cinema, we could be rendering hundres of images. Keeping
  // all the canvases around can hog memory so we will conserve it.
  // For example, if we rendered 360 images at 1024^2, all the canvases
  // would consume 7GB of space. Not good on the GPU, where resources
  // are limited.
  //
  const int render_size = this->Plots.size();
  int batch_start = 0;
  while(batch_start < render_size)
  {
    int batch_end = std::min(this->BatchSize + batch_start, render_size);
    auto begin = this->Plots.begin() + batch_start;
    auto end = this->Plots.begin() + batch_end;

    std::vector<vtkm::rendering_new::Plot> current_batch(begin, end);

    for(auto plot : current_batch)
      plot.GetCanvas().Clear();

    const int numRenderers = this->Renderers.size();
    auto renderer = this->Renderers.begin();

    // render order is enforced inside add
    // Order is:
    // 1) surfaces
    // 2) meshes
    // 3) volume

    // if we have both surfaces/mesh and volumes
    // we need to synchronize depths so that volume
    // only render to the max depth
    bool synch_depths = false;

    int opaqueRenderers = numRenderers;
    if(this->HasVolume)
    {
      opaqueRenderers -= 1;
    }

    //
    // pass 1: opaque geometry
    //
    for(int i = 0; i < opaqueRenderers; ++i)
    {
      if(i == opaqueRenderers - 1)
      {
        (*renderer)->SetDoComposite(true);
      }
      else
      {
        (*renderer)->SetDoComposite(false);
      }

      //(*renderer)->SetPlots(current_batch);
      (*renderer)->Update(current_batch);

      //(*renderer)->ClearPlots();

      synch_depths = true;
      renderer++;
    }

    //
    // pass 2: volume
    //
    if(this->HasVolume)
    {
      if(synch_depths)
      {
        SynchDepths(current_batch);
      }
      (*renderer)->SetDoComposite(true);
      //(*renderer)->SetPlots(current_batch);
      (*renderer)->Update(current_batch);

      //current_batch  = (*renderer)->GetPlots();
      //(*renderer)->ClearPlots();
    }

    if(do_once)
    {
      // gather color tables and other information for
      // annotations
      for (const auto& r : this->Renderers)
      {
        if((*r).GetHasColorTable())
        {
          ranges.push_back((*r).GetScalarRange());
          field_names.push_back((*r).GetFieldName());
          color_tables.push_back((*r).GetColorTable());
        }
      }
      do_once = false;
    }

    // render screen annotations last and save
    for(std::size_t i = 0; i < current_batch.size(); ++i)
    {
      current_batch[i].RenderWorldAnnotations();
      current_batch[i].RenderScreenAnnotations(field_names, ranges, color_tables);
      current_batch[i].RenderBackground();
      current_batch[i].Save();
    }

    batch_start = batch_end;
  } // while
}

void
Scene::AddRenderer(vtkm::rendering_new::Renderer *renderer)
{
  bool is_volume = this->IsVolume(renderer);
  bool is_mesh = this->IsMesh(renderer);

  if(is_volume)
  {
    if(this->HasVolume)
    {
      throw vtkm::cont::ErrorBadValue("Scenes only support a single volume plot");
    }

    this->HasVolume = true;
    // make sure that the volume render is last
    this->Renderers.push_back(renderer);
  }
  else if(is_mesh)
  {
    // make sure that the mesh plot is last
    // and before the volume plot
    if(this->HasVolume)
    {
      if(this->Renderers.size() == 1)
      {
        this->Renderers.push_front(renderer);
      }
      else
      {
        auto it = this->Renderers.end();
        it--;
        it--;
        this->Renderers.insert(it,renderer);
      }
    }
    else
    {
      this->Renderers.push_back(renderer);
    }
  }
  else
  {
    this->Renderers.push_front(renderer);
  }
}

bool
Scene::IsMesh(vtkm::rendering_new::Renderer *renderer)
{
  bool is_mesh = false;

  if(dynamic_cast<vtkm::rendering_new::MeshRenderer*>(renderer) != nullptr)
  {
    is_mesh = true;
  }
  return is_mesh;
}

bool
Scene::IsVolume(vtkm::rendering_new::Renderer *renderer)
{
  bool is_volume = false;

  if(dynamic_cast<vtkm::rendering_new::VolumeRenderer*>(renderer) != nullptr)
  {
    is_volume = true;
  }
  return is_volume;
}


void Scene::SynchDepths(std::vector<vtkm::rendering_new::Plot> &plots)
{
  for (auto plot : plots)
    plot.SyncDepth();
}
#endif

}
} // namespace vtkm::rendering_new
