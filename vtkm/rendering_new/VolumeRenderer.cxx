//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <memory>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering_new/VolumeRenderer.h>
#include <vtkm/rendering_new/compositing/Compositor.h>
#include <vtkm/thirdparty/diy/diy.h>

#ifdef VTKM_ENABLE_MPI
#include <mpi.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#endif

#include <vtkm/cont/ColorTable.h>
#include <vtkm/rendering/ConnectivityProxy.h>
#include <vtkm/rendering/raytracing/Camera.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <vtkm/rendering/raytracing/VolumeRendererStructured.h>
#include <vtkm/rendering_new/compositing/PartialCompositor.h>

#include <vtkm/rendering_new/compositing/VolumePartial.h>

#define VTKH_OPACITY_CORRECTION 10.f

namespace vtkm
{
namespace rendering_new
{

namespace detail
{

struct VisOrdering
{
  int m_rank;
  int m_domain_index;
  int m_order;
  float m_minz;
};

struct DepthOrder
{
  inline bool operator()(const VisOrdering& lhs, const VisOrdering& rhs)
  {
    return lhs.m_minz < rhs.m_minz;
  }
};

struct RankOrder
{
  inline bool operator()(const VisOrdering& lhs, const VisOrdering& rhs)
  {
    if (lhs.m_rank < rhs.m_rank)
    {
      return true;
    }
    else if (lhs.m_rank == rhs.m_rank)
    {
      return lhs.m_domain_index < rhs.m_domain_index;
    }
    return false;
  }
};

vtkm::cont::ArrayHandle<vtkm::Vec4f_32> convert_table(const vtkm::cont::ColorTable& colorTable)
{

  constexpr vtkm::Float32 conversionToFloatSpace = (1.0f / 255.0f);

  vtkm::cont::ArrayHandle<vtkm::Vec4ui_8> temp;

  {
    vtkm::cont::ScopedRuntimeDeviceTracker tracker(vtkm::cont::DeviceAdapterTagSerial{});
    colorTable.Sample(1024, temp);
  }

  vtkm::cont::ArrayHandle<vtkm::Vec4f_32> color_map;
  color_map.Allocate(1024);
  auto portal = color_map.WritePortal();
  auto colorPortal = temp.ReadPortal();
  for (vtkm::Id i = 0; i < 1024; ++i)
  {
    auto color = colorPortal.Get(i);
    vtkm::Vec4f_32 t(color[0] * conversionToFloatSpace,
                     color[1] * conversionToFloatSpace,
                     color[2] * conversionToFloatSpace,
                     color[3] * conversionToFloatSpace);
    portal.Set(i, t);
  }
  return color_map;
}

class VolumeWrapper
{
protected:
  const vtkm::cont::DataSet& DataSet;
  vtkm::Range ScalarRange;
  std::string FieldName;
  vtkm::Float32 SampleDist;
  vtkm::cont::ArrayHandle<vtkm::Vec4f_32> ColorMap;

public:
  VolumeWrapper() = delete;

  VolumeWrapper(const vtkm::cont::DataSet& data_set)
    : DataSet(data_set)
  {
  }

  virtual ~VolumeWrapper() {}

  void sample_distance(const vtkm::Float32& distance) { this->SampleDist = distance; }

  void field(const std::string& field_name) { this->FieldName = field_name; }

  void scalar_range(const vtkm::Range& range) { this->ScalarRange = range; }

  void color_map(vtkm::cont::ArrayHandle<vtkm::Vec4f_32>& color_map) { this->ColorMap = color_map; }

  virtual void render(const vtkm::rendering::Camera& camera,
                      vtkm::rendering::CanvasRayTracer& canvas,
                      std::vector<VolumePartial<float>>& partials) = 0;
};

void vtkm_to_partials(vtkm::rendering::PartialVector32& vtkm_partials,
                      std::vector<VolumePartial<float>>& partials)
{
  const int num_vecs = vtkm_partials.size();
  std::vector<int> offsets;
  offsets.reserve(num_vecs);

  int total_size = 0;
  for (int i = 0; i < num_vecs; ++i)
  {
    const int size = vtkm_partials[i].PixelIds.GetNumberOfValues();
    offsets.push_back(total_size);
    total_size += size;
  }

  partials.resize(total_size);

  for (int i = 0; i < num_vecs; ++i)
  {
    const int size = vtkm_partials[i].PixelIds.GetNumberOfValues();
    auto pixel_ids = vtkm_partials[i].PixelIds.ReadPortal();
    auto distances = vtkm_partials[i].Distances.ReadPortal();
    auto colors = vtkm_partials[i].Buffer.Buffer.ReadPortal();

    const int offset = offsets[i];
#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int p = 0; p < size; ++p)
    {
      VolumePartial<float>& partial = partials[offset + p];
      partial.Pixel[0] = colors.Get(p * 4 + 0);
      partial.Pixel[1] = colors.Get(p * 4 + 1);
      partial.Pixel[2] = colors.Get(p * 4 + 2);
      partial.Alpha = colors.Get(p * 4 + 3);
      partial.PixelId = pixel_ids.Get(p);
      partial.Depth = distances.Get(p);
    }
  }
}

class UnstructuredWrapper : public VolumeWrapper
{
  vtkm::rendering::ConnectivityProxy m_tracer;

public:
  UnstructuredWrapper(const vtkm::cont::DataSet& data_set)
    : VolumeWrapper(data_set)
    , m_tracer(data_set, "") //DRP: Fix me!!  fieldname required in ConnectivityTracer
  {
  }

  virtual void render(const vtkm::rendering::Camera& camera,
                      vtkm::rendering::CanvasRayTracer& canvas,
                      std::vector<VolumePartial<float>>& partials) override
  {
    const vtkm::cont::CoordinateSystem& coords = this->DataSet.GetCoordinateSystem();

    vtkm::rendering::raytracing::Camera rayCamera;
    vtkm::rendering::raytracing::Ray<vtkm::Float32> rays;
    vtkm::Int32 width = (vtkm::Int32)canvas.GetWidth();
    vtkm::Int32 height = (vtkm::Int32)canvas.GetHeight();

    rayCamera.SetParameters(camera, width, height);

    rayCamera.CreateRays(rays, coords.GetBounds());
    rays.Buffers.at(0).InitConst(0.f);
    vtkm::rendering::raytracing::RayOperations::MapCanvasToRays(rays, camera, canvas);

    m_tracer.SetSampleDistance(this->SampleDist);
    m_tracer.SetColorMap(this->ColorMap);
    m_tracer.SetScalarField(this->FieldName);
    m_tracer.SetScalarRange(this->ScalarRange);

    vtkm::rendering::PartialVector32 vtkm_partials;
    vtkm_partials = m_tracer.PartialTrace(rays);

    vtkm_to_partials(vtkm_partials, partials);
  }
};

class StructuredWrapper : public VolumeWrapper
{
public:
  StructuredWrapper(const vtkm::cont::DataSet& data_set)
    : VolumeWrapper(data_set)
  {
  }
  virtual void render(const vtkm::rendering::Camera& camera,
                      vtkm::rendering::CanvasRayTracer& canvas,
                      std::vector<VolumePartial<float>>& partials) override
  {
    const vtkm::cont::UnknownCellSet& cellset = this->DataSet.GetCellSet();
    const vtkm::cont::Field& field = this->DataSet.GetField(this->FieldName);
    const vtkm::cont::CoordinateSystem& coords = this->DataSet.GetCoordinateSystem();

    vtkm::rendering::raytracing::Camera rayCamera;
    vtkm::rendering::raytracing::Ray<vtkm::Float32> rays;
    vtkm::Int32 width = (vtkm::Int32)canvas.GetWidth();
    vtkm::Int32 height = (vtkm::Int32)canvas.GetHeight();
    rayCamera.SetParameters(camera, width, height);

    rayCamera.CreateRays(rays, coords.GetBounds());
    rays.Buffers.at(0).InitConst(0.f);
    vtkm::rendering::raytracing::RayOperations::MapCanvasToRays(rays, camera, canvas);

    vtkm::rendering::raytracing::VolumeRendererStructured tracer;
    tracer.SetSampleDistance(this->SampleDist);
    tracer.SetData(
      coords, field, cellset.AsCellSet<vtkm::cont::CellSetStructured<3>>(), this->ScalarRange);
    tracer.SetColorMap(this->ColorMap);


    tracer.Render(rays);

    // Convert the rays to partial composites
    const int ray_size = rays.NumRays;
    // partials use the max distance
    auto depths = rays.MaxDistance.ReadPortal();
    auto pixel_ids = rays.PixelIdx.ReadPortal();
    auto colors = rays.Buffers.at(0).Buffer.ReadPortal();

    // TODO: better way? we could do this in parallel if we
    // don't check the alpha
    partials.reserve(ray_size);
    for (int i = 0; i < ray_size; ++i)
    {
      const int offset = i * 4;
      float alpha = colors.Get(offset + 3);
      if (alpha < 0.001f)
        continue;
      VolumePartial<float> partial;
      partial.Pixel[0] = colors.Get(offset + 0);
      partial.Pixel[1] = colors.Get(offset + 1);
      partial.Pixel[2] = colors.Get(offset + 2);
      partial.Alpha = alpha;
      partial.PixelId = pixel_ids.Get(i);
      partial.Depth = depths.Get(i);
      partials.push_back(std::move(partial));
    }
  }
};

void partials_to_canvas(std::vector<VolumePartial<float>>& partials,
                        const vtkm::rendering::Camera& camera,
                        vtkm::rendering::CanvasRayTracer& canvas)
{

  // partial depths are in world space but the canvas depths
  // are in image space. We have to find the intersection
  // point to project it into image space to get the correct
  // depths for annotations
  vtkm::Id width = canvas.GetWidth();
  vtkm::Id height = canvas.GetHeight();
  vtkm::Matrix<vtkm::Float32, 4, 4> projview =
    vtkm::MatrixMultiply(camera.CreateProjectionMatrix(width, height), camera.CreateViewMatrix());

  const vtkm::Vec3f_32 origin = camera.GetPosition();

  float fov_y = camera.GetFieldOfView();
  float fov_x = fov_y;
  if (width != height)
  {
    vtkm::Float32 fovyRad = fov_y * vtkm::Pi_180f();
    vtkm::Float32 verticalDistance = vtkm::Tan(0.5f * fovyRad);
    vtkm::Float32 aspectRatio = vtkm::Float32(width) / vtkm::Float32(height);
    vtkm::Float32 horizontalDistance = aspectRatio * verticalDistance;
    vtkm::Float32 fovxRad = 2.0f * vtkm::ATan(horizontalDistance);
    fov_x = fovxRad / vtkm::Pi_180f();
  }

  vtkm::Vec3f_32 look = camera.GetLookAt() - origin;
  vtkm::Normalize(look);
  vtkm::Vec3f_32 up = camera.GetViewUp();

  const vtkm::Float32 thx = tanf((fov_x * vtkm::Pi_180f()) * .5f);
  const vtkm::Float32 thy = tanf((fov_y * vtkm::Pi_180f()) * .5f);
  vtkm::Vec3f_32 ru = vtkm::Cross(look, up);

  vtkm::Normalize(ru);
  vtkm::Vec3f_32 rv = vtkm::Cross(ru, look);
  vtkm::Normalize(rv);
  vtkm::Vec3f_32 delta_x = ru * (2 * thx / (float)width);
  vtkm::Vec3f_32 delta_y = ru * (2 * thy / (float)height);

  vtkm::Float32 zoom = camera.GetZoom();
  if (zoom > 0)
  {
    delta_x[0] = delta_x[0] / zoom;
    delta_x[1] = delta_x[1] / zoom;
    delta_x[2] = delta_x[2] / zoom;
    delta_y[0] = delta_y[0] / zoom;
    delta_y[1] = delta_y[1] / zoom;
    delta_y[2] = delta_y[2] / zoom;
  }

  const int size = partials.size();
  auto colors = canvas.GetColorBuffer().WritePortal();
  auto depths = canvas.GetDepthBuffer().WritePortal();

#ifdef VTKH_OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int p = 0; p < size; ++p)
  {
    const int pixel_id = partials[p].PixelId;
    const int i = pixel_id % width;
    const int j = pixel_id / width;

    vtkm::Vec3f_32 dir;
    dir = look + delta_x * ((2.f * float(i) - float(width)) / 2.0f) +
      delta_y * ((2.f * float(j) - float(height)) / 2.0f);
    vtkm::Normalize(dir);

    const float world_depth = partials[p].Depth;

    vtkm::Vec3f_32 pos = origin + world_depth * dir;
    vtkm::Vec4f_32 point(pos[0], pos[1], pos[2], 1.f);
    vtkm::Vec4f_32 newpoint;
    newpoint = vtkm::MatrixMultiply(projview, point);

    // don't push it all the way(.49 instead of .5) so that
    // subtle differences allow bounding box annotations don't
    // draw in front of the back
    const float image_depth = 0.5f * (newpoint[2] / newpoint[3]) + 0.49f;

    vtkm::Vec4f_32 color;
    color[0] = partials[p].Pixel[0];
    color[1] = partials[p].Pixel[1];
    color[2] = partials[p].Pixel[2];
    color[3] = partials[p].Alpha;

    vtkm::Vec4f_32 inColor = colors.Get(pixel_id);
    // We crafted the rendering so that all new colors are in front
    // of the colors that exist in the canvas
    // if transparency exists, all alphas have been pre-multiplied
    vtkm::Float32 alpha = (1.f - color[3]);
    color[0] = color[0] + inColor[0] * alpha;
    color[1] = color[1] + inColor[1] * alpha;
    color[2] = color[2] + inColor[2] * alpha;
    color[3] = inColor[3] * alpha + color[3];

    colors.Set(pixel_id, color);
    depths.Set(pixel_id, image_depth);
  }
}

} //  namespace detail

VolumeRenderer::VolumeRenderer()
{
  typedef vtkm::rendering::MapperVolume TracerType;
  this->Tracer = std::make_shared<TracerType>();
  this->Mapper = this->Tracer;
  this->Tracer->SetCompositeBackground(false);
  //
  // add some default opacity to the color table
  //

  std::cout << __FILE__ << " " << __LINE__ << "  Do we need this????" << std::endl;
  //this->ColorTable.AddPointAlpha(0.0f, .02);
  //this->ColorTable.AddPointAlpha(.0f, .5);
  this->NumSamples = 100.f;
  this->HasUnstructured = false;
}

VolumeRenderer::~VolumeRenderer()
{
  ClearWrappers();
}

void VolumeRenderer::Update(vtkm::rendering_new::Plot& plot)
{
  //DRP: Logger
  //VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  //DRP: Logger
  //VTKH_DATA_ADD("device", GetCurrentDevice());
  long long int in_cells = this->Input->GetNumberOfCells();
  //VTKH_DATA_ADD("input_cells", in_cells);
  //VTKH_DATA_ADD("input_domains", this->Input->GetNumberOfDomains());
  int in_topo_dims;
  bool in_structured = this->Input->IsStructured(in_topo_dims);
  if (in_structured)
  {
    //VTKH_DATA_ADD("in_topology", "structured");
  }
  else
  {
    //VTKH_DATA_ADD("in_topology", "unstructured");
  }
#endif

  PreExecute(plot);
  DoExecute(plot);
  PostExecute(plot);

  //DRP: Logger
  //VTKH_DATA_CLOSE();
}

void VolumeRenderer::CorrectOpacity()
{
  const float correction_scalar = VTKH_OPACITY_CORRECTION;
  float samples = this->NumSamples;

  float ratio = correction_scalar / samples;
  vtkm::cont::ColorTable corrected;
  corrected = this->GetColorTable().MakeDeepCopy();
  int num_points = corrected.GetNumberOfPointsAlpha();
  for (int i = 0; i < num_points; i++)
  {
    vtkm::Vec<vtkm::Float64, 4> point;
    corrected.GetPointAlpha(i, point);
    point[1] = 1. - vtkm::Pow((1. - point[1]), double(ratio));
    corrected.UpdatePointAlpha(i, point);
  }

  this->CorrectedColorTable = corrected;
}

void VolumeRenderer::DoExecute(vtkm::rendering_new::Plot& plot)
{
  bool localHasOnePartition = this->Actor.GetDataSet().GetNumberOfPartitions();
  bool globalHasOnePartition = localHasOnePartition;
#ifdef VTKM_ENABLE_MPI
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  vtkmdiy::mpi::all_reduce(
    comm, localHasOnePartition, globalHasOnePartition, std::logical_and<bool>());
#endif

  if (globalHasOnePartition && !this->HasUnstructured)
  {
    // Danger: this logic only works if there is exactly one per rank
    RenderOneDomainPerRank(plot);
  }
  else
  {
    RenderMultipleDomainsPerRank(plot);
  }
}

void VolumeRenderer::RenderOneDomainPerRank(vtkm::rendering_new::Plot& plot)
{
  if (this->Mapper.get() == 0)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class";
    throw vtkm::cont::ErrorBadValue(msg);
  }

  this->Tracer->SetSampleDistance(this->SampleDist);

  int num_domains = static_cast<int>(this->Actor.GetDataSet().GetNumberOfPartitions());
  if (num_domains > 1)
  {
    throw vtkm::cont::ErrorBadValue("RenderOneDomainPerRank: this should never happend.");
  }

  for (int dom = 0; dom < num_domains; ++dom)
  {
    vtkm::cont::DataSet data_set = this->Actor.GetDataSet().GetPartition(dom);

    if (!data_set.HasField(this->GetFieldName()))
    {
      continue;
    }

    const vtkm::cont::UnknownCellSet& cellset = data_set.GetCellSet();
    const vtkm::cont::Field& field = data_set.GetField(this->GetFieldName());
    const vtkm::cont::CoordinateSystem& coords = data_set.GetCoordinateSystem();

    if (cellset.GetNumberOfCells() == 0)
      continue;

    this->Mapper->SetActiveColorTable(this->CorrectedColorTable);

    Plot::vtkmCanvas& canvas = plot.GetCanvas();
    const auto& camera = plot.GetCamera();
    this->Mapper->SetCanvas(&canvas);
    this->Mapper->RenderCells(
      cellset, coords, field, this->CorrectedColorTable, camera, this->Actor.GetScalarRange());
  }

  if (this->DoComposite)
    this->Composite(plot);
}

void VolumeRenderer::RenderMultipleDomainsPerRank(vtkm::rendering_new::Plot& plot)
{
  // We are treating this as the most general case
  // where we could have a mix of structured and
  // unstructured data sets. There are zero
  // assumptions

  // this might be smaller than the input since
  // it is possible for cell sets to be empty
  const int num_domains = this->Wrappers.size();

  vtkm::cont::ArrayHandle<vtkm::Vec4f_32> color_map =
    detail::convert_table(this->CorrectedColorTable);
  vtkm::cont::ArrayHandle<vtkm::Vec4f_32> color_map2 = detail::convert_table(this->GetColorTable());

  // render/domain/result
  std::vector<std::vector<VolumePartial<float>>> render_partials;
  render_partials.resize(num_domains);

  for (int i = 0; i < num_domains; ++i)
  {
    detail::VolumeWrapper* wrapper = this->Wrappers[i];
    wrapper->sample_distance(this->SampleDist);
    wrapper->color_map(color_map);
    wrapper->field(this->GetFieldName());
    wrapper->scalar_range(this->GetScalarRange());

    Plot::vtkmCanvas& canvas = plot.GetCanvas();
    const auto& camera = plot.GetCamera();
    wrapper->render(camera, canvas, render_partials[i]);
  }

  PartialCompositor<VolumePartial<float>> compositor;
#ifdef VTKM_ENABLE_MPI
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MPI_Comm comm = vtkmdiy::mpi::mpi_cast(diy_comm.handle());
  compositor.set_comm_handle(comm);
#endif

  // composite
  std::vector<VolumePartial<float>> res;
  compositor.composite(render_partials, res);
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() == 0)
  {
    detail::partials_to_canvas(res, plot.GetCamera(), plot.GetCanvas());
  }
}

void VolumeRenderer::PreExecute(vtkm::rendering_new::Plot& plot)
{
  Renderer::PreExecute(plot);

  CorrectOpacity();

  vtkm::Vec<vtkm::Float32, 3> extent;
  extent[0] = static_cast<vtkm::Float32>(this->Actor.GetSpatialBounds().X.Length());
  extent[1] = static_cast<vtkm::Float32>(this->Actor.GetSpatialBounds().Y.Length());
  extent[2] = static_cast<vtkm::Float32>(this->Actor.GetSpatialBounds().Z.Length());
  vtkm::Float32 dist = vtkm::Magnitude(extent) / this->NumSamples;
  this->SampleDist = dist;
}

void VolumeRenderer::PostExecute(vtkm::rendering_new::Plot& vtkmNotUsed(plot))
{
  // do nothing and override compositing since
  // we already did it
}

void VolumeRenderer::SetNumberOfSamples(const int num_samples)
{
  if (num_samples < 1)
  {
    throw vtkm::cont::ErrorBadValue("Volume rendering samples must be greater than 0");
  }
  this->NumSamples = num_samples;
}

std::shared_ptr<vtkm::rendering::Canvas> VolumeRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

float VolumeRenderer::FindMinDepth(const vtkm::rendering::Camera& camera,
                                   const vtkm::Bounds& bounds) const
{

  vtkm::Vec<vtkm::Float64, 3> center = bounds.Center();
  vtkm::Vec<vtkm::Float64, 3> fcenter;
  fcenter[0] = static_cast<vtkm::Float32>(center[0]);
  fcenter[1] = static_cast<vtkm::Float32>(center[1]);
  fcenter[2] = static_cast<vtkm::Float32>(center[2]);
  vtkm::Vec<vtkm::Float32, 3> pos = camera.GetPosition();
  vtkm::Float32 dist = vtkm::Magnitude(fcenter - pos);
  return dist;
}

void VolumeRenderer::Composite(vtkm::rendering_new::Plot& plot)
{
  this->Compositor->SetCompositeMode(Compositor::VIS_ORDER_BLEND);
  FindVisibilityOrdering(plot);

  float* color_buffer = &(plot.GetCanvas().GetColorBuffer().WritePortal().GetArray()[0][0]);
  float* depth_buffer = plot.GetCanvas().GetDepthBuffer().WritePortal().GetArray();
  int height = plot.GetCanvas().GetHeight();
  int width = plot.GetCanvas().GetWidth();

  this->Compositor->AddImage(color_buffer, depth_buffer, width, height, this->VisibilityOrders[0]);

  Image result = this->Compositor->Composite();
#ifdef VTKM_ENABLE_MPI
  if (vtkm::cont::EnvironmentTracker::GetCommunicator().rank() == 0)
  {
#endif
    ImageToCanvas(result, plot.GetCanvas(), true);
#ifdef VTKM_ENABLE_MPI
  }
#endif
  this->Compositor->ClearImages();
}

void VolumeRenderer::DepthSort(int num_domains,
                               std::vector<float>& min_depths,
                               std::vector<int>& local_vis_order)
{
  if (min_depths.size() != static_cast<std::size_t>(num_domains))
  {
    throw vtkm::cont::ErrorBadValue("min depths size does not equal the number of domains");
  }
  if (local_vis_order.size() != static_cast<std::size_t>(num_domains))
  {
    throw vtkm::cont::ErrorBadValue("local vis order not equal to number of domains");
  }
#ifdef VTKM_ENABLE_MPI
  int root = 0;
  auto diy_comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  MPI_Comm comm = vtkmdiy::mpi::mpi_cast(diy_comm.handle());
  int num_ranks = diy_comm.size();
  int rank = diy_comm.rank();
  int* domain_counts = NULL;
  int* domain_offsets = NULL;
  int* vis_order = NULL;
  float* depths = NULL;

  if (rank == root)
  {
    domain_counts = new int[num_ranks];
    domain_offsets = new int[num_ranks];
  }

  MPI_Gather(&num_domains, 1, MPI_INT, domain_counts, 1, MPI_INT, root, comm);

  int depths_size = 0;
  if (rank == root)
  {
    //scan for dispacements
    domain_offsets[0] = 0;
    for (int i = 1; i < num_ranks; ++i)
    {
      domain_offsets[i] = domain_offsets[i - 1] + domain_counts[i - 1];
    }

    for (int i = 0; i < num_ranks; ++i)
    {
      depths_size += domain_counts[i];
    }

    depths = new float[depths_size];
  }

  MPI_Gatherv(&min_depths[0],
              num_domains,
              MPI_FLOAT,
              depths,
              domain_counts,
              domain_offsets,
              MPI_FLOAT,
              root,
              comm);

  if (rank == root)
  {
    std::vector<detail::VisOrdering> order;
    order.resize(depths_size);

    for (int i = 0; i < num_ranks; ++i)
    {
      for (int c = 0; c < domain_counts[i]; ++c)
      {
        int index = domain_offsets[i] + c;
        order[index].m_rank = i;
        order[index].m_domain_index = c;
        order[index].m_minz = depths[index];
      }
    }

    std::sort(order.begin(), order.end(), detail::DepthOrder());

    for (int i = 0; i < depths_size; ++i)
    {
      order[i].m_order = i;
    }

    std::sort(order.begin(), order.end(), detail::RankOrder());

    vis_order = new int[depths_size];
    for (int i = 0; i < depths_size; ++i)
    {
      vis_order[i] = order[i].m_order;
    }
  }

  MPI_Scatterv(vis_order,
               domain_counts,
               domain_offsets,
               MPI_INT,
               &local_vis_order[0],
               num_domains,
               MPI_INT,
               root,
               comm);

  if (rank == root)
  {
    delete[] domain_counts;
    delete[] domain_offsets;
    delete[] vis_order;
    delete[] depths;
  }
#else

  std::vector<detail::VisOrdering> order;
  order.resize(num_domains);

  for (int i = 0; i < num_domains; ++i)
  {
    order[i].m_rank = 0;
    order[i].m_domain_index = i;
    order[i].m_minz = min_depths[i];
  }
  std::sort(order.begin(), order.end(), detail::DepthOrder());

  for (int i = 0; i < num_domains; ++i)
  {
    order[i].m_order = i;
  }

  std::sort(order.begin(), order.end(), detail::RankOrder());

  for (int i = 0; i < num_domains; ++i)
  {
    local_vis_order[i] = order[i].m_order;
  }
#endif
}

void VolumeRenderer::FindVisibilityOrdering(vtkm::rendering_new::Plot& plot)
{
  const int num_domains = static_cast<int>(this->Actor.GetDataSet().GetNumberOfPartitions());
  this->VisibilityOrders.resize(num_domains);

  //
  // In order for parallel volume rendering to composite correctly,
  // we nee to establish a visibility ordering to pass to IceT.
  // We will transform the data extents into camera space and
  // take the minimum z value. Then sort them while keeping
  // track of rank, then pass the list in.
  //
  std::vector<float> min_depths;
  min_depths.resize(num_domains);

  const vtkm::rendering::Camera& camera = plot.GetCamera();
  for (int dom = 0; dom < num_domains; ++dom)
  {
    vtkm::Bounds bounds =
      this->Actor.GetDataSet().GetPartition(dom).GetCoordinateSystem().GetBounds();
    min_depths[dom] = FindMinDepth(camera, bounds);
  }

  this->DepthSort(num_domains, min_depths, this->VisibilityOrders);
}

void VolumeRenderer::SetInput(const vtkm::rendering::Actor& actor)
{
  Renderer::SetInput(actor);
  ClearWrappers();

  this->HasUnstructured = false;
  for (const auto& ds : this->Actor.GetDataSet().GetPartitions())
  {
    const vtkm::cont::UnknownCellSet& cellset = ds.GetCellSet();
    if (cellset.GetNumberOfCells() == 0)
    {
      continue;
    }

    const vtkm::cont::CoordinateSystem& coords = ds.GetCoordinateSystem();
    using Uniform = vtkm::cont::ArrayHandleUniformPointCoordinates;
    using DefaultHandle = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
    using Rectilinear =
      vtkm::cont::ArrayHandleCartesianProduct<DefaultHandle, DefaultHandle, DefaultHandle>;
    bool structured = coords.GetData().IsType<Uniform>() || coords.GetData().IsType<Rectilinear>();

    if (structured)
    {
      this->Wrappers.push_back(new detail::StructuredWrapper(ds));
    }
    else
    {
      this->HasUnstructured = true;
      this->Wrappers.push_back(new detail::UnstructuredWrapper(ds));
    }
  }
}

void VolumeRenderer::ClearWrappers()
{
  const int num_wrappers = this->Wrappers.size();
  for (int i = 0; i < num_wrappers; ++i)
  {
    delete this->Wrappers[i];
  }
  this->Wrappers.clear();
}

std::string VolumeRenderer::GetName() const
{
  return "vtkm::rendering_new::VolumeRenderer";
}


}
} // namespace vtkm::rendering_new
