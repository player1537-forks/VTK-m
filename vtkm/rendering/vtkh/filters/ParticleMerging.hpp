#ifndef VTK_H_PARTICLE_MERGING_HPP
#define VTK_H_PARTICLE_MERGING_HPP

#include <vtkm/rendering/vtkm_rendering_export.h>
//#include <vtkh/vtkh.hpp>
#include <vtkm/rendering/vtkh/filters/Filter.hpp>
#include <vtkm/rendering/vtkh/DataSet.hpp>

namespace vtkh
{

class VTKM_RENDERING_EXPORT ParticleMerging : public vtkh::Filter
{
public:
  ParticleMerging();
  virtual ~ParticleMerging();
  std::string GetName() const override;
  void SetField(const std::string &field_name);
  void SetRadius(const vtkm::Float64 radius);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
  vtkm::Float64 m_radius;
};

} //namespace vtkh
#endif
