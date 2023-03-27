#ifndef VTK_H_RENDERER_LINES_HPP
#define VTK_H_RENDERER_LINES_HPP

#include <vtkm/rendering/vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/vtkm_rendering_export.h>

namespace vtkh {

class VTKM_RENDERING_EXPORT LineRenderer : public Renderer
{
public:
  LineRenderer();
  virtual ~LineRenderer();
  std::string GetName() const override;
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);
  void PreExecute() override;
  void SetRadius(vtkm::Float32 radius);
private:
  bool m_radius_set;
  vtkm::Float32 m_radius;

};

} // namespace vtkh
#endif
