//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtkm_rendering_compositing_PayloadCompositor_h
#define vtkm_rendering_compositing_PayloadCompositor_h

#include <sstream>
#include <vtkm/rendering/vtkm_rendering_export.h>
#include <vtkm/rendering/vtkh/compositing/PayloadImage.hpp>

namespace vtkh
{

class VTKM_RENDERING_EXPORT PayloadCompositor
{
public:
    PayloadCompositor();

    void ClearImages();

    void AddImage(PayloadImage &image);

    PayloadImage Composite();
protected:
    std::vector<PayloadImage>  m_images;
};

};

#endif //vtkm_rendering_compositing_PayloadCompoistor_h
