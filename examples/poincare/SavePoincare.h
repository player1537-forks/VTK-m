#ifndef SavePoincare_h
#define SavePoincare_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>

#include "XGCParameters.h"


void
SavePoincare(const vtkm::cont::DataSet& ds,
             const XGCParameters& xgcParams,
             std::map<std::string, std::vector<std::string>>& args);
#endif
