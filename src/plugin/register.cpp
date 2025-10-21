#include <FECore/FECoreKernel.h>

// #include "VFM.h"
#include "task/vfm_task.h"

FECORE_EXPORT unsigned int GetSDKVersion()
{
    return FE_SDK_VERSION;
}

FECORE_EXPORT void PluginInitialize(FECoreKernel &febio)
{
    FECoreKernel::SetInstance(&febio);
    REGISTER_FECORE_CLASS(VFMTask, "VFM");
}

FECORE_EXPORT void PluginCleanup()
{
    // No resources to release yet.
}
