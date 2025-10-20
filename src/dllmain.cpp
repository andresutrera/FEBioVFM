// Minimal entry points required by FEBio to discover this plugin. We only
// register the VFM task for now so the engine can be invoked via
//   febio -i file.feb -task="VFM" VFMData.feb

#include <FECore/FECoreKernel.h>

// #include "VFM.h"
#include "task/vfm_task.h"

FECORE_EXPORT unsigned int GetSDKVersion()
{
    return FE_SDK_VERSION;
}

FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
    FECoreKernel::SetInstance(&febio);
    REGISTER_FECORE_CLASS(VFMTask, "VFM2");
}

FECORE_EXPORT void PluginCleanup()
{
    // No resources to release yet.
}

