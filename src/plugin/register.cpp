#include <FECore/FECoreKernel.h>
#include "version.h"
#include "task/vfm_task.h"

FECORE_EXPORT unsigned int GetSDKVersion()
{
    return FE_SDK_VERSION;
}

FECORE_EXPORT void GetPluginVersion(int &major, int &minor, int &patch)
{
    major = VERSION;
    minor = SUBVERSION;
    patch = SUBSUBVERSION;
}

FECORE_EXPORT void PluginInitialize(FECoreKernel &febio)
{
    FECoreKernel::SetInstance(&febio);
    REGISTER_FECORE_CLASS(VFMTask, "VFM");
}
