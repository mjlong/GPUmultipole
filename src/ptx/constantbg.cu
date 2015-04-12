#include <optix_world.h>
#include "commonStructs.h"
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );

RT_PROGRAM void miss()
{
#if defined(__MANY__)
  prd.out = 1;
#endif
}
