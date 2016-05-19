
/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PATHTOF_H
#define PBRT_INTEGRATORS_PATHTOF_H
#include "stdafx.h"

// integrators/pathtof.h*
#include "pbrt.h"
#include "integrator.h"

// PathToFIntegrator Declarations
class PathToFIntegrator : public SamplerIntegrator {
  public:
    // PathToFIntegrator Public Methods
    IntegrationResult Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;
	PathToFIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                   std::shared_ptr<Sampler> sampler)
        : SamplerIntegrator(camera, sampler), maxDepth(maxDepth) {}

  private:
    // PathToFIntegrator Private Data
    const int maxDepth;
};

PathToFIntegrator *CreatePathToFIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);

#endif  // PBRT_INTEGRATORS_PATHTOF_H
