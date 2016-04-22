
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

#ifndef PBRT_INTEGRATORS_TOFPATH_H
#define PBRT_INTEGRATORS_TOFPATH_H
#include "stdafx.h"

// integrators/tofpath.h*
#include "pbrt.h"
#include "integrator.h"

// ToFPathIntegrator Declarations
class ToFPathIntegrator : public SamplerIntegrator {
  public:
    // ToFPathIntegrator Public Methods
    IntegrationResult Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;
	ToFPathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                   std::shared_ptr<Sampler> sampler)
        : SamplerIntegrator(camera, sampler), maxDepth(maxDepth) {}

  private:
    // ToFPathIntegrator Private Data
    const int maxDepth;
};

ToFPathIntegrator *CreateToFPathIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);

#endif  // PBRT_INTEGRATORS_TOFPATH_H
