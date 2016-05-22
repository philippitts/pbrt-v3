/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2014 Microsoft Corporation

Author: Phil Pitts
*/

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_DIRECTTOF_H
#define PBRT_INTEGRATORS_DIRECTTOF_H
#include "stdafx.h"

// integrators/directtof.h*
#include "pbrt.h"
#include "integrator.h"
#include "integrators/directlighting.h"
#include "scene.h"

// DirectToFIntegrator Declarations
class DirectToFIntegrator : public SamplerIntegrator {
public:
	// GroundTruthIntegrator Public Methods
	DirectToFIntegrator(LightStrategy strategy,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler)
		: SamplerIntegrator(camera, sampler),
		strategy(strategy) {}
	IntegrationResult Li(const RayDifferential &ray, const Scene &scene,
		Sampler &sampler, MemoryArena &arena, int depth) const;
	void Preprocess(const Scene &scene, Sampler &sampler);

private:
	// DirectToFIntegrator Private Methods
	const LightStrategy strategy;
	std::vector<int> nLightSamples;
};

DirectToFIntegrator *CreateDirectToFIntegrator(
	const ParamSet &params, std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera);


#endif // PBRT_INTEGRATORS_DIRECTTOF_H