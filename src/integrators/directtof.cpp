/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2014 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

// integrators/directtof.cpp*
#include "integrators/directtof.h"
#include "paramset.h"
#include "stats.h"

// DirectToFIntegrator Method Definitions
void DirectToFIntegrator::Preprocess(const Scene &scene,
	Sampler &sampler) {
	if (strategy == LightStrategy::UniformSampleAll) {
		// Compute number of samples to use for each light
		for (const auto &light : scene.lights)
			nLightSamples.push_back(sampler.RoundCount(light->nSamples));

		// Request samples for sampling all lights
		for (size_t i = 0; i < scene.lights.size(); ++i) {
			sampler.Request2DArray(nLightSamples[i]);
			sampler.Request2DArray(nLightSamples[i]);
		}
	}
}

IntegrationResult DirectToFIntegrator::Li(const RayDifferential &ray,
	const Scene &scene, Sampler &sampler,
	MemoryArena &arena, int depth) const {
	ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.f);
	// Find closest ray intersection or return nothing on a miss
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		return L;
	}

	// Compute distance from camera to surface interaction
	Float distance = Distance(ray.o, isect.p);
	// Compute scattering functions for surface interaction
	isect.ComputeScatteringFunctions(ray, arena);
	if (!isect.bsdf)
		return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
	Vector3f wo = isect.wo;
	// Compute emitted light if ray hit an area light source
	L += isect.Le(wo);
	if (scene.lights.size() > 0) {
		// Compute direct lighting and distance from surface to light
		Float lightDistance = 0.;
		if (strategy == LightStrategy::UniformSampleAll)
			L += UniformSampleAllLights(isect, scene, arena, sampler,
				nLightSamples, false, &lightDistance);
		else
			L += UniformSampleOneLight(isect, scene, arena, sampler, false, &lightDistance);
		distance += lightDistance;
	}
	return IntegrationResult(L, HistogramSample{ L.y(), distance });
}

DirectToFIntegrator *CreateDirectToFIntegrator(
	const ParamSet &params, std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera) {
	LightStrategy strategy;
	std::string st = params.FindOneString("strategy", "all");
	if (st == "one")
		strategy = LightStrategy::UniformSampleOne;
	else if (st == "all")
		strategy = LightStrategy::UniformSampleAll;
	else {
		Warning(
			"Strategy \"%s\" for direct lighting unknown. "
			"Using \"all\".",
			st.c_str());
		strategy = LightStrategy::UniformSampleAll;
	}
	return new DirectToFIntegrator(strategy, camera, sampler);
}
