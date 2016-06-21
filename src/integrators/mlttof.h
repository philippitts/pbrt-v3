
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

#ifndef PBRT_INTEGRATORS_MLTTOF_H
#define PBRT_INTEGRATORS_MLTTOF_H
#include "stdafx.h"

// integrators/mlttof.h*
#include "pbrt.h"
#include "integrator.h"
#include "sampler.h"
#include "spectrum.h"
#include "integrators/mlt.h"

// MLT ToF Declarations
class MLTToFIntegrator : public Integrator {
public:
	// MLTToFIntegrator Public Methods
	MLTToFIntegrator(std::shared_ptr<const Camera> camera, int maxDepth,
		int nBootstrap, int nChains, int mutationsPerPixel,
		Float sigma, Float largeStepProbability)
		: camera(camera),
		maxDepth(maxDepth),
		nBootstrap(nBootstrap),
		nChains(nChains),
		mutationsPerPixel(mutationsPerPixel),
		sigma(sigma),
		largeStepProbability(largeStepProbability) {}
	void Render(const Scene &scene);
	HistogramSample Sample(const Scene &scene, MemoryArena &arena,
		const std::unique_ptr<Distribution1D> &lightDistr,
		MLTSampler &sampler, int k, Point2f *pRaster);

private:
	// MLTToFIntegrator Private Data
	std::shared_ptr<const Camera> camera;
	const int maxDepth;
	const int nBootstrap;
	const int nChains;
	const int mutationsPerPixel;
	const Float sigma, largeStepProbability;
};

MLTToFIntegrator *CreateMLTToFIntegrator(const ParamSet &params,
	std::shared_ptr<const Camera> camera);

#endif  // PBRT_INTEGRATORS_MLTTOF_H
