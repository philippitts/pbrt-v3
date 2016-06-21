
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

#ifndef PBRT_INTEGRATORS_BDPTTOF_H
#define PBRT_INTEGRATORS_BDPTTOF_H
#include "stdafx.h"

// integrators/bdpttof.h*
#include "pbrt.h"
#include "integrators/bdpt.h"
#include "histogram.h"

// BDPT ToF Declarations
class BDPTToFIntegrator : public Integrator {
public:
	// BDPTToFIntegrator Public Methods
	BDPTToFIntegrator(std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera, int maxDepth,
		bool visualizeStrategies, bool visualizeWeights)
		: sampler(sampler),
		camera(camera),
		maxDepth(maxDepth),
		visualizeStrategies(visualizeStrategies),
		visualizeWeights(visualizeWeights) {}
	void Render(const Scene &scene);

private:
	// BDPTToFIntegrator Private Data
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<const Camera> camera;
	const int maxDepth;
	const bool visualizeStrategies;
	const bool visualizeWeights;
};

Spectrum ConnectBDPTToF(const Scene &scene, Vertex *lightVertices,
	Vertex *cameraVertices, int s, int t,
	const Distribution1D &lightDistr, const Camera &camera,
	Sampler &sampler, Point2f *pRaster, HistogramSample& sample,
	Float *misWeight = nullptr);

Float PathLength(Vertex* vertices, int lastIndex);

BDPTToFIntegrator *CreateBDPTToFIntegrator(const ParamSet &params,
	std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera);

#endif  // PBRT_INTEGRATORS_BDPTTOF_H
