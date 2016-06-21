
/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

// integrators/mlttof.cpp*
#include "integrators/mlttof.h"
#include "integrators/bdpttof.h"
#include "scene.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "camera.h"
#include "stats.h"
#include "filters/box.h"
#include "paramset.h"
#include "sampling.h"
#include "progressreporter.h"
#include "integrationresult.h"

STAT_TIMER("Time/Rendering", renderingTime);
STAT_PERCENT("Integrator/Acceptance rate", acceptedMutations, totalMutations);

// MLTSampler Constants
static const int cameraStreamIndex = 0;
static const int lightStreamIndex = 1;
static const int connectionStreamIndex = 2;
static const int nSampleStreams = 3;

// MLT ToF Method Definitions
HistogramSample MLTToFIntegrator::Sample(const Scene &scene, MemoryArena &arena,
	const std::unique_ptr<Distribution1D> &lightDistr,
	MLTSampler &sampler, int depth, Point2f *pRaster) {
	sampler.StartStream(cameraStreamIndex);
	// Determine the number of available strategies and pick a specific one
	int s, t, nStrategies;
	if (depth == 0) {
		nStrategies = 1;
		s = 0;
		t = 2;
	}
	else {
		nStrategies = depth + 2;
		s = std::min((int)(sampler.Get1D() * nStrategies), nStrategies - 1);
		t = nStrategies - s;
	}

	// Generate a camera subpath with exactly _t_ vertices
	Vertex *cameraVertices = arena.Alloc<Vertex>(t);
	Bounds2f sampleBounds = (Bounds2f)camera->film->GetSampleBounds();
	*pRaster = sampleBounds.Lerp(sampler.Get2D());
	if (GenerateCameraSubpath(scene, sampler, arena, t, *camera, *pRaster,
		cameraVertices) != t)
		return HistogramSample();

	// Generate a light subpath with exactly _s_ vertices
	sampler.StartStream(lightStreamIndex);
	Vertex *lightVertices = arena.Alloc<Vertex>(s);
	if (GenerateLightSubpath(scene, sampler, arena, s, cameraVertices[0].time(),
		*lightDistr, lightVertices) != s)
		return HistogramSample();

	// Execute connection strategy and return the radiance estimate
	sampler.StartStream(connectionStreamIndex);
	HistogramSample sample = ConnectBDPTToF(scene, lightVertices,
		cameraVertices, s, t, *lightDistr, *camera, sampler, pRaster);
	sample.L *= nStrategies;
	return sample;
}

void MLTToFIntegrator::Render(const Scene &scene) {
	ProfilePhase p(Prof::IntegratorRender);
	std::unique_ptr<Distribution1D> lightDistr =
		ComputeLightPowerDistribution(scene);
	// Generate bootstrap samples and compute normalization constant $b$
	int nBootstrapSamples = nBootstrap * (maxDepth + 1);
	std::vector<Float> bootstrapWeights(nBootstrapSamples, 0);
	if (scene.lights.size() > 0) {
		ProgressReporter progress(nBootstrap / 256,
			"Generating bootstrap paths");
		std::vector<MemoryArena> bootstrapThreadArenas(MaxThreadIndex());
		int chunkSize = Clamp(nBootstrap / 128, 1, 8192);
		ParallelFor([&](int i) {
			// Generate _i_th bootstrap sample
			MemoryArena &arena = bootstrapThreadArenas[ThreadIndex];
			for (int depth = 0; depth <= maxDepth; ++depth) {
				int rngIndex = i * (maxDepth + 1) + depth;
				MLTSampler sampler(mutationsPerPixel, rngIndex, sigma,
					largeStepProbability, nSampleStreams);
				Point2f pRaster;
				bootstrapWeights[rngIndex] =
					Sample(scene, arena, lightDistr, sampler, depth, &pRaster).L.y();
				arena.Reset();
			}
			if ((i + 1 % 256) == 0) progress.Update();
		}, nBootstrap, chunkSize);
		progress.Done();
	}
	Distribution1D bootstrap(&bootstrapWeights[0], nBootstrapSamples);
	Float b = bootstrap.funcInt * (maxDepth + 1);

	// Run _nChains_ Markov chains in parallel
	Film &film = *camera->film;
	int64_t nTotalMutations =
		(int64_t)mutationsPerPixel * (int64_t)film.GetSampleBounds().Area();
	if (scene.lights.size() > 0) {
		StatTimer timer(&renderingTime);
		const int progressFrequency = 32768;
		ProgressReporter progress(nTotalMutations / progressFrequency,
			"Rendering");
		ParallelFor([&](int i) {
			int64_t nChainMutations =
				std::min((i + 1) * nTotalMutations / nChains, nTotalMutations) -
				i * nTotalMutations / nChains;
			// Follow {i}th Markov chain for _nChainMutations_
			MemoryArena arena;

			// Select initial state from the set of bootstrap samples
			RNG rng(i);
			int bootstrapIndex = bootstrap.SampleDiscrete(rng.UniformFloat());
			int depth = bootstrapIndex % (maxDepth + 1);

			// Initialize local variables for selected state
			MLTSampler sampler(mutationsPerPixel, bootstrapIndex, sigma,
				largeStepProbability, nSampleStreams);
			Point2f pCurrent;
			HistogramSample current =
				Sample(scene, arena, lightDistr, sampler, depth, &pCurrent);

			// Run the Markov chain for _nChainMutations_ steps
			for (int64_t j = 0; j < nChainMutations; ++j) {
				sampler.StartIteration();
				Point2f pProposed;
				HistogramSample proposed =
					Sample(scene, arena, lightDistr, sampler, depth, &pProposed);
				// Compute acceptance probability for proposed sample
				Float accept = std::min((Float)1, proposed.L.y() / current.L.y());

				// Splat both current and proposed samples to _film_
				Spectrum storage;
				if (accept > 0) {
					storage = proposed.L;
					proposed.L *= accept / proposed.L.y();
					film.AddSplat(pProposed, 
						IntegrationResult(proposed.L, proposed));
					proposed.L = storage;
				}
				storage = current.L;
				current.L *= (1 - accept) / current.L.y();
				film.AddSplat(pCurrent,
					IntegrationResult(current.L, current));
				current.L = storage;

				// Accept or reject the proposal
				if (rng.UniformFloat() < accept) {
					pCurrent = pProposed;
					current = proposed;
					sampler.Accept();
					++acceptedMutations;
				}
				else
					sampler.Reject();
				++totalMutations;
				if ((i * nTotalMutations / nChains + j) % progressFrequency ==
					0)
					progress.Update();
				arena.Reset();
			}
		}, nChains);
		progress.Done();
	}

	// Store final image computed with MLT
	camera->film->WriteImage(b / mutationsPerPixel);
}

MLTToFIntegrator *CreateMLTToFIntegrator(const ParamSet &params,
	std::shared_ptr<const Camera> camera) {
	int maxDepth = params.FindOneInt("maxdepth", 5);
	int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
	int64_t nChains = params.FindOneInt("chains", 1000);
	int mutationsPerPixel = params.FindOneInt("mutationsperpixel", 100);
	Float largeStepProbability =
		params.FindOneFloat("largestepprobability", 0.3f);
	Float sigma = params.FindOneFloat("sigma", .01f);
	if (PbrtOptions.quickRender) {
		mutationsPerPixel = std::max(1, mutationsPerPixel / 16);
		nBootstrap = std::max(1, nBootstrap / 16);
	}
	return new MLTToFIntegrator(camera, maxDepth, nBootstrap, nChains,
		mutationsPerPixel, sigma, largeStepProbability);
}