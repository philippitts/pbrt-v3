/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

// integrators/bdpttof.cpp*
#include "integrators/bdpttof.h"
#include "films/image.h"
#include "sampler.h"
#include "integrator.h"
#include "stats.h"
#include "filters/box.h"
#include "paramset.h"
#include "progressreporter.h"

STAT_TIMER("Time/Rendering", renderingTime);
STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
STAT_INT_DISTRIBUTION("Integrator/Intersections", bounces);
STAT_FLOAT_DISTRIBUTION("Integrator/Path length", pathLength);

// BDPT ToF Method Definitions
void BDPTToFIntegrator::Render(const Scene &scene) {
	ProfilePhase p(Prof::IntegratorRender);
	// Compute _lightDistr_ for sampling lights proportional to power
	std::unique_ptr<Distribution1D> lightDistr =
		ComputeLightPowerDistribution(scene);

	// Partition the image into tiles
	Film *film = camera->film;
	const Bounds2i sampleBounds = film->GetSampleBounds();
	const Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	const int nXTiles = (sampleExtent.x + tileSize - 1) / tileSize;
	const int nYTiles = (sampleExtent.y + tileSize - 1) / tileSize;
	ProgressReporter reporter(nXTiles * nYTiles, "Rendering");

	// Allocate buffers for debug visualization
	const int bufferCount = (1 + maxDepth) * (6 + maxDepth) / 2;
	std::vector<std::unique_ptr<Film>> weightFilms(bufferCount);
	if (visualizeStrategies || visualizeWeights) {
		for (int depth = 0; depth <= maxDepth; ++depth) {
			for (int s = 0; s <= depth + 2; ++s) {
				int t = depth + 2 - s;
				if (t == 0 || (s == 1 && t == 1)) continue;

				char filename[32];
				snprintf(filename, sizeof(filename),
					"bdpt_d%02i_s%02i_t%02i.exr", depth, s, t);

				weightFilms[BufferIndex(s, t)] = std::unique_ptr<Film>(new ImageFilm(
					film->fullResolution,
					Bounds2f(Point2f(0, 0), Point2f(1, 1)),
					std::unique_ptr<Filter>(CreateBoxFilter(ParamSet())),
					film->diagonal * 1000, filename, 1.f));
			}
		}
	}

	// Render and write the output image to disk
	if (scene.lights.size() > 0) {
		StatTimer timer(&renderingTime);
		ParallelFor2D([&](const Point2i tile) {
			// Render a single tile using BDPT
			MemoryArena arena;
			int seed = tile.y * nXTiles + tile.x;
			std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);
			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
			int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
			int y0 = sampleBounds.pMin.y + tile.y * tileSize;
			int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
			Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
			std::unique_ptr<FilmTile> filmTile =
				camera->film->GetFilmTile(tileBounds);
			for (Point2i pPixel : tileBounds) {
				tileSampler->StartPixel(pPixel);
				do {
					// Generate a single sample using BDPT
					Point2f pFilm = (Point2f)pPixel + tileSampler->Get2D();

					// Trace the camera and light subpaths
					Vertex *cameraVertices = arena.Alloc<Vertex>(maxDepth + 2);
					Vertex *lightVertices = arena.Alloc<Vertex>(maxDepth + 1);
					int nCamera = GenerateCameraSubpath(
						scene, *tileSampler, arena, maxDepth + 2, *camera,
						pFilm, cameraVertices);
					int nLight = GenerateLightSubpath(
						scene, *tileSampler, arena, maxDepth + 1,
						cameraVertices[0].time(), *lightDistr, lightVertices);

					std::vector<HistogramSample> samples((nCamera + 1) * (nLight + 1));

					// Execute all BDPT connection strategies
					Spectrum L(0.f);
					HistogramSample sample;
					for (int t = 1; t <= nCamera; ++t) {
						for (int s = 0; s <= nLight; ++s) {
							int depth = t + s - 2;
							if ((s == 1 && t == 1) || depth < 0 ||
								depth > maxDepth)
								continue;
							// Execute the $(s, t)$ connection strategy and
							// update _L_
							Point2f pFilmNew = pFilm;
							Float misWeight = 0.f;
							Spectrum Lpath = ConnectBDPTToF(
								scene, lightVertices, cameraVertices, s, t,
								*lightDistr, *camera, *tileSampler, &pFilmNew,
								sample, &misWeight);
							if (visualizeStrategies || visualizeWeights) {
								Spectrum value;
								if (visualizeStrategies)
									value =
									misWeight == 0 ? 0 : Lpath / misWeight;
								if (visualizeWeights) value = Lpath;
								weightFilms[BufferIndex(s, t)]->AddSplat(
									pFilmNew, value);
							}
							if (t != 1) {
								L += Lpath;
								samples[t * (nLight + 1) + s] = sample;
							}
							else {
								film->AddSplat(pFilmNew, 
									IntegrationResult(Lpath, sample));
							}
						}
					}
					filmTile->AddSample(pFilm, IntegrationResult(L, samples));
					arena.Reset();
				} while (tileSampler->StartNextSample());
			}
			film->MergeFilmTile(std::move(filmTile));
			reporter.Update();
		}, Point2i(nXTiles, nYTiles));
		reporter.Done();
	}
	film->WriteImage(1.0f / sampler->samplesPerPixel);

	// Write buffers for debug visualization
	if (visualizeStrategies || visualizeWeights) {
		const Float invSampleCount = 1.0f / sampler->samplesPerPixel;
		for (size_t i = 0; i < weightFilms.size(); ++i)
			if (weightFilms[i]) weightFilms[i]->WriteImage(invSampleCount);
	}
}

Spectrum ConnectBDPTToF(const Scene &scene, Vertex *lightVertices,
	Vertex *cameraVertices, int s, int t,
	const Distribution1D &lightDistr, const Camera &camera,
	Sampler &sampler, Point2f *pRaster, HistogramSample& sample,
	Float *misWeightPtr) {
	Spectrum L(0.f);
	// Ignore invalid connections related to infinite area lights
	if (t > 1 && s != 0 && cameraVertices[t - 1].type == VertexType::Light)
		return Spectrum(0.f);

	// Perform connection and write contribution to _L_
	Vertex sampled;
	if (s == 0) {
		// Interpret the camera subpath as a complete path
		const Vertex &pt = cameraVertices[t - 1];
		if (pt.IsLight()) L = pt.Le(scene, cameraVertices[t - 2]) * pt.beta;
		sample.pathLength = PathLength(cameraVertices, t - 1);
	}
	else if (t == 1) {
		// Sample a point on the camera and connect it to the light subpath
		const Vertex &qs = lightVertices[s - 1];
		if (qs.IsConnectible()) {
			VisibilityTester vis;
			Vector3f wi;
			Float pdf;
			Spectrum Wi = camera.Sample_Wi(qs.GetInteraction(), sampler.Get2D(),
				&wi, &pdf, pRaster, &vis);
			if (pdf > 0 && !Wi.IsBlack()) {
				// Initialize dynamically sampled vertex and _L_ for $t=1$ case
				sampled = Vertex::CreateCamera(&camera, vis.P1(), Wi / pdf);
				L = qs.beta * qs.f(sampled) * vis.Tr(scene, sampler) *
					sampled.beta;
				if (qs.IsOnSurface()) L *= AbsDot(wi, qs.ns());
				sample.pathLength = PathLength(lightVertices, s - 1);
				sample.pathLength += Distance(qs.p(), sampled.p());
			}
		}
	}
	else if (s == 1) {
		// Sample a point on a light and connect it to the camera subpath
		const Vertex &pt = cameraVertices[t - 1];
		if (pt.IsConnectible()) {
			Float lightPdf;
			VisibilityTester vis;
			Vector3f wi;
			Float pdf;
			int lightNum =
				lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
			const std::shared_ptr<Light> &light = scene.lights[lightNum];
			Spectrum lightWeight = light->Sample_Li(
				pt.GetInteraction(), sampler.Get2D(), &wi, &pdf, &vis);
			if (pdf > 0 && !lightWeight.IsBlack()) {
				EndpointInteraction ei(vis.P1(), light.get());
				sampled =
					Vertex::CreateLight(ei, lightWeight / (pdf * lightPdf), 0);
				sampled.pdfFwd = sampled.PdfLightOrigin(scene, pt, lightDistr);
				L = pt.beta * pt.f(sampled) * sampled.beta;
				if (pt.IsOnSurface()) L *= AbsDot(wi, pt.ns());
				// Only check visibility if the path would carry radiance.
				if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
				sample.pathLength = PathLength(cameraVertices, t - 1);
				sample.pathLength += Distance(pt.p(), sampled.p());
			}
		}
	}
	else {
		// Handle all other bidirectional connection cases
		const Vertex &qs = lightVertices[s - 1], &pt = cameraVertices[t - 1];
		if (qs.IsConnectible() && pt.IsConnectible()) {
			L = qs.beta * qs.f(pt) * pt.f(qs) * pt.beta;
			if (!L.IsBlack()) L *= G(scene, sampler, qs, pt);
			sample.pathLength = Distance(qs.p(), pt.p());
			sample.pathLength += PathLength(cameraVertices, t - 1);
			sample.pathLength += PathLength(lightVertices, s - 1);
		}
	}

	++totalPaths;
	if (L.IsBlack()) ++zeroRadiancePaths;
	ReportValue(bounces, s + t - 2);

	// Compute MIS weight for connection strategy
	Float misWeight =
		L.IsBlack() ? 0.f : MISWeight(scene, lightVertices, cameraVertices,
			sampled, s, t, lightDistr);
	L *= misWeight;
	if (misWeightPtr) *misWeightPtr = misWeight;
	sample.L = L;
	return L;
}

Float PathLength(Vertex* vertices, int lastIndex) {
	Float result = 0.f;
	for (int i = 0; i < lastIndex; ++i) {
		result += Distance(vertices[i].p(), vertices[i + 1].p());
	}
	return result;
}

BDPTToFIntegrator *CreateBDPTToFIntegrator(const ParamSet &params,
	std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera) {
	int maxDepth = params.FindOneInt("maxdepth", 5);
	bool visualizeStrategies = params.FindOneBool("visualizestrategies", false);
	bool visualizeWeights = params.FindOneBool("visualizeweights", false);

	if ((visualizeStrategies || visualizeWeights) && maxDepth > 5) {
		Warning(
			"visualizestrategies/visualizeweights was enabled, limiting "
			"maxdepth to 5");
		maxDepth = 5;
	}

	return new BDPTToFIntegrator(sampler, camera, maxDepth, visualizeStrategies,
		visualizeWeights);
}
