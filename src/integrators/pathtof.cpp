
/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

#include "integrators/pathtof.h"
#include "scene.h"
#include "interaction.h"
#include "paramset.h"
#include "bssrdf.h"
#include "stats.h"
#include "histogram.h"

#include <queue>

STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);

// ToFPathIntegrator Method Definitions
IntegrationResult PathToFIntegrator::Li(const RayDifferential &r, const Scene &scene,
                            Sampler &sampler, MemoryArena &arena,
                            int depth) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum directL(0.f), indirectL(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
	Float pathLength = 0.f;

	std::queue<HistogramSample> histogram;
	
    int bounces;
    for (bounces = 0;; ++bounces) {
        // Find next path vertex and accumulate contribution

		Spectrum localIndirectL(0.f);

        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

        // Possibly add emitted light at intersection
        if (bounces == 0 || specularBounce) {
            // Add emitted light at path vertex or from the environment
            if (foundIntersection)
				directL += beta * isect.Le(-ray.d);
            else
                for (const auto &light : scene.lights)
					directL += beta * light->Le(ray);
        }

        // Terminate path if ray escaped or _maxDepth_ was reached
        if (!foundIntersection || bounces >= maxDepth) break;

		// Update the length of the path traced so far
		pathLength += (ray.o - isect.p).Length();

        // Compute scattering functions and skip over medium boundaries
        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
            ray = isect.SpawnRay(ray.d);
            bounces--;
            continue;
        }

		Float lightPathLength;

        // Sample illumination from lights to find path contribution.
        // (But skip this for perfectly specular BSDFs.)
        if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >
            0) {
            ++totalPaths;
            Spectrum Ld =
                beta * UniformSampleOneLight(isect, scene, arena, sampler, false, &lightPathLength);
            if (Ld.IsBlack()) ++zeroRadiancePaths;
            Assert(Ld.y() >= 0.f);
			localIndirectL += Ld;
        }

        // Sample BSDF to get new path direction
        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) break;
        beta *= f * AbsDot(wi, isect.shading.n) / pdf;
        Assert(beta.y() >= 0.f);
        Assert(std::isinf(beta.y()) == false);
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        ray = isect.SpawnRay(wi);

        // Subsurface scattering is currently unsupported
		Assert(!(isect.bssrdf && (flags & BSDF_TRANSMISSION)));

		// Possibly terminate the path with Russian roulette
		if (bounces > 3) {
			Float q = std::max((Float).05, 1 - beta.y());
			if (sampler.Get1D() < q) break;
			beta /= 1 - q;
			Assert(std::isinf(beta.y()) == false);
		}

		// Add the histogram sample to the queue
		HistogramSample histSample;
		histSample.pathLength = pathLength + lightPathLength;
		histSample.L = Spectrum(directL + localIndirectL).y();
		histogram.push(HistogramSample(histSample));

		indirectL += localIndirectL;
    }
    ReportValue(pathLength, bounces);
    return IntegrationResult(Spectrum(directL + indirectL), histogram);
}

PathToFIntegrator *CreatePathToFIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new PathToFIntegrator(maxDepth, camera, sampler);
}
