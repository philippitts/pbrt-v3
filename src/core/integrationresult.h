
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

#ifndef PBRT_CORE_INTEGRATIONRESULT_H
#define PBRT_CORE_INTEGRATIONRESULT_H
#include "stdafx.h"

// core/integrationresult.h*
#include "pbrt.h"
#include "histogram.h"

#include <queue>

class IntegrationResult {
public:
	IntegrationResult() {}
	IntegrationResult(const IntegrationResult& src) :
		L(src.L),
		histogramSamples(src.histogramSamples) {}
	IntegrationResult(Spectrum &L) : L(L) {}
	IntegrationResult(Spectrum &L, std::vector<HistogramSample> &samples)
		: L(L), histogramSamples(samples) {}
	IntegrationResult(Spectrum &L, std::queue<HistogramSample> &samples) : L(L) {
		histogramSamples.resize(samples.size());
		for (size_t i = 0; !samples.empty(); i++) {
			histogramSamples[i] = samples.front();
			samples.pop();
		}
	}
	IntegrationResult(Spectrum &L, HistogramSample& sample) : L(L) {
		histogramSamples = { sample };
	}

	Spectrum L;
	std::vector<HistogramSample> histogramSamples;
};


#endif  // PBRT_CORE_INTEGRATIONRESULT_H