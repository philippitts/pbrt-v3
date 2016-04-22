
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
	IntegrationResult() { samples = nullptr; nSamples = 0; }
	IntegrationResult(const IntegrationResult& src);
	IntegrationResult(Spectrum &L);
	IntegrationResult(Spectrum &L, std::queue<HistogramSample> &samples);
	~IntegrationResult();

	Spectrum L;
	size_t nSamples;
	HistogramSample* samples;
};


#endif  // PBRT_CORE_INTEGRATIONRESULT_H