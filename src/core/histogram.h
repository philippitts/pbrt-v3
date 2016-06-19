
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

#ifndef PBRT_CORE_HISTOGRAM_H
#define PBRT_CORE_HISTOGRAM_H
#include "stdafx.h"

// core/histogram.h*
#include "pbrt.h"
#include "spectrum.h"

class Histogram {
public:
	// Histogram Public Methods
	Histogram() : binSize(0) { }
	Histogram(Float binSize, Float maxDistance) : binSize(binSize) {
		if (binSize <= 0) Severe("Illegal histogram bin size");
		bins.resize(maxDistance / binSize);
	}

	// Histogram Public Methods
	Float binSize;
	std::vector<Spectrum> bins;

};

struct HistogramSample {
	HistogramSample() { distance = 0.f; }
	HistogramSample(Spectrum L, Float distance) : L(L), distance(distance) { }

	Spectrum L;
	Float distance;
};

#endif  // PBRT_CORE_HISTOGRAM_H
