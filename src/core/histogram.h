
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

class HistogramSample {
public:
	HistogramSample() {
		distance = 0.f;
	}

	HistogramSample(Float distance, Spectrum radiance) { 
		this->distance = distance;
		this->radiance = radiance; 
	}

	Float distance;
	Spectrum radiance;
};

#endif  // PBRT_CORE_HISTOGRAM_H
