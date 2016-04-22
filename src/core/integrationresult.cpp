
/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

// core/integrationresult.cpp*
#include "integrationresult.h"

IntegrationResult::IntegrationResult(const IntegrationResult& src) {
	L = src.L;
	nSamples = src.nSamples;
	samples = new HistogramSample[nSamples];

	for (size_t i = 0; i < nSamples; i++) {
		samples[i] = src.samples[i];
	}
}

IntegrationResult::IntegrationResult(Spectrum &L) {
	this->L = L;
	this->nSamples = 0;
	this->samples = nullptr;
}

IntegrationResult::IntegrationResult(Spectrum &L, std::queue<HistogramSample> &samples) {
	this->L = L;
	this->nSamples = samples.size();
	this->samples = new HistogramSample[nSamples];

	for (size_t i = 0; i < nSamples; i++) {
		this->samples[i] = samples.front();
		samples.pop();
	}
}

IntegrationResult::~IntegrationResult() {
	if (samples != nullptr) {
		delete samples;
		samples = nullptr;
		nSamples = 0;
	}
}