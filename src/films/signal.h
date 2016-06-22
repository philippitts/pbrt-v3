
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

#ifndef PBRT_FILMS_SIGNAL_H
#define PBRT_FILMS_SIGNAL_H

#define SPEED_LIGHT 299792458
#define M_PI 3.14159265358979323846f

#include "stdafx.h"

// films/signal.h*
#include "pbrt.h"
#include "film.h"
#include "paramset.h"

// SignalTilePixel Declarations
class SignalTilePixel {
public:
	// SignalTilePixel Public Methods
	SignalTilePixel() { filterWeightSum = 0.f; }
	void Initialize(size_t nFrequencies, size_t nPhases) {
		values = std::vector<Float>(nFrequencies * nPhases);
	}

	// SignalTilePixel Public Members
	std::vector<Float> values;
	Float filterWeightSum;
};

// SignalFilmTile Declarations
class SignalFilmTile : public FilmTile {
public:
	// SignalFilmTile Public Methods
	SignalFilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
		const Float *filterTable, int filterTableSize, 
		std::vector<Float>& frequencies, std::vector<Float>& phases);
	void AddSample(const Point2f &pFilm, const IntegrationResult &integration,
		Float sampleWeight);
	SignalTilePixel &GetPixel(const Point2i &p);

private:
	// SignalFilmTile Private Methods
	std::vector<SignalTilePixel> pixels;
	std::vector<Float> frequencies;
	std::vector<Float> phases;
};

// SignalFilm Declarations
class SignalFilm : public Film {
public:
	// SignalFilm Public Methods
	SignalFilm(const Point2i &resolution, const Bounds2f &cropWindow,
		std::unique_ptr<Filter> filter, Float diagonal,
		const std::string &filename, Float scale, 
		std::vector<Float>& frequencies, std::vector<Float>& phases);

	std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);
	void MergeFilmTile(std::unique_ptr<FilmTile> tile);
	void SetImage(const Spectrum *img) const;
	void AddSplat(const Point2f &p, const IntegrationResult &v);
	void WriteImage(Float splatScale);

private:
	// Film Private Data
	class Pixel {
	public:
		Pixel() { filterWeightSum = 0; }
		void Initialize(size_t nFrequencies, size_t nPhases) {
			values = std::vector<Float>(nFrequencies * nPhases);
			splatValues = std::vector<Float>(nFrequencies * nPhases);
		}

		std::vector<Float> values;
		std::vector<Float> splatValues;
		Float filterWeightSum;
	};
	std::unique_ptr<Pixel[]> pixels;
	std::vector<Float> frequencies;
	std::vector<Float> phases;

	// Film Private Methods
	Pixel &GetPixel(const Point2i &p) {
		Assert(InsideExclusive(p, croppedPixelBounds));
		int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
		int offset = (p.x - croppedPixelBounds.pMin.x) +
			(p.y - croppedPixelBounds.pMin.y) * width;
		return pixels[offset];
	}
};

inline float GetKernel(float frequency, float phase, float pathLength) {
	return cos(4 * M_PI * frequency * pathLength / SPEED_LIGHT + phase);
}

SignalFilm *CreateSignalFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

#endif  // PBRT_FILMS_SIGNAL_H
