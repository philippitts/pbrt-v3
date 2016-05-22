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

#ifndef PBRT_FILMS_GROUNDTRUTH_H
#define PBRT_FILMS_GROUNDTRUTH_H
#include "stdafx.h"

// films/groundtruth.h*
#include "pbrt.h"
#include "film.h"
#include "histogram.h"
#include "paramset.h"

// GroundTruthTilePixel Declarations
struct GroundTruthTilePixel {
	HistogramSample value;
	Float filterWeightSum;
};

// GroundTruthFilmTile Declarations
class GroundTruthFilmTile : public FilmTile {
public:
	// GroundTruthFilmTile Public Methods
	GroundTruthFilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
		const Float *filterTable, int filterTableSize);
	void AddSample(const Point2f &pFilm, const IntegrationResult &integration,
		Float sampleWeight = 1.);
	GroundTruthTilePixel &GetPixel(const Point2i &p);
	Bounds2i GetPixelBounds() const { return pixelBounds; }

private:
	// GroundTruthFilmTile Private Methods
	std::vector<GroundTruthTilePixel> pixels;
};

// GroundTruthFilm Declarations
class GroundTruthFilm : public Film {
public:
	// GroundTruthFilm Public Methods
	GroundTruthFilm(const Point2i &resolution, const Bounds2f &cropWindow,
		std::unique_ptr<Filter> filter, Float diagonal,
		const std::string &filename, Float scale);

	std::unique_ptr<GroundTruthFilmTile> GetFilmTile(const Bounds2i &sampleBounds);
	void MergeFilmTile(std::unique_ptr<GroundTruthFilmTile> tile);
	void SetImage(const Spectrum *img) const;
	void AddSplat(const Point2f &p, const IntegrationResult &v);
	void WriteImage(Float splatScale = 1);

private:
	// Film Private Data
	class Pixel {
	public:
		Pixel() { filterWeightSum = 0; }

		HistogramSample value;
		HistogramSample splatValue;
		Float filterWeightSum;
	};
	std::unique_ptr<Pixel[]> pixels;

	// Film Private Methods
	Pixel &GetPixel(const Point2i &p) {
		Assert(InsideExclusive(p, croppedPixelBounds));
		int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
		int offset = (p.x - croppedPixelBounds.pMin.x) +
			(p.y - croppedPixelBounds.pMin.y) * width;
		return pixels[offset];
	}
};

GroundTruthFilm *CreateGroundTruthFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

#endif  // PBRT_FILMS_GROUNDTRUTH_H
