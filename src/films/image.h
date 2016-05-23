
/*
	pbrt source code is Copyright(c) 1998-2015
						Matt Pharr, Greg Humphreys, and Wenzel Jakob.

	This file is part of pbrt.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are
	met:

	- Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.

	- Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
	IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
	TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
	PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_FILMS_IMAGE_H
#define PBRT_FILMS_IMAGE_H
#include "stdafx.h"

 // films/image.h*
#include "pbrt.h"
#include "film.h"

// ImageFilmTilePixel Declarations
struct ImageTilePixel {
	Spectrum contribSum = 0.f;
	Float filterWeightSum = 0.f;
};

// ImageFilm Declarations
class ImageFilm : public Film {
public:
	// ImageFilm Public Methods
	ImageFilm(const Point2i &resolution, const Bounds2f &cropWindow,
		std::unique_ptr<Filter> filter, Float diagonal,
		const std::string &filename, Float scale);
	std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);
	void MergeFilmTile(std::unique_ptr<FilmTile> tile);
	void SetImage(const Spectrum *img) const;
	void AddSplat(const Point2f &p, const IntegrationResult &v);
	void WriteImage(Float splatScale);

private:
	// ImageFilm Private Data
	struct Pixel {
		Pixel() { xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0; }
		Float xyz[3];
		Float filterWeightSum;
		AtomicFloat splatXYZ[3];
		Float pad;
	};
	std::unique_ptr<Pixel[]> pixels;

	// ImageFilm Private Methods
	Pixel &GetPixel(const Point2i &p) {
		Assert(InsideExclusive(p, croppedPixelBounds));
		int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
		int offset = (p.x - croppedPixelBounds.pMin.x) +
			(p.y - croppedPixelBounds.pMin.y) * width;
		return pixels[offset];
	}
};

class ImageFilmTile : public FilmTile {
public:
	// ImageFilmTile Public Methods
	ImageFilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
		const Float *filterTable, int filterTableSize);
	void AddSample(const Point2f &pFilm, const IntegrationResult &integration,
		Float sampleWeight);
	ImageTilePixel &GetPixel(const Point2i &p);

private:
	// FilmTile Private Data
	std::vector<ImageTilePixel> pixels;
};

ImageFilm *CreateImageFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

#endif  // PBRT_FILMS_IMAGE_H
