
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

#ifndef PBRT_CORE_FILM_H
#define PBRT_CORE_FILM_H
#include "stdafx.h"

// core/film.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "integrationresult.h"

// Film Declarations
class Film {
  public:
    // Film Public Methods
    Film(const Point2i &resolution, const Bounds2f &cropWindow,
         std::unique_ptr<Filter> filter, Float diagonal,
         const std::string &filename, Float scale);
    Bounds2i GetSampleBounds() const;
    Bounds2f GetPhysicalExtent() const;
	
	virtual std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds) = 0;
	virtual void MergeFilmTile(std::unique_ptr<FilmTile> tile) = 0;
	virtual void SetImage(const Spectrum *img) const = 0;
	virtual void AddSplat(const Point2f &p, const IntegrationResult &v) = 0;
	virtual void WriteImage(Float splatScale = 1.) = 0;

    // Film Public Data
    const Point2i fullResolution;
    const Float diagonal;
    std::unique_ptr<Filter> filter;
    const std::string filename;
    Bounds2i croppedPixelBounds;

protected:
	// Film Protected Data
	static PBRT_CONSTEXPR int filterTableWidth = 16;
	Float filterTable[filterTableWidth * filterTableWidth];
	std::mutex mutex;
	const Float scale;
};

// FilmTile Declarations
class FilmTile {
public:
	// FilmTile Public Methods
	FilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
		const Float *filterTable, int filterTableSize) :
		pixelBounds(pixelBounds),
		filterRadius(filterRadius),
		invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
		filterTable(filterTable),
		filterTableSize(filterTableSize) { }
	Bounds2i GetPixelBounds() const { return pixelBounds; }
	
	virtual void AddSample(const Point2f &pFilm, const IntegrationResult &integration,
		Float sampleWeight = 1.) = 0;

protected:
	// FilmTile Protected Data
	const Bounds2i pixelBounds;
	const Vector2f filterRadius, invFilterRadius;
	const Float *filterTable;
	const int filterTableSize;
};

#endif  // PBRT_CORE_FILM_H
