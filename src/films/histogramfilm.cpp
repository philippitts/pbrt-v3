/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

#include "films/histogramfilm.h"
#include "stats.h"

// HistogramFilm Method Definitions
HistogramFilm::HistogramFilm(const Point2i &resolution, const Bounds2f &cropWindow,
	std::unique_ptr<Filter> filter, Float diagonal, const std::string &filename,
	Float scale, Float binSize, Float maxPathLength, Float minL) : 
	Film(resolution, cropWindow, std::move(filter), diagonal, filename, scale),
	binSize(binSize),
	maxPathLength(maxPathLength),
	minL(minL) {
	pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
	for (int i = 0; i < croppedPixelBounds.Area(); i++) {
		pixels[i].Initialize(binSize, maxPathLength);
	}
}

std::unique_ptr<FilmTile> HistogramFilm::GetFilmTile(const Bounds2i &sampleBounds) {
	// Bound image pixels that samples in _sampleBounds_ contribute to
	Vector2f halfPixel = Vector2f(0.5f, 0.5f);
	Bounds2f floatBounds = (Bounds2f)sampleBounds;
	Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
	Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
		Point2i(1, 1);
	Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
	return std::unique_ptr<HistogramFilmTile>(new HistogramFilmTile(
		tilePixelBounds, filter->radius, filterTable, filterTableWidth,
		binSize, maxPathLength));
}

void HistogramFilm::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
	ProfilePhase p(Prof::MergeFilmTile);
	std::lock_guard<std::mutex> lock(mutex);

	HistogramFilmTile *histogramTile = static_cast<HistogramFilmTile*>(tile.get());
	if (histogramTile == nullptr) {
		Warning("Skipping alien film tile in MergeFilmTile");
		return;
	}

	for (Point2i pixel : histogramTile->GetPixelBounds()) {
		// Merge _pixel_ into _Film::pixels_
		const HistogramTilePixel &tilePixel = histogramTile->GetPixel(pixel);
		Pixel &mergePixel = GetPixel(pixel);

		if (tilePixel.histogram.binSize != mergePixel.histogram.binSize ||
			tilePixel.histogram.bins.size() != mergePixel.histogram.bins.size()) {
			Severe("HistogramFilm histograms have different sizes");
		}

		for (size_t i = 0; i < tilePixel.histogram.bins.size(); ++i) {
			mergePixel.histogram.bins[i] += tilePixel.histogram.bins[i];
		}
		mergePixel.filterWeightSum += tilePixel.filterWeightSum;
	}
}

void HistogramFilm::SetImage(const Spectrum *img) const {
	// TODO: Stub - maybe in the future histogram images could be loaded?
	// It doesn't make sense to just load radiance from images.
}

void HistogramFilm::AddSplat(const Point2f &p, const IntegrationResult &v) {
	if (v.L.HasNaNs()) {
		Warning("Film ignoring splatted spectrum with NaN values");
		return;
	}
	ProfilePhase pp(Prof::SplatFilm);
	if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
	Pixel &pixel = GetPixel((Point2i)p);
	
	size_t nBins = pixel.splatHistogram.bins.size();
	for (auto sample : v.histogramSamples) {
		size_t binIdx = (size_t)(sample.pathLength / binSize);
		if (binIdx < nBins) {
			pixel.splatHistogram.bins[binIdx] += sample.L;
		}
	}
}

void HistogramFilm::WriteImage(Float splatScale) {
	FILE* fp = fopen(filename.c_str(), "w");
	if (!fp) Severe("HistogramFilm file %s could not be opened", filename);

	for (Point2i p : croppedPixelBounds) {
		Pixel &pixel = GetPixel(p);
		if (binSize != pixel.histogram.binSize ||
			pixel.histogram.binSize != pixel.splatHistogram.binSize ||
			pixel.histogram.bins.size() != pixel.splatHistogram.bins.size()) {
			Severe("HistogramFilm histograms have different sizes");
		}

		Float invFilterWt = 0;
		size_t nBins = pixel.histogram.bins.size();
		if (pixel.filterWeightSum != 0) invFilterWt = (Float)1 / pixel.filterWeightSum;

		bool isFirst = true;
		for (size_t i = 0; i < nBins; ++i) {
			Histogram& histogram = pixel.histogram;

			Float rgb[3];
			histogram.bins[i].ToRGB(rgb);
			if (pixel.filterWeightSum != 0) {
				Float invWt = (Float)1 / pixel.filterWeightSum;
				rgb[0] = std::max((Float)0, rgb[0] * invWt);
				rgb[1] = std::max((Float)0, rgb[1] * invWt);
				rgb[2] = std::max((Float)0, rgb[2] * invWt);
			}

			Float splatRGB[3];
			pixel.splatHistogram.bins[i].ToRGB(splatRGB);

			rgb[0] += splatScale * splatRGB[0];
			rgb[1] += splatScale * splatRGB[1];
			rgb[2] += splatScale * splatRGB[2];

			rgb[0] *= scale;
			rgb[1] *= scale;
			rgb[2] *= scale;

			float L = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];

			if (L >= minL) {
				if (isFirst) {
					fprintf(fp, "# %d %d ", p.x, p.y);
					isFirst = false;
				}
				fprintf(fp, "%f %f ", pixel.histogram.binSize * i, L);
			}
		}
	}

	fclose(fp);
}

HistogramFilmTile::HistogramFilmTile(const Bounds2i &pixelBounds, 
	const Vector2f &filterRadius, const Float *filterTable, int filterTableSize, 
	Float binSize, Float maxPathLength)
	: FilmTile(pixelBounds, filterRadius, filterTable, filterTableSize) {
	pixels = std::vector<HistogramTilePixel>(std::max(0, pixelBounds.Area()));
	for (int i = 0; i < pixelBounds.Area(); ++i) {
		pixels[i].Initialize(binSize, maxPathLength);
	}
}

void HistogramFilmTile::AddSample(const Point2f &pFilm, const IntegrationResult &integration,
	Float sampleWeight) {
	// Compute sample's raster bounds
	Point2f pFilmDiscrete = pFilm - Vector2f(0.5f, 0.5f);
	Point2i p0 = (Point2i)Ceil(pFilmDiscrete - filterRadius);
	Point2i p1 =
		(Point2i)Floor(pFilmDiscrete + filterRadius) + Point2i(1, 1);
	p0 = Max(p0, pixelBounds.pMin);
	p1 = Min(p1, pixelBounds.pMax);

	// Loop over filter support and add sample to pixel arrays

	// Precompute $x$ and $y$ filter table offsets
	int *ifx = ALLOCA(int, p1.x - p0.x);
	for (int x = p0.x; x < p1.x; ++x) {
		Float fx = std::abs((x - pFilmDiscrete.x) * invFilterRadius.x *
			filterTableSize);
		ifx[x - p0.x] = std::min((int)std::floor(fx), filterTableSize - 1);
	}
	int *ify = ALLOCA(int, p1.y - p0.y);
	for (int y = p0.y; y < p1.y; ++y) {
		Float fy = std::abs((y - pFilmDiscrete.y) * invFilterRadius.y *
			filterTableSize);
		ify[y - p0.y] = std::min((int)std::floor(fy), filterTableSize - 1);
	}
	for (int y = p0.y; y < p1.y; ++y) {
		for (int x = p0.x; x < p1.x; ++x) {
			// Evaluate filter value at $(x,y)$ pixel
			int offset = ify[y - p0.y] * filterTableSize + ifx[x - p0.x];
			Float filterWeight = filterTable[offset];

			// Update pixel histogram
			HistogramTilePixel &pixel = GetPixel(Point2i(x, y));
			pixel.filterWeightSum += filterWeight;

			int nBins = pixel.histogram.bins.size();
			for (auto sample : integration.histogramSamples) {
				int binIdx = (int)(sample.pathLength / pixel.histogram.binSize);
				if (binIdx < nBins) {
					pixel.histogram.bins[binIdx] += sample.L * sampleWeight * filterWeight;
				}
			}
		}
	}
}

HistogramTilePixel& HistogramFilmTile::GetPixel(const Point2i &p) {
	Assert(InsideExclusive(p, pixelBounds));
	int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
	int offset =
		(p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
	return pixels[offset];
}

HistogramFilm *CreateHistogramFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
	// Intentionally use FindOneString() rather than FindOneFilename() here
	// so that the rendered image is left in the working directory, rather
	// than the directory the scene file lives in.
	std::string filename = params.FindOneString("filename", "");
	if (PbrtOptions.imageFile != "") {
		if (filename != "") {
			Warning(
				"Output filename supplied on command line, \"%s\", ignored "
				"due to filename provided in scene description file, \"%s\".",
				PbrtOptions.imageFile.c_str(), filename.c_str());
		}
		else
			filename = PbrtOptions.imageFile;
	}
	if (filename == "") filename = "pbrt.dat";

	int xres = params.FindOneInt("xresolution", 1280);
	int yres = params.FindOneInt("yresolution", 720);
	if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
	if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
	Bounds2f crop(Point2f(0, 0), Point2f(1, 1));
	int cwi;
	const Float *cr = params.FindFloat("cropwindow", &cwi);
	if (cr && cwi == 4) {
		crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
		crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
		crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
		crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
	}
	else if (cr)
		Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);

	Float scale = params.FindOneFloat("scale", 1.);
	Float diagonal = params.FindOneFloat("diagonal", 35.);
	Float binSize = params.FindOneFloat("binsize", 0.1);
	Float maxPathLength = params.FindOneFloat("maxpathlength", 10.);
	Float minL = params.FindOneFloat("minL", 0.0001);

	return new HistogramFilm(Point2i(xres, yres), crop, std::move(filter), diagonal,
		filename, scale, binSize, maxPathLength, minL);
}