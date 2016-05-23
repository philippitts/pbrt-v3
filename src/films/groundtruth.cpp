/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

#include "films/groundtruth.h"
#include "stats.h"

// GroundTruth Method Definitions
GroundTruthFilm::GroundTruthFilm(const Point2i &resolution, const Bounds2f &cropWindow,
	std::unique_ptr<Filter> filter, Float diagonal, const std::string &filename,
	Float scale) :
	Film(resolution, cropWindow, std::move(filter), diagonal, filename, scale) {
	pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
}

std::unique_ptr<FilmTile> GroundTruthFilm::GetFilmTile(const Bounds2i &sampleBounds) {
	// Bound image pixels that samples in _sampleBounds_ contribute to
	Vector2f halfPixel = Vector2f(0.5f, 0.5f);
	Bounds2f floatBounds = (Bounds2f)sampleBounds;
	Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
	Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
		Point2i(1, 1);
	Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
	return std::unique_ptr<GroundTruthFilmTile>(new GroundTruthFilmTile(
		tilePixelBounds, filter->radius, filterTable, filterTableWidth));
}

void GroundTruthFilm::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
	ProfilePhase p(Prof::MergeFilmTile);
	std::lock_guard<std::mutex> lock(mutex);

	GroundTruthFilmTile* gtTile = static_cast<GroundTruthFilmTile*>(tile.get());
	if (gtTile == nullptr) {
		Warning("Skipping alien film tile in MergeFilmTile");
		return;
	}

	for (Point2i pixel : gtTile->GetPixelBounds()) {
		// Merge _pixel_ into _Film::pixels_
		const GroundTruthTilePixel &tilePixel = gtTile->GetPixel(pixel);
		Pixel &mergePixel = GetPixel(pixel);

		mergePixel.value.L += tilePixel.value.L;
		mergePixel.value.distance += tilePixel.value.distance;
		mergePixel.filterWeightSum += tilePixel.filterWeightSum;
		mergePixel.nContribs += tilePixel.nContribs;
	}
}

void GroundTruthFilm::SetImage(const Spectrum *img) const {
	// TODO: Stub - maybe in the future ground truth images could be loaded?
	// It doesn't make sense to just load radiance from images.
}

void GroundTruthFilm::AddSplat(const Point2f &p, const IntegrationResult &v) {
	if (v.L.HasNaNs()) {
		Warning("Film ignoring splatted spectrum with NaN values");
		return;
	}
	if (v.histogramSamples.empty()) {
		Warning("Film ignoring splatted integration result with no samples");
		return;
	}
	if (v.histogramSamples.size() > 1) {
		Warning("Film ignoring extra values of splatted integration result");
	}

	ProfilePhase pp(Prof::SplatFilm);
	if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
	Pixel &pixel = GetPixel((Point2i)p);

	pixel.splatValue.distance += v.histogramSamples[0].distance;
	pixel.splatValue.L += v.histogramSamples[0].L;
}

void GroundTruthFilm::WriteImage(Float splatScale) {
	FILE* fp = fopen(filename.c_str(), "w");
	if (!fp) Severe("GroundTruthFilm file %s could not be opened", filename);

	for (Point2i p : croppedPixelBounds) {
		Pixel &pixel = GetPixel(p);

		Float invFilterWt = 0;
		if (pixel.filterWeightSum != 0) invFilterWt = (Float)1 / pixel.filterWeightSum;

		Float invNContribs = 0;
		if (pixel.nContribs != 0) invNContribs = (Float)1 / pixel.nContribs;

		Float L = pixel.value.L * invFilterWt;
		Float distance = pixel.value.distance * invNContribs;

		L += pixel.splatValue.L * splatScale;
		distance += pixel.splatValue.distance * splatScale;

		fprintf(fp, "# %d %d %f %f ", p.x, p.y, distance, L);
	}

	fclose(fp);
}

GroundTruthFilmTile::GroundTruthFilmTile(const Bounds2i &pixelBounds,
	const Vector2f &filterRadius, const Float *filterTable, int filterTableSize)
	: FilmTile(pixelBounds, filterRadius, filterTable, filterTableSize) {
	pixels = std::vector<GroundTruthTilePixel>(std::max(0, pixelBounds.Area()));
}

void GroundTruthFilmTile::AddSample(const Point2f &pFilm, const IntegrationResult &integration,
	Float sampleWeight) {
	if (integration.histogramSamples.empty()) {
		Warning("Film ignoring added integration result with no samples");
		return;
	}
	if (integration.histogramSamples.size() > 1) {
		Warning("Film ignoring extra values of added integration result");
	}

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

			// Update pixel value
			GroundTruthTilePixel &pixel = GetPixel(Point2i(x, y));
			pixel.value.distance += integration.histogramSamples[0].distance * sampleWeight;
			pixel.value.L += integration.histogramSamples[0].L * sampleWeight * filterWeight;
			pixel.filterWeightSum += filterWeight;
			pixel.nContribs++;
		}
	}
}

GroundTruthTilePixel& GroundTruthFilmTile::GetPixel(const Point2i &p) {
	Assert(InsideExclusive(p, pixelBounds));
	int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
	int offset =
		(p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
	return pixels[offset];
}

GroundTruthFilm *CreateGroundTruthFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
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

	return new GroundTruthFilm(Point2i(xres, yres), crop, std::move(filter), diagonal,
		filename, scale);
}