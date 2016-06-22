/*
This file is part of the Time-of-Flight Tracer program. It is not part of
the original PBRT source distribution. See the included license file for
more information.

Copyright(c) 2016 Microsoft Corporation

Author: Phil Pitts
*/

#include "stdafx.h"

#include "films/signal.h"
#include "stats.h"

// SignalFilm Method Definitions
SignalFilm::SignalFilm(const Point2i &resolution, const Bounds2f &cropWindow,
	std::unique_ptr<Filter> filter, Float diagonal, const std::string &filename,
	Float scale, std::vector<Float>& frequencies, std::vector<Float>& phases) :
	Film(resolution, cropWindow, std::move(filter), diagonal, filename, scale),
	frequencies(frequencies),
	phases(phases) {
	pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
	for (int i = 0; i < croppedPixelBounds.Area(); i++) {
		pixels[i].Initialize(frequencies.size(), phases.size());
	}
}

std::unique_ptr<FilmTile> SignalFilm::GetFilmTile(const Bounds2i &sampleBounds) {
	// Bound image pixels that samples in _sampleBounds_ contribute to
	Vector2f halfPixel = Vector2f(0.5f, 0.5f);
	Bounds2f floatBounds = (Bounds2f)sampleBounds;
	Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
	Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
		Point2i(1, 1);
	Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
	return std::unique_ptr<SignalFilmTile>(new SignalFilmTile(
		tilePixelBounds, filter->radius, filterTable, filterTableWidth,
		frequencies, phases));
}

void SignalFilm::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
	ProfilePhase p(Prof::MergeFilmTile);
	std::lock_guard<std::mutex> lock(mutex);

	SignalFilmTile *signalTile = static_cast<SignalFilmTile*>(tile.get());
	if (signalTile == nullptr) {
		Warning("Skipping alien film tile in MergeFilmTile");
		return;
	}

	for (Point2i pixel : signalTile->GetPixelBounds()) {
		// Merge _pixel_ into _Film::pixels_
		const SignalTilePixel &tilePixel = signalTile->GetPixel(pixel);
		Pixel &mergePixel = GetPixel(pixel);

		if (tilePixel.values.size() != mergePixel.values.size()) {
			Severe("SignalFilm value buffers are different sizes");
		}

		for (size_t i = 0; i < tilePixel.values.size(); ++i) {
			mergePixel.values[i] += tilePixel.values[i];
		}

		mergePixel.filterWeightSum += tilePixel.filterWeightSum;
	}
}

void SignalFilm::SetImage(const Spectrum *img) const {
	// TODO: Stub - maybe in the future histogram images could be loaded?
	// It doesn't make sense to just load radiance from images.
}

void SignalFilm::AddSplat(const Point2f &p, const IntegrationResult &v) {
	if (v.L.HasNaNs()) {
		Warning("Film ignoring splatted spectrum with NaN values");
		return;
	}
	ProfilePhase pp(Prof::SplatFilm);
	if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
	Pixel &pixel = GetPixel((Point2i)p);

	for (auto sample : v.histogramSamples) {
		for (size_t i = 0; i < frequencies.size(); ++i) {
			for (size_t j = 0; j < phases.size(); ++j) {
				float kernel = 
					GetKernel(frequencies[i], phases[j], sample.pathLength);
				size_t idx = i * phases.size() + j;
				pixel.splatValues[idx] += sample.L.y() * kernel;
			}
		}
	}
}

void SignalFilm::WriteImage(Float splatScale) {
	FILE* fp = fopen(filename.c_str(), "w");
	if (!fp) Severe("SignalFilm file %s could not be opened", filename);

	for (Point2i p : croppedPixelBounds) {
		Pixel &pixel = GetPixel(p);
		if (pixel.splatValues.size() != pixel.values.size()) {
			Severe("SignalFilm value buffers are different sizes");
		}

		Float invWt = 0.f;
		for (size_t i = 0; i < pixel.values.size(); ++i) {
			if (pixel.filterWeightSum != 0) {
				invWt = (Float)1 / pixel.filterWeightSum;
				pixel.values[i] *= invWt;
			}

			pixel.values[i] += splatScale * pixel.splatValues[i];
			pixel.values[i] *= scale;

			fprintf(fp, "%f ", pixel.values[i]);
		}
	}
	fclose(fp);
}

SignalFilmTile::SignalFilmTile(const Bounds2i &pixelBounds,
	const Vector2f &filterRadius, const Float *filterTable, int filterTableSize,
	std::vector<Float>& frequencies, std::vector<Float>& phases)
	: FilmTile(pixelBounds, filterRadius, filterTable, filterTableSize),
	frequencies(frequencies), phases(phases) {
	pixels = std::vector<SignalTilePixel>(std::max(0, pixelBounds.Area()));
	for (int i = 0; i < pixelBounds.Area(); i++) {
		pixels[i].Initialize(frequencies.size(), phases.size());
	}
}

void SignalFilmTile::AddSample(const Point2f &pFilm, 
	const IntegrationResult &integration, Float sampleWeight) {
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
			SignalTilePixel &pixel = GetPixel(Point2i(x, y));
			pixel.filterWeightSum += filterWeight;

			for (auto sample : integration.histogramSamples) {
				for (size_t i = 0; i < frequencies.size(); ++i) {
					for (size_t j = 0; j < phases.size(); ++j) {
						float kernel =
							GetKernel(frequencies[i], phases[j], sample.pathLength);
						size_t idx = i * phases.size() + j;
						pixel.values[idx] += sample.L.y() * kernel;
					}
				}
			}
		}
	}
}

SignalTilePixel& SignalFilmTile::GetPixel(const Point2i &p) {
	Assert(InsideExclusive(p, pixelBounds));
	int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
	int offset =
		(p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
	return pixels[offset];
}

SignalFilm *CreateSignalFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
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

	int freqi;
	std::vector<Float> frequencies;
	const Float *freq = params.FindFloat("frequencies", &freqi);
	if (freq && freqi != 0) frequencies = std::vector<Float>(freq, freq + freqi);

	int pi;
	std::vector<Float> phases;
	const Float *p = params.FindFloat("phases", &pi);
	if (p && pi != 0) phases = std::vector<Float>(p, p + pi);

	return new SignalFilm(Point2i(xres, yres), crop, std::move(filter), diagonal,
		filename, scale, frequencies, phases);
}