
/*
    pbrt source code is Copyright(c) 1998-2016
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

// core/film.h*
#include "pbrt.h"
#include "api.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "stats.h"
#include "parallel.h"
#include "imageio.h"
#include <sys/stat.h>

//////////////////////
// PrISE-3D Updates //
//////////////////////
#include "tools/childprocess.hpp"
//////////////////////////
// End PrISE-3D Updates //
//////////////////////////

namespace pbrt {


class FilmTilePixelMoN {

    public:

        unsigned k; // number of means clusters
        unsigned index; // keep track of index used
        
        std::vector<Float> rvalues; // store sum of r lightness
        std::vector<Float> gvalues; // store sum of g lightness
        std::vector<Float> bvalues; // store sum of b lightness

        std::vector<Float> weights; // store sum of lightness
        std::vector<unsigned> counters; // number of elements

        FilmTilePixelMoN(unsigned _k = 9){
            k = _k;
            index = 0;

            rvalues = std::vector<Float>(k);
            gvalues = std::vector<Float>(k);
            bvalues = std::vector<Float>(k);

            weights = std::vector<Float>(k);

            counters = std::vector<unsigned>(k);

            for (int i = 0; i < k; i++){
                rvalues[i] = Float(0.);
                gvalues[i] = Float(0.);
                bvalues[i] = Float(0.);
                weights[i] = Float(0.);
                counters[i] = 0;
            }
        }

        void add(Spectrum sample, Float weight){
            
            Float rgb[3];
            sample.ToRGB(rgb);

            rvalues.at(index) += rgb[0];
            gvalues.at(index) += rgb[1];
            bvalues.at(index) += rgb[2];
            weights.at(index) += weight;
            counters.at(index) += 1;
            
            index += 1;

            if (index >= k)
                index = 0;
        }
};

// FilmTilePixel Declarations
struct FilmTilePixel {
    Spectrum contribSum = 0.f;
    Float filterWeightSum = 0.f;
};

// Film Declarations
class Film {
  public:
    // Film Public Methods
    Film(const Point2i &resolution, const Bounds2f &cropWindow,
         std::unique_ptr<Filter> filter, Float diagonal,
         const std::string &filename, Float scale,
         Float maxSampleLuminance = Infinity);
    Bounds2i GetSampleBounds() const;
    Bounds2f GetPhysicalExtent() const;
    std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);
    void MergeFilmTile(std::unique_ptr<FilmTile> tile);

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    void ApplyDL();
    Float getMaxZBuffer();
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////
    
    void SetImage(const Spectrum *img) const;
    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    void LoadRawlsImage(const std::string filename) const;
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

    void AddSplat(const Point2f &p, Spectrum v);
    void WriteImage(Float splatScale = 1);

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    void WriteImageTemp(int index, Float splatScale = 1);
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

    void Clear();

    // Film Public Data
    const Point2i fullResolution;
    const Float diagonal;
    std::unique_ptr<Filter> filter;
    const std::string filename;
    Bounds2i croppedPixelBounds;

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    std::unique_ptr<ChildProcess> child_process; // python model interaction
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    // Film public Methods
    void UpdateAdditionals(const Point2i &p, const Ray &ray) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        
        zbuffer[offset] = ray.tMax;
        normals[offset] = ray.nn;
    }
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

  private:
    
    struct PixelMoN {

        PixelMoN(unsigned _k = 9) { 
            k = _k;
            index = 0;

            rvalues = std::vector<Float>(k);
            gvalues = std::vector<Float>(k);
            bvalues = std::vector<Float>(k);

            weights = std::vector<Float>(k);

            counters = std::vector<unsigned>(k);

            for (int i = 0; i < k; i++){
                rvalues[i] = Float(0.);
                gvalues[i] = Float(0.);
                bvalues[i] = Float(0.);
                weights[i] = Float(0.);
                counters[i] = 0;
            }
        }

        AtomicFloat splatXYZ[3];
        Float filterWeightSum;

        unsigned k; // number of means clusters
        unsigned index; // keep track of index used
        
        std::vector<Float> rvalues; // store sum of r lightness
        std::vector<Float> gvalues; // store sum of g lightness
        std::vector<Float> bvalues; // store sum of b lightness

        std::vector<Float> weights; // store sum of lightness
        std::vector<unsigned> counters; // number of elements
        Float xyz[3];

        void append(FilmTilePixelMoN tile) {

            for (int i = 0; i < k; i++){
                rvalues[i] += tile.rvalues[i];
                gvalues[i] += tile.gvalues[i];
                bvalues[i] += tile.bvalues[i];

                weights[i] += tile.weights[i];

                counters[i] += tile.counters[i];                
            }
        }

        void estimateRGB() {

            Float rgb[3];

            rgb[0] = estimate(rvalues);
            rgb[1] = estimate(rvalues);
            rgb[2] = estimate(rvalues);

            Float weight = estimate(weights);

            rgb[0] /= weight;
            rgb[1] /= weight;
            rgb[2] /= weight;

            Float currentXYZ[3];

            RGBToXYZ(rgb, currentXYZ);

            xyz[0] = currentXYZ[0];
            xyz[1] = currentXYZ[1];
            xyz[2] = currentXYZ[2];
        }

        void estimateWeight(Float weight) const{
            
            weight = estimate(weights);
        }

        Float estimate(std::vector<Float> cvalues) const{

            std::vector<Float> means;

            unsigned nElements = cvalues.size();

            for (unsigned i = 0; i < nElements; i++){
                means.push_back(cvalues[i] / counters[i]);
            }

            std::sort(means.begin(), means.end()); 

            // compute median from means
            if (nElements % 2 == 1){
                return means[int(nElements/2)];
            }
            else{
                int k_mean = int(nElements/2);
                return (means[k_mean - 1] + means[k_mean]) / 2;
            }
        }
    };

    // Film Private Data
    struct Pixel {
        Pixel() { xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0; }
        Float xyz[3];
        Float filterWeightSum;
        AtomicFloat splatXYZ[3];
        Float pad;
    };
    std::unique_ptr<PixelMoN[]> pixels;
    static PBRT_CONSTEXPR int filterTableWidth = 32; // TODO check if really 32 ?
    Float filterTable[filterTableWidth * filterTableWidth];
    std::mutex mutex;
    const Float scale;
    const Float maxSampleLuminance;

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    std::unique_ptr<Float[]> zbuffer;
    std::unique_ptr<Normal3f[]> normals;
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

    // Film Private Methods
    // Pixel &GetPixel(const Point2i &p) {
    //     CHECK(InsideExclusive(p, croppedPixelBounds));
    //     int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
    //     int offset = (p.x - croppedPixelBounds.pMin.x) +
    //                  (p.y - croppedPixelBounds.pMin.y) * width;
    //     return pixels[offset];
    // }

    PixelMoN &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    Float &GetBufferPoint(const Point2i &p) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        return zbuffer[offset];
    }

    Normal3f &GetNormalPoint(const Point2i &p) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        return normals[offset];
    }
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////
};

class FilmTile {
  public:

    // FilmTile Public Methods
    FilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
             const Float *filterTable, int filterTableSize,
             Float maxSampleLuminance)
        : pixelBounds(pixelBounds),
          filterRadius(filterRadius),
          invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
          filterTable(filterTable),
          filterTableSize(filterTableSize),
          maxSampleLuminance(maxSampleLuminance) {
        pixels = std::vector<FilmTilePixelMoN>(std::max(0, pixelBounds.Area()));
    }

    void AddSample(const Point2f &pFilm, Spectrum L,
                   Float sampleWeight = 1.) {
        ProfilePhase _(Prof::AddFilmSample);
        if (L.y() > maxSampleLuminance)
            L *= maxSampleLuminance / L.y();
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

                // Update pixel values with filtered sample contribution
                // FilmTilePixel &pixel = GetPixel(Point2i(x, y));
                // pixel.contribSum += L * sampleWeight * filterWeight;
                // pixel.filterWeightSum += filterWeight;
                FilmTilePixelMoN &pixel = GetPixel(Point2i(x, y));
                pixel.add(L * sampleWeight * filterWeight, filterWeight); // TODO : check if use of weight

            }
        }
    }
    FilmTilePixelMoN &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }
    const FilmTilePixelMoN &GetPixel(const Point2i &p) const {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    Bounds2i GetPixelBounds() const { return pixelBounds; }

  private:
    // FilmTile Private Data
    const Bounds2i pixelBounds;
    const Vector2f filterRadius, invFilterRadius;
    const Float *filterTable;
    const int filterTableSize;
    std::vector<FilmTilePixelMoN> pixels;
    const Float maxSampleLuminance;
    friend class Film;
};

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

}  // namespace pbrt

#endif  // PBRT_CORE_FILM_H
