
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


// core/film.cpp*
#include "film.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"
#include <fstream>
#include <unistd.h>
#include <bits/stdc++.h>
#include "progressreporter.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);

// Film Method Definitions
Film::Film(const Point2i &resolution, const Bounds2f &cropWindow,
           std::unique_ptr<Filter> filt, Float diagonal,
           const std::string &filename, Float scale, Float maxSampleLuminance)
    : fullResolution(resolution),
      diagonal(diagonal * .001),
      filter(std::move(filt)),
      filename(filename),
      scale(scale),
      maxSampleLuminance(maxSampleLuminance) {
    // Compute film image bounds
    croppedPixelBounds =
        Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                         std::ceil(fullResolution.y * cropWindow.pMin.y)),
                 Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                         std::ceil(fullResolution.y * cropWindow.pMax.y)));
    LOG(INFO) << "Created film with full resolution " << resolution <<
        ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
        croppedPixelBounds;

    // Allocate film image storage
    pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);

    //////////////////////
    // PrISE-3D Updates //
    //////////////////////
    // zbuffer = std::unique_ptr<Float[]>(new Float[croppedPixelBounds.Area()]);
    // normals = std::unique_ptr<Normal3f[]>(new Normal3f[croppedPixelBounds.Area()]); 

    // for (int i = 0; i < croppedPixelBounds.Area(); ++i){
    //     zbuffer[i] = 0; // default value
    // }

    // // initialize child process
    // if (PbrtOptions.nn_path.length() > 0){
    //     char *const argv[] = {
    //         const_cast<char*>("python"),
    //         const_cast<char*>(PbrtOptions.nn_path.c_str()),
    //         NULL
    //     };

    //     child_process = std::unique_ptr<ChildProcess>(
    //                 new ChildProcess(
    //                     std::string("python"),
    //                     argv
    //                     )
    //                 );
    // }
    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////

    filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

    // Precompute filter weight table
    int offset = 0;
    for (int y = 0; y < filterTableWidth; ++y) {
        for (int x = 0; x < filterTableWidth; ++x, ++offset) {
            Point2f p;
            p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
            p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
            filterTable[offset] = filter->Evaluate(p);
        }
    }
}

//////////////////////
// PrISE-3D Updates //
//////////////////////
void Film::LoadRawlsImage(const std::string filename) const {

    // only one read buffer used for the whole function
    std::ifstream rf(filename, std::ios::out | std::ios::binary);

    if(!rf) {
      std::cout << "Cannot open file!" << std::endl;
    }

    std::string line; 
    unsigned nbChanels, width, height;
    char c; // to get space of new line char

    // READ IHDR info
    bool ihdrBegin = false;

    while (!ihdrBegin && std::getline(rf, line)) { 

        if (line.find(std::string("IHDR")) != std::string::npos){
            ihdrBegin = true;
            std::getline(rf, line); // avoid data size line

            rf.read((char *) &width, sizeof(unsigned));
            rf.get(c);
            rf.read((char *) &height, sizeof(unsigned));
            rf.get(c);
            rf.read((char *) &nbChanels, sizeof(unsigned));
            rf.get(c);
        }
    }

    bool dataBegin = false;

    // READ DATA info
    // case of data chunck begin
    while (!dataBegin && std::getline(rf, line)) { 

        if (line.find(std::string("DATA")) != std::string::npos){
            dataBegin = true;
        }
    }

    // resize buffer if necessary
    float* buffer = new float[height * width * nbChanels];

    std::getline(rf, line);
    unsigned size = std::stoi(line);
    float chanelValue;

    for(unsigned y = 0; y < height; y++){

        for(unsigned x = 0; x < width; x++) {

            for(unsigned j = 0; j < nbChanels; j++){
                rf.read((char *) &chanelValue, sizeof(float));  
                
                buffer[nbChanels * width * y + nbChanels * x + j] = chanelValue; 
            } 

            // go to next line
            rf.get(c);
        }
    }

    rf.close();

    if(!rf.good()) {
        std::cout << "Error occurred at writing time!" << std::endl;
    }

    // for each values into buffer set current film buffer
    int nPixels = croppedPixelBounds.Area();
    int offset = 0;

    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        p.xyz[0] = buffer[offset];
        p.xyz[1] = buffer[offset + 1];
        p.xyz[2] = buffer[offset + 2];

        offset++;
    }
}
//////////////////////////
// End PrISE-3D Updates //
//////////////////////////

Bounds2i Film::GetSampleBounds() const {
    Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
                               Vector2f(0.5f, 0.5f) - filter->radius),
                         Ceil(Point2f(croppedPixelBounds.pMax) -
                              Vector2f(0.5f, 0.5f) + filter->radius));
    return (Bounds2i)floatBounds;
}

Bounds2f Film::GetPhysicalExtent() const {
    Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);

    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);

    return std::unique_ptr<FilmTile>(new FilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void Film::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int c = 0; c < 3; ++c)
            pixel.splatXYZ[c] = pixel.xyz[c] = 0;
        pixel.filterWeightSum = 0;
    }
}

void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->pixelBounds;
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) mergePixel.xyz[i] += xyz[i];
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

//////////////////////
// PrISE-3D Updates //
//////////////////////
// Float Film::getMaxZBuffer(){

//     // TODO : To improve...
//     Float maxZBuffer = 0;

//     for (int i = 0; i < croppedPixelBounds.Area(); i++){
//         if (zbuffer[i] > maxZBuffer && zbuffer[i] != Infinity){
//             maxZBuffer = Float(zbuffer[i]);
//         }
//     }

//     return maxZBuffer;
// };


// void Film::ApplyDL() {

//     VLOG(1) << "Apply DL and merge film on tiles using " << PbrtOptions.DLConfidence << "%";
//     std::lock_guard<std::mutex> lock(mutex);

//     // build model input
//     // expected (1, 7, 32, 32)
//     // 1. mean pixel values [preprocessing: log + normalization ((x_i - min_Log(X))/(max_Log(X) - min_Log(X)))]
//     // 2. normals values [preprocessing: normalization ((x_i + 1) * 0.5)]
//     // 3. zbuffer values [preprocessing: normalization (x_i / max(X)) ]

//     // check tile size
//     int nChannels = 7;
//     int hSize = 32;
//     int wSize = 32;
//     int nChannelValues = hSize * wSize;
//     unsigned nValues = nChannels * nChannelValues;

//     // store image width
//     unsigned imageWidth = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
//     unsigned imageHeight = croppedPixelBounds.pMax.y - croppedPixelBounds.pMin.y;
//     unsigned totalNChannelValues = imageWidth * imageHeight;

//     // the whole image buffer
//     std::unique_ptr<Float[]> rgbLogBased(new Float[3 * croppedPixelBounds.Area()]);

//     // get Max Z buffer information data
//     Float maxZBuffer = getMaxZBuffer();
//     Float splatScale = 1.; // by default

//     unsigned nWeightTiles = ceil(fullResolution.x / (float)wSize); 
//     unsigned nHeightTiles = ceil(fullResolution.y / (float)hSize);

//     float minLogRGB = 0.;  
//     float maxLogRGB = std::numeric_limits<double>::min(); // max log value of whole film image

//     // convert each pixel of film into rgb log based value
//     int offset = 0;

   
//     //ParallelFor2D([&](Point2i p) { 
//     for (Point2i p : croppedPixelBounds) {
//         // Convert pixel XYZ color to RGB
//         Pixel &pixel = GetPixel(p);

//         // int offset = (p.x - croppedPixelBounds.pMin.x) + (p.y - croppedPixelBounds.pMin.y) * imageWidth;

//         Float rgb[3];
//         XYZToRGB(pixel.xyz, rgb);

//         // Normalize pixel with weight sum
//         Float filterWeightSum = pixel.filterWeightSum;
//         if (filterWeightSum != 0) {
//             Float invWt = (Float)1 / filterWeightSum;
//             rgbLogBased[offset] = std::max((Float)0, rgb[0] * invWt);
//             rgbLogBased[offset + 1 * totalNChannelValues] = std::max((Float)0, rgb[1] * invWt);
//             rgbLogBased[offset + 2 * totalNChannelValues] = std::max((Float)0, rgb[2] * invWt);
//         }

//         // Add splat value at pixel
//         Float splatRGB[3];
//         Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
//                             pixel.splatXYZ[2]};
//         XYZToRGB(splatXYZ, splatRGB);

//         rgbLogBased[offset] += splatScale * splatRGB[0];
//         rgbLogBased[offset + 1 * totalNChannelValues] += splatScale * splatRGB[1];
//         rgbLogBased[offset + 2 * totalNChannelValues] += splatScale * splatRGB[2];

//         // Scale pixel value by _scale_
//         rgbLogBased[offset] *= scale;
//         rgbLogBased[offset + 1 * totalNChannelValues] *= scale;
//         rgbLogBased[offset + 2 * totalNChannelValues] *= scale;

//         // apply log on it (avoid negative values and -inf issue)
//         // TODO : use of log(x_i + 1) or adding epsilon to x_i ?
//         for (int i = 0; i < 3; i++){
//             rgbLogBased[offset + i * totalNChannelValues] = (Float)(log10(rgbLogBased[offset + i * totalNChannelValues] + 1));

//             // get max log based value
//             if (rgbLogBased[offset + i * totalNChannelValues] > maxLogRGB)
//                 maxLogRGB = rgbLogBased[offset + i * totalNChannelValues];
//         }

//         ++offset;
//     }//, croppedPixelBounds);

//     // update reporter base on tile denoising process
//     ProgressReporter denoising(nHeightTiles * nWeightTiles, "NN Denoising");
//     {   
//         // for each tile found
//         for (int j = 0; j < nHeightTiles; j++){
//             for (int i = 0; i < nWeightTiles; i++){

//                 unsigned wStart = i * wSize;
//                 unsigned hStart = j * hSize;

//                 unsigned finalWStart = wStart;
//                 unsigned finalHStart = hStart;

//                 unsigned wEnd = wStart + wSize;
//                 unsigned hEnd = hStart + hSize;

//                 // check and avoid out of bounds
//                 if (hEnd >= fullResolution.y && wEnd >= fullResolution.x){

//                     finalWStart = fullResolution.x - wSize;
//                     finalHStart = fullResolution.y - hSize;
//                 } else if (hEnd >= fullResolution.y){

//                     finalWStart = wStart;
//                     finalHStart = fullResolution.y - hSize;
//                 }else if (wEnd >= fullResolution.x){
                    
//                     finalWStart = fullResolution.x - wSize;
//                     finalHStart = hStart;
//                 }
                
//                 // final valid tile end for each axis
//                 wEnd = finalWStart + wSize;
//                 hEnd = finalHStart + hSize;

//                 // construct tile bounds based on this
//                 Bounds2i tileBounds(Point2i(finalWStart, finalHStart), Point2i(wEnd, hEnd));

//                 // prepare input values
//                 std::vector<float> inputValues(nValues);

//                 int pixelCounter = 0;

//                 // extract input model data
//                 for (Point2i pixel : tileBounds) {

//                     // Merge _pixel_ into _Film::pixels_
//                     Normal3f &normal = GetNormalPoint(pixel);
//                     Float &bufferPoint = GetBufferPoint(pixel);

//                     int offset = (pixel.x - croppedPixelBounds.pMin.x) + (pixel.y - croppedPixelBounds.pMin.y) * imageWidth;
                    
//                     // update and norm log based values
//                     float logDiff = maxLogRGB - minLogRGB;
//                     inputValues.at(pixelCounter) = (float)(rgbLogBased[offset] - minLogRGB) / logDiff;
//                     inputValues.at(pixelCounter + 1 * nChannelValues) = (float)(rgbLogBased[offset + 1 * totalNChannelValues] - minLogRGB) / logDiff;
//                     inputValues.at(pixelCounter + 2 * nChannelValues) = (float)(rgbLogBased[offset + 2 * totalNChannelValues] - minLogRGB) / logDiff;
                    
//                     // fill normal values and normalize
//                     inputValues.at(pixelCounter + 3 * nChannelValues) = (float)((normal.x + 1) * 0.5);
//                     inputValues.at(pixelCounter + 4 * nChannelValues) = (float)((normal.y + 1) * 0.5);
//                     inputValues.at(pixelCounter + 5 * nChannelValues) = (float)((normal.z + 1) * 0.5);

//                     // fill zBuffer value
//                     // normalize zBuffer value
//                     if (bufferPoint != Infinity)
//                         inputValues.at(pixelCounter + 6 * nChannelValues) = (float)(bufferPoint / maxZBuffer);
//                     else
//                         inputValues.at(pixelCounter + 6 * nChannelValues) = 1.;

//                     pixelCounter++;
//                 }

//                 // Send data to model !
//                 child_process->write_n_float32(&inputValues[0], nValues);

//                 unsigned status = 0;
//                 unsigned nOutputValues = nChannelValues * 3;

//                 std::vector<float> output_array(nOutputValues);

//                 // Read
//                 int code = child_process->read_n_float32(&output_array[0], nOutputValues);
//                 if (code) {
//                     std::cerr << "film.cpp: Error when reading float array" << std::endl;
//                 }

//                 // Check magic characters
//                 char c0 = child_process->read_char();
//                 char c1 = child_process->read_char();
//                 if (c0 == 'x' && c1 == '\n') {
//                     status = 0;
//                 } else {
//                     std::cerr << "film.cpp: magic characters don't match: ["<< c0 <<"] ["<< c1 <<"]" << std::endl;
//                     status = 1;
//                 }


//                 pixelCounter = 0;

//                 // denormalize (normalization * maxLog + exp) output model data
//                 // extract input model data
//                 for (Point2i pixel : tileBounds) {

//                     // Merge _pixel_ into _Film::pixels_
//                     Pixel &mergePixel = GetPixel(pixel);

//                     Float rgbValues[3];

//                     Float splatRGB[3];
//                     Float splatXYZ[3] = {mergePixel.splatXYZ[0], 
//                                     mergePixel.splatXYZ[1], 
//                                     mergePixel.splatXYZ[2]};
//                     XYZToRGB(splatXYZ, splatRGB);

//                     for (int i = 0; i < 3; ++i){
                        
//                         // reverse normalization
//                         // y = x_i * (max(X) - min(X)) + min(X)
//                         rgbValues[i] = output_array.at(pixelCounter + i * nChannelValues) * (maxLogRGB - minLogRGB) + minLogRGB;
                        
//                         // reverse log10 function (and remove 1 which was added to log function previously)
//                         rgbValues[i] = pow(10, rgbValues[i]) - 1;

//                         // reverse scale
//                         // Scale pixel value by _scale_
//                         rgbValues[i] /= scale;

//                         // reverse splat scale
//                         rgbValues[i] -= splatScale * splatRGB[i];

//                         // rescale based on sum weight contribution
//                         rgbValues[i] *= mergePixel.filterWeightSum;
//                     }

//                     // convert to XYZ
//                     Float xyzValues[3];
//                     RGBToXYZ(rgbValues, xyzValues);

//                     // merge final values using weighted sum
//                     // affects model data with percent of confidence
//                     for (int i = 0; i < 3; ++i){
//                         mergePixel.xyz[i] =  mergePixel.xyz[i] * (1. - PbrtOptions.DLConfidence) + xyzValues[i] * PbrtOptions.DLConfidence;
//                     }
                    
//                     pixelCounter++;
//                 }
            
//             denoising.Update();
            
//             }
//         }
//     }

//     denoising.Done();
// }
//////////////////////////
// End PrISE-3D Updates //
//////////////////////////

void Film::SetImage(const Spectrum *img) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        img[i].ToXYZ(p.xyz);
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}

void Film::AddSplat(const Point2f &p, Spectrum v) {
    ProfilePhase pp(Prof::SplatFilm);

    if (v.HasNaNs()) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with NaN values "
                                   "at (%f, %f)", p.x, p.y);
        return;
    } else if (v.y() < 0.) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with negative "
                                   "luminance %f at (%f, %f)", v.y(), p.x, p.y);
        return;
    } else if (std::isinf(v.y())) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with infinite "
                                   "luminance at (%f, %f)", p.x, p.y);
        return;
    }

    Point2i pi = Point2i(Floor(p));
    if (!InsideExclusive(pi, croppedPixelBounds)) return;
    if (v.y() > maxSampleLuminance)
        v *= maxSampleLuminance / v.y();
    Float xyz[3];
    v.ToXYZ(xyz);
    Pixel &pixel = GetPixel(pi);
    for (int i = 0; i < 3; ++i) pixel.splatXYZ[i].Add(xyz[i]);
}

void Film::WriteImage(Float splatScale) {
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[3 * offset]);

        // Normalize pixel with weight sum
        Float filterWeightSum = pixel.filterWeightSum;
        if (filterWeightSum != 0) {
            Float invWt = (Float)1 / filterWeightSum;
            rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
            rgb[3 * offset + 1] =
                std::max((Float)0, rgb[3 * offset + 1] * invWt);
            rgb[3 * offset + 2] =
                std::max((Float)0, rgb[3 * offset + 2] * invWt);
        }

        // Add splat value at pixel
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                             pixel.splatXYZ[2]};
        XYZToRGB(splatXYZ, splatRGB);
        rgb[3 * offset] += splatScale * splatRGB[0];
        rgb[3 * offset + 1] += splatScale * splatRGB[1];
        rgb[3 * offset + 2] += splatScale * splatRGB[2];

        // Scale pixel value by _scale_
        rgb[3 * offset] *= scale;
        rgb[3 * offset + 1] *= scale;
        rgb[3 * offset + 2] *= scale;
        ++offset;
    }

    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
    pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}


void Film::WriteImageTemp(int index, Float splatScale) {
    
    //////////////////////
    // PrISE-3D Updates //
    //////////////////////

    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[3 * offset]);

        // Normalize pixel with weight sum
        Float filterWeightSum = pixel.filterWeightSum;
        if (filterWeightSum != 0) {
            Float invWt = (Float)1 / filterWeightSum;
            rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
            rgb[3 * offset + 1] =
                std::max((Float)0, rgb[3 * offset + 1] * invWt);
            rgb[3 * offset + 2] =
                std::max((Float)0, rgb[3 * offset + 2] * invWt);
        }

        // Add splat value at pixel
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                             pixel.splatXYZ[2]};
        XYZToRGB(splatXYZ, splatRGB);
        rgb[3 * offset] += splatScale * splatRGB[0];
        rgb[3 * offset + 1] += splatScale * splatRGB[1];
        rgb[3 * offset + 2] += splatScale * splatRGB[2];

        // Scale pixel value by _scale_
        rgb[3 * offset] *= scale;
        rgb[3 * offset + 1] *= scale;
        rgb[3 * offset + 2] *= scale;
        ++offset;
    }

    // Write RGB image

    // define delimiter to split image name
    std::string delimiter = ".";
    std::string output_folder = PbrtOptions.folder;
    
    // find prefix and postfix information from `filename`
    std::string filename_prefix = filename.substr(0, filename.find(delimiter));
    std::string filename_postfix = filename.substr(filename.find(delimiter), filename.length());

    // create custom image
    std::string indexStr(std::to_string(index));

    while(indexStr.length() < PbrtOptions.digits){
        indexStr = "0" + indexStr;
    }

    // build folder
    std::string folder_image = std::string(output_folder + "/" + filename_prefix);
    std::string temp_filename= output_folder + "/" + filename_prefix + "/" + filename_prefix+ "-S" + std::to_string(PbrtOptions.samples) + "-" + indexStr + filename_postfix;
    
    // TODO : improve (recursively create folders)
    mkdir(output_folder.c_str(), 0775);
    mkdir(folder_image.c_str(), 0775);

    pbrt::WriteImage(temp_filename, &rgb[0], croppedPixelBounds, fullResolution);

    // check if necessary to write zbuffer (only at first image)
    // if(PbrtOptions.zbuffer && index == 0){
    //     std::string zbuffer_filename= output_folder + "/" + filename_prefix + "/" + filename_prefix+ "-S" + std::to_string(PbrtOptions.samples) + "-zbuffer" + filename_postfix;        

    //     std::unique_ptr<Float[]> zbufferFloat(new Float[croppedPixelBounds.Area()]);
    //     int offset = 0;

    //     Float maxValue = Float(0.);

    //     for (Point2i p : croppedPixelBounds) {
    //         // Convert pixel XYZ color to RGB
    //         Float &bufferPoint = GetBufferPoint(p);

    //         if (bufferPoint > maxValue && bufferPoint != Infinity){
    //             maxValue = Float(bufferPoint);
    //         }

    //         // Scale pixel value by _scale_
    //         zbufferFloat[offset] = bufferPoint;
    //         ++offset;
    //     }

    //     // for each Infinity value update value using max current value of zbuffer
    //     for(int i = 0; i < croppedPixelBounds.Area(); ++i){

    //         if (zbufferFloat[i] == Infinity){
    //             zbufferFloat[i] = Float(maxValue);
    //         }
    //     }

    //     pbrt::WriteImage(zbuffer_filename, &zbufferFloat[0], croppedPixelBounds, fullResolution, 1);
    // }   

    // // check if necessary to write normals (only at first image)
    // if(PbrtOptions.normals && index == 0){
    //     std::string normals_filename= output_folder + "/" + filename_prefix + "/" + filename_prefix+ "-S" + std::to_string(PbrtOptions.samples) + "-normals" + filename_postfix;

    //     std::unique_ptr<Float[]> normalsFloat(new Float[3 * croppedPixelBounds.Area()]);
    //     int offset = 0;
    //     for (Point2i p : croppedPixelBounds) {
    //         // Convert pixel XYZ color to RGB
    //         Normal3f &normal = GetNormalPoint(p);

    //         // Scale pixel value by _scale_
    //         normalsFloat[3 * offset] = (normal.x + 1) * 0.5;
    //         normalsFloat[3 * offset + 1] = (normal.y + 1) * 0.5;
    //         normalsFloat[3 * offset + 2] = (normal.z + 1) * 0.5;
    //         ++offset;
    //     }

    //     pbrt::WriteImage(normals_filename, &normalsFloat[0], croppedPixelBounds, fullResolution);
    // }

    //////////////////////////
    // End PrISE-3D Updates //
    //////////////////////////
}

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    std::string filename;
    if (PbrtOptions.imageFile != "") {
        filename = PbrtOptions.imageFile;
        std::string paramsFilename = params.FindOneString("filename", "");
        if (paramsFilename != "")
            Warning(
                "Output filename supplied on command line, \"%s\" is overriding "
                "filename provided in scene description file, \"%s\".",
                PbrtOptions.imageFile.c_str(), paramsFilename.c_str());
    } else
        filename = params.FindOneString("filename", "pbrt.exr");

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop;
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);
    else
        crop = Bounds2f(Point2f(Clamp(PbrtOptions.cropWindow[0][0], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][0], 0, 1)),
                        Point2f(Clamp(PbrtOptions.cropWindow[0][1], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][1], 0, 1)));

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                   Infinity);
    return new Film(Point2i(xres, yres), crop, std::move(filter), diagonal,
                    filename, scale, maxSampleLuminance);
}

}  // namespace pbrt
