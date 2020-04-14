
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
#include <bits/stdc++.h> 

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
    zbuffer = std::unique_ptr<Float[]>(new Float[croppedPixelBounds.Area()]);
    normals = std::unique_ptr<Normal3f[]>(new Normal3f[croppedPixelBounds.Area()]); 

    for (int i = 0; i < croppedPixelBounds.Area(); ++i){
        zbuffer[i] = 0; // default value
    }
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
void Film::ApplyDL(FilmTile* tile, torch::jit::script::Module module) {

    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Apply DL and merge film tile" << tile->pixelBounds << " using " << PbrtOptions.DLConfidence << "%";
    std::lock_guard<std::mutex> lock(mutex);

    // build model input
    // expected (1, 7, 32, 32)
    // 1. mean pixel values (preprocessing: log + normalization)
    // 2. normals values (preprocessing: normalization)
    // 3. zbuffer values (preprocessing: normalization)

    // check tile size
    int nChannels = 7;
    int tileSize = 32;
    int nChannelValues = tileSize * tileSize;

    Bounds2i bounds = tile->GetPixelBounds();

    int widthExtra = (bounds[1].x - bounds[0].x) % tileSize;
    int heightExtra = (bounds[1].y - bounds[0].y) % tileSize;

    int widthBounds = (int)(widthExtra / 2);
    int heightBounds = (int)(heightExtra / 2);

    std::vector<float> inputValues(nChannels * nChannelValues);
    std::cout << "Input values " << nChannels * nChannelValues << std::endl;

    // store max values for mean input and zBuffer
    Float maxLogXYZ = 0;
    Float maxZbuffer = 0;

    int pixelCounter = 0;
    double epsilon = std::numeric_limits<double>::epsilon();

    // extract input model data
    for (Point2i pixel : tile->GetPixelBounds()) {
        
        // check out of bounds for model input (get the center of tile because tile is not really of size 32, exemple: 34 x 32 for interpolation)
        if (pixel.x < (widthBounds + bounds[0].x) || pixel.x > (bounds[1].x - widthBounds - 1))
            continue;

        if (pixel.y < (heightBounds + bounds[0].y) || pixel.y > (bounds[1].y - heightBounds - 1))
            continue;

        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Normal3f &normal = GetNormalPoint(pixel);
        Float &bufferPoint = GetBufferPoint(pixel);

        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);

        // fill XYZ values
        Float filterSum = mergePixel.filterWeightSum + tilePixel.filterWeightSum;

        for (int i = 0; i < 3; ++i){
            // normalize channel and apply log on it
            Float currentChannel = (float)(log10((mergePixel.xyz[i] + xyz[i] + epsilon) / filterSum));
            inputValues.at(pixelCounter + i * nChannelValues) = currentChannel;
            
            if (currentChannel > maxLogXYZ){
                maxLogXYZ = currentChannel;
            }
        }
        
        // fill normal values and normalize
        inputValues.at(pixelCounter + 3 * nChannelValues) = (float)((normal.x + 1) * 0.5);
        inputValues.at(pixelCounter + 4 * nChannelValues) = (float)((normal.y + 1) * 0.5);
        inputValues.at(pixelCounter + 5 * nChannelValues) = (float)((normal.z + 1) * 0.5);

        // fill zBuffer value
        inputValues.at(pixelCounter + 6 * nChannelValues) = (float)(bufferPoint);

        pixelCounter++;
    }

    // normalize model input data if necessary
    pixelCounter = 0;

    for (Point2i pixel : tile->GetPixelBounds()) {

        // check out of bounds for model input (get the center of tile because tile is not really of size 32, exemple: 34 x 32 for interpolation)
        if (pixel.x < (widthBounds + bounds[0].x) || pixel.x > (bounds[1].x - widthBounds - 1))
            continue;

        if (pixel.y < (heightBounds + bounds[0].y) || pixel.y > (bounds[1].y - heightBounds - 1))
            continue;

        for (int i = 0; i < 3; ++i){
            // normalize channel and apply log on it
            inputValues.at(pixelCounter + i * nChannelValues) /= ((float)maxLogXYZ);
        }

        // normalize zBuffer value
        inputValues.at(pixelCounter + 6 * nChannelValues) /= ((float)maxZbuffer);

        pixelCounter++;
    }

    // Create a vector of inputs for load model
    std::vector<torch::jit::IValue> inputs;
    
    auto inputData = torch::tensor(inputValues);
    auto inputDataReshaped = inputData.view({1, nChannels, tileSize, tileSize});
  
    // add input data to input IValue vector for model
    inputs.push_back(inputDataReshaped);

    // Execute the model and turn its output into a tensor.
    at::Tensor output = module.forward(inputs).toTensor();

    // Read output data
    std::vector<float> xv;
    
    size_t x_size = output.numel();

    // reading output data from buffer
    auto outputData = static_cast<float*>(output.data_ptr<float>());
    for(size_t i = 0; i < x_size; i++)
    {
      xv.push_back(outputData[i]);
    }

    // TODO :
    // renormalize (normalization * maxLog + exp) output model data
    // affects model data with percent of confidence
    // reinit filterWidth after applying output on film normalization
}
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
    if(PbrtOptions.zbuffer && index == 0){
        std::string zbuffer_filename= output_folder + "/" + filename_prefix + "/" + filename_prefix+ "-S" + std::to_string(PbrtOptions.samples) + "-zbuffer" + filename_postfix;        

        std::unique_ptr<Float[]> zbufferFloat(new Float[croppedPixelBounds.Area()]);
        int offset = 0;

        Float maxValue = Float(0.);

        for (Point2i p : croppedPixelBounds) {
            // Convert pixel XYZ color to RGB
            Float &bufferPoint = GetBufferPoint(p);

            if (bufferPoint > maxValue && bufferPoint != Infinity){
                maxValue = Float(bufferPoint);
            }

            // Scale pixel value by _scale_
            zbufferFloat[offset] = bufferPoint;
            ++offset;
        }

        // for each Infinity value update value using max current value of zbuffer
        for(int i = 0; i < croppedPixelBounds.Area(); ++i){

            if (zbufferFloat[i] == Infinity){
                zbufferFloat[i] = Float(maxValue);
            }
        }

        pbrt::WriteImage(zbuffer_filename, &zbufferFloat[0], croppedPixelBounds, fullResolution, 1);
    }   

    // check if necessary to write normals (only at first image)
    if(PbrtOptions.normals && index == 0){
        std::string normals_filename= output_folder + "/" + filename_prefix + "/" + filename_prefix+ "-S" + std::to_string(PbrtOptions.samples) + "-normals" + filename_postfix;

        std::unique_ptr<Float[]> normalsFloat(new Float[3 * croppedPixelBounds.Area()]);
        int offset = 0;
        for (Point2i p : croppedPixelBounds) {
            // Convert pixel XYZ color to RGB
            Normal3f &normal = GetNormalPoint(p);

            // Scale pixel value by _scale_
            normalsFloat[3 * offset] = (normal.x + 1) * 0.5;
            normalsFloat[3 * offset + 1] = (normal.y + 1) * 0.5;
            normalsFloat[3 * offset + 2] = (normal.z + 1) * 0.5;
            ++offset;
        }

        pbrt::WriteImage(normals_filename, &normalsFloat[0], croppedPixelBounds, fullResolution);
    }

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
