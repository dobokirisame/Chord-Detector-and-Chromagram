#include "Chromagram.h"
#include <algorithm>

void Chromagram::generateNoteFrequencies() const noexcept {
    for (size_t i = 0; i < kNotesCount; i++) {
        kNoteFrequencies[i] = kRreferenceFrequency * pow (2,(static_cast<float>(i) / kNotesCount));
    }
}

Chromagram::Chromagram (size_t frameSize, size_t samplingFrequency)
     : mNumSamplesSinceLastCalculation(0),
       mChromaCalculationInterval(kBufferSize/2),
       mChromaReady(false),
       mSamplingFrequency(samplingFrequency) {
    generateNoteFrequencies();
    setupFFT();
    mBuffer.resize(kBufferSize);
    mChromagram.resize(kNotesCount, 0.0);
    magnitudeSpectrum.resize((kBufferSize/2)+1);
    makeHammingWindow();
    setInputAudioFrameSize (frameSize);
}

Chromagram::~Chromagram() {
    // free the Kiss FFT configuration
    free(mFftConfig);
}

void Chromagram::processAudioFrame (double* inputAudioFrame) {
    // create a vector
    std::vector<double> v;
    // use the array to initialise it
    v.assign (inputAudioFrame, inputAudioFrame + mInputAudioFrameSize);
    // process the vector
    processAudioFrame (v);
}

void Chromagram::processAudioFrame (std::vector<double> inputAudioFrame)
{
    // our default state is that the chroma is not ready
    mChromaReady = false;
    // downsample the input audio frame by 4
    downSampleFrame(inputAudioFrame);
    // move samples back
    for (size_t i = 0; i < kBufferSize - downSampledAudioFrameSize; i++) {
        mBuffer[i] = mBuffer[i + downSampledAudioFrameSize];
    }
    // add new samples to buffer
    for (size_t i = (kBufferSize - downSampledAudioFrameSize), n=0; i < kBufferSize; i++, n++)
    {
        mBuffer[i] = downsampledInputAudioFrame[n];
    }
    // add number of samples from calculation
    mNumSamplesSinceLastCalculation += mInputAudioFrameSize;
        
    // if we have had enough samples
    if (mNumSamplesSinceLastCalculation >= mChromaCalculationInterval) {
        // calculate the chromagram
        calculateChromagram();
        // reset num samples counter
        mNumSamplesSinceLastCalculation = 0;
    }
}

void Chromagram::setInputAudioFrameSize(size_t frameSize) {
    mInputAudioFrameSize = frameSize;
    downsampledInputAudioFrame.resize (mInputAudioFrameSize / 4);
    downSampledAudioFrameSize = downsampledInputAudioFrame.size();
}

void Chromagram::setSamplingFrequency(size_t samplingFrequency) {
    mSamplingFrequency = samplingFrequency;
}

void Chromagram::setChromaCalculationInterval (int numSamples) {
    mChromaCalculationInterval = numSamples;
}

std::vector<double> Chromagram::getChromagram() {
    return mChromagram;
}

bool Chromagram::isReady() {
    return mChromaReady;
}

void Chromagram::setupFFT() {
    mFftConfig = kiss_fft_alloc (kBufferSize, 0, nullptr, nullptr);
}


void Chromagram::calculateChromagram() {
    calculateMagnitudeSpectrum();
    double divisorRatio = (static_cast<double>(mSamplingFrequency) / 4.0) / static_cast<double>(kBufferSize);
    for (size_t n = 0; n < kNotesCount; n++) {
        double chromaSum = 0.0;
        for (size_t octave = 1; octave <= kOctavesCount; octave++) {
            double noteSum = 0.0;
            for (size_t harmonic = 1; harmonic <= kHarmonicsCount; harmonic++) {
                size_t centerBin = static_cast<size_t>(round ((kNoteFrequencies.at(n) * octave * harmonic) / divisorRatio));
                size_t minBin = centerBin - (kBinsCount * harmonic);
                size_t maxBin = centerBin + (kBinsCount * harmonic);
                auto maxVal = std::max_element(magnitudeSpectrum.begin() + static_cast<int>(minBin),
                                               magnitudeSpectrum.begin() + static_cast<int>(maxBin));
                noteSum += (*maxVal / static_cast<double>(harmonic));
            }
            chromaSum += noteSum;
        }
        mChromagram[n] = chromaSum;
    }
    mChromaReady = true;
}

void Chromagram::calculateMagnitudeSpectrum() {
    for (size_t i = 0;i < kBufferSize; i++) {
        mFftInputData[i].r = static_cast<float>(mBuffer.at(i) * window.at(i));
        mFftInputData[i].i = 0.0;
    }
    
    // execute kiss fft
    kiss_fft (mFftConfig, mFftInputData.data(), mFftOutputData.data());
    
    // compute first (N/2)+1 mag values
    for (size_t i = 0; i < (kBufferSize / 2) + 1; i++)
    {
        magnitudeSpectrum[i] = sqrt (pow (mFftOutputData[i].r, 2) + pow (mFftOutputData[i].i, 2));
        magnitudeSpectrum[i] = sqrt (magnitudeSpectrum[i]);
    }
}

void Chromagram::downSampleFrame (std::vector<double> inputAudioFrame) {
    std::vector<double> filteredFrame(mInputAudioFrameSize);
    double b0,b1,b2,a1,a2;
    double x1,x2,y1,y2;
    b0 = 0.2929;
    b1 = 0.5858;
    b2 = 0.2929;
    a1 = -0.0000;
    a2 = 0.1716;
    x1 = 0;
    x2 = 0;
    y1 = 0;
    y2 = 0;
    
    for (size_t i = 0; i < mInputAudioFrameSize; i++) {
        filteredFrame[i] = inputAudioFrame.at(i) * b0 + x1 * b1 + x2 * b2 - y1 * a1 - y2 * a2;
        x2 = x1;
        x1 = inputAudioFrame.at(i);
        y2 = y1;
        y1 = filteredFrame.at(i);
    }
    for (size_t i = 0; i < mInputAudioFrameSize/4; i++) {
        downsampledInputAudioFrame[i] = filteredFrame[i * 4];
    }
}

void Chromagram::makeHammingWindow() {
    window.clear();
    window.reserve(kBufferSize);
    for (size_t i = 0; i < kBufferSize;i++) {
        const auto windowValue = 0.54 - 0.46 * cos (2*M_PI*(static_cast<double>(i)/static_cast<double>(kBufferSize)));
        window.emplace_back(windowValue);
    }
}

double Chromagram::round (double val) {
    return floor (val + 0.5);
}
