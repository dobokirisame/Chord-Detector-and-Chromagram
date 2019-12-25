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
    mMagnitudeSpectrum.resize((kBufferSize/2)+1);
    makeHammingWindow();
    setInputAudioFrameSize (frameSize);
}

Chromagram::~Chromagram() {
    free(mFftConfig);
}

void Chromagram::processAudioFrame (double* inputAudioFrame) {
    std::vector<double> v;
    v.assign (inputAudioFrame, inputAudioFrame + mInputAudioFrameSize);
    processAudioFrame(v);
}

void Chromagram::processAudioFrame (std::vector<double> inputAudioFrame) {
    //TODO: this shit must be refactored
    // it's just soooo bad
    mChromaReady = false;
    auto downsampledInputAudioFrame = std::move(downSampleFrame(inputAudioFrame));
    // move samples back
    for (size_t i = 0; i < kBufferSize - downsampledInputAudioFrame.size(); i++) {
        mBuffer[i] = mBuffer.at(i + downsampledInputAudioFrame.size());
    }
    // add new samples to buffer
    for (size_t i = (kBufferSize - downsampledInputAudioFrame.size()), n=0; i < kBufferSize; i++, n++)
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
}

void Chromagram::setSamplingFrequency(size_t samplingFrequency) {
    mSamplingFrequency = samplingFrequency;
}

void Chromagram::setChromaCalculationInterval (int numSamples) {
    mChromaCalculationInterval = numSamples;
}

std::vector<double> Chromagram::chromagram() {
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
                const size_t centerBin = static_cast<size_t>(std::round(kNoteFrequencies.at(n)*octave*harmonic/divisorRatio));
                const size_t minBin = centerBin - (kBinsCount * harmonic);
                const size_t maxBin = centerBin + (kBinsCount * harmonic);
                const auto maxVal = std::max_element(mMagnitudeSpectrum.begin() + static_cast<int>(minBin),
                                               mMagnitudeSpectrum.begin() + static_cast<int>(maxBin));
                noteSum += (*maxVal / static_cast<double>(harmonic));
            }
            chromaSum += noteSum;
        }
        mChromagram[n] = chromaSum;
    }
    mChromaReady = true;
}

void Chromagram::calculateMagnitudeSpectrum() {
    for (size_t i = 0; i < kBufferSize; i++) {
        mFftInputData[i] = {static_cast<float>(mBuffer.at(i) * mHammingWindow.at(i)), 0.0};
    }
    kiss_fft (mFftConfig, mFftInputData.data(), mFftOutputData.data());
    for (size_t i = 0; i < (kBufferSize / 2) + 1; i++) {
        mMagnitudeSpectrum[i] = pow(pow(mFftOutputData.at(i).r, 2) + pow(mFftOutputData.at(i).i, 2), -4);
    }
}

std::vector<double> &&Chromagram::downSampleFrame(std::vector<double> inputAudioFrame) const {
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
    
    for (size_t i=0; i<mInputAudioFrameSize; i++) {
        filteredFrame[i] = inputAudioFrame.at(i)*b0+x1*b1+x2*b2-y1*a1-y2*a2;
        x2 = x1;
        x1 = inputAudioFrame.at(i);
        y2 = y1;
        y1 = filteredFrame.at(i);
    }
    std::vector<double>  downsampledInputAudioFrame;
    downsampledInputAudioFrame.reserve(mInputAudioFrameSize/4);
    for (size_t i = 0; i < downsampledInputAudioFrame.capacity(); i++) {
        downsampledInputAudioFrame.emplace_back(filteredFrame.at(i*4));
    }
    return std::move(downsampledInputAudioFrame);
}

void Chromagram::makeHammingWindow() {
    for (size_t i = 0; i < mHammingWindow.size();i++) {
        mHammingWindow[i] = 0.54-0.46*cos(2*M_PI*static_cast<double>(i)/static_cast<double>(kBufferSize));
    }
}
