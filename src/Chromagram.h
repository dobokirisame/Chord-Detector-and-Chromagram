#ifndef CHROMAGRAM_H
#define CHROMAGRAM_H

#include <math.h>
#include <vector>

#ifdef USE_KISS_FFT
#include "kiss_fft.h"
#endif
#include <array>
namespace {
constexpr size_t kNotesCount(12);
constexpr size_t kBufferSize(8192);
constexpr size_t kHarmonicsCount(2);
constexpr size_t kOctavesCount(2);
constexpr size_t kBinsCount(2);
constexpr double kRreferenceFrequency(130.81278265); // Note C3 MIDI number 48 Lenght 1mm
std::array<double, 12> kNoteFrequencies;
/**
  http://us-satellite.net/nasa/endeavor/resources/mathdocs/Mathematics%20of%20SoundY.pdf
*/
} //namespace

//=======================================================================
/** A class for calculating a Chromagram from input audio
 * in a real-time context */
class Chromagram
{
    
public:
    /** Constructor
     * @param frameSize the input audio frame size 
     * @param fs the sampling frequency
     */
    Chromagram (size_t frameSize, size_t mSamplingFrequency);

    /** Destructor */
    ~Chromagram();
    
    /** Process a single audio frame. This will determine whether enough samples
     * have been accumulated and if so, will calculate the chromagram
     * @param inputAudioFrame an array containing the input audio frame. This should be
     * the length indicated by the input audio frame size passed to the constructor
     * @see setInputAudioFrameSize
     */
    void processAudioFrame (double* inputAudioFrame);
    
    /** Process a single audio frame. This will determine whether enough samples
     * have been accumulated and if so, will calculate the chromagram
     * @param inputAudioFrame a vector containing the input audio frame. This should be
     * the length indicated by the input audio frame size passed to the constructor
     * @see setInputAudioFrameSize
     */
    void processAudioFrame (std::vector<double> inputAudioFrame);
    
    /** Sets the input audio frame size
     * @param frameSize the input audio frame size
     */
    void setInputAudioFrameSize (size_t frameSize);
    
    /** Set the sampling frequency of the input audio
     * @param fs the sampling frequency in Hz
     */
    void setSamplingFrequency(size_t samplingFrequency);
    
    /** Set the interval at which the chromagram is calculated. As the algorithm requires
     * a significant amount of audio to be accumulated, it may be desirable to have the algorithm
     * not calculate the chromagram at every new audio frame. This function allows you to set the 
     * interval at which the chromagram will be calculated, specified in the number of samples at
     * the audio sampling frequency 
     * @param numSamples the number of samples that the algorithm will receive before calculating a new chromagram
     */
    void setChromaCalculationInterval (int numSamples);
    
    /** @returns the chromagram vector */
    std::vector<double> chromagram();
    
    /** @returns true if a new chromagram vector has been calculated at the current iteration. This should
     * be called after processAudioFrame
     */
    bool isReady();
    
private:
    void generateNoteFrequencies() const noexcept;
    void setupFFT();
    void calculateChromagram();
    void calculateMagnitudeSpectrum();
    std::vector<double> &&downSampleFrame(std::vector<double> inputAudioFrame) const;
    void makeHammingWindow();
private:
    std::array<double, kBufferSize> mHammingWindow;
    std::vector<double> mBuffer;
    std::vector<double> mMagnitudeSpectrum;
    std::vector<double> mChromagram;
    size_t mInputAudioFrameSize;
    size_t downSampledAudioFrameSize;
    int mNumSamplesSinceLastCalculation;
    int mChromaCalculationInterval;
    bool mChromaReady;
    size_t mSamplingFrequency;
#ifdef USE_KISS_FFT
    kiss_fft_cfg mFftConfig;
    std::array<kiss_fft_cpx, kBufferSize> mFftInputData;
    std::array<kiss_fft_cpx, kBufferSize> mFftOutputData;
#endif
};

#endif
