/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// this is an abstract class for sound modification modules
#pragma once

#include <map>
#include <string>


/**
 * A base class for realtime sound modification
 */
class modbase {
public:
    modbase() {
        sample_rate_ = 48000; // default samplerate
        num_channels_ = 1; // should manually guarantee that num_channels equals to the size of bufferData
    };

    virtual ~modbase() {}; // make sure there is not any leaks

    /** 
     * in-place approach for realtime processing
     * 
     * we should manually guarantee we have enough data in bufferData
     * 
     * @param bufferData contains a 2D array, where the first points to channel, the second points to data
     * @param num_samples total number of available samples per channel
     */
    virtual void processBlock(float *const * bufferData, int num_samples) = 0;

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    virtual void setParams(std::map<std::string, float> params) = 0;

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    virtual void getParams(std::map<std::string, float> &params) = 0;

    /**
     * indicate whether the output data is ready
     * @return whether the output data is ready
     */
    virtual bool outputReady() {return true;}

protected:
    int sample_rate_; //!< sample rate
    int num_channels_; //!< number of channels
};
/** \example main.cc
 * This is an example of this class
 */


/**
 * A base class for offline sound modification
 */
class modbase_offline {
public:
    modbase_offline() {
        sample_rate_ = 48000; // default samplerate
        num_channels_ = 1; // should manually guarantee that num_channels equals to the size of bufferData
        num_res_ = 0;
    };
    virtual ~modbase_offline() {};

    /**
     * process the input data for non-realtime process (such as time-stretching)
     * @param inData contains a 2D array, where the first points to channel, the second points to data
     * @param num_in_samples total number of available samples per channel
     */
    virtual void processInData (float *const * inData, int num_in_samples) = 0; // for non-real-time process

    /**
     * get the output data after processInData()
     * @param outData returns a 2D array, where the first points to channel, the second points to data
     * @param num_out_samples total number of available samples per channel to retreive, if
     *                this is larger than what is available, outData only got what is available
     */
    virtual void getOutData(float *const * outData, int num_out_samples) = 0;

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    virtual void setParams(std::map<std::string, float> params) = 0;

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    virtual void getParams(std::map<std::string, float> &params) = 0;

    /**
     * @return number of output samples available at this time
     */
    virtual int getOutSamples() const {return num_res_;}

    /**
     * indicate whether the output data is ready
     * @return whether the output data is ready
     */
    virtual bool outputReady() {return true;}
    
protected:
    int sample_rate_; //!< sample rate
    int num_channels_; //!< number of channels
    int num_res_; //!< number of resulting samples after one process
};
/** \example main.cc
 * This is an example of this class
 */

/**
 * A base class for sound analyzer
 */
class modbase_analyzer {
public:
    modbase_analyzer() {
        sample_rate_ = 48000; // default samplerate
        num_channels_ = 1; // should manually guarantee that num_channels equals to the size of bufferData
        num_res_ = 0;
    };
    virtual ~modbase_analyzer() {};

    /**
     * process the input data for non-realtime process (such as time-stretching)
     * @param inData contains a 2D array, where the first points to channel, the second points to data
     * @param num_in_samples total number of available samples per channel
     */
    virtual void processInData (float *const * inData, int num_in_samples) = 0; // for non-real-time process

    /**
     * get the output data after processInData()
     * @param outData returns a 2D array of features, where the first points to channel, the second points to feature
     * @param num_out_samples total number of available symbols per channel to retreive, if
     *                this is larger than what is available, outData only got what is available
     */
    virtual void getOutData(float *const * outData, int num_out_symbols) = 0;

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    virtual void setParams(std::map<std::string, float> params) = 0;

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    virtual void getParams(std::map<std::string, float> &params) = 0;

    /**
     * @return number of output samples available at this time
     */
    virtual int getOutSymbols() const {return num_res_;}

    /**
     * @return a global scalar measurement for this metrics
     */
    virtual float getScalarMeasurement() const = 0;

    /**
     * indicate whether the output data is ready
     * @return whether the output data is ready
     */
    virtual bool outputReady() {return true;}
    
protected:
    int sample_rate_; //!< sample rate
    int num_channels_; //!< number of channels
    int num_res_; //!< number of resulting symbols after one process
};

/**
 * A base class for sound monitoring
 */
class modbase_meter {
public:
    modbase_meter() {
        sample_rate_ = 48000; // default samplerate
        num_channels_ = 1; // should manually guarantee that num_channels equals to the size of bufferData
    };

    virtual ~modbase_meter() {}; // make sure there is not any leaks

    /**
     * @return a global scalar measurement for this metrics
     */
    virtual float getScalarMeasurement() const = 0;

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    virtual void setParams(std::map<std::string, float> params) = 0;

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    virtual void getParams(std::map<std::string, float> &params) = 0;

    /** 
     * in-place approach for realtime processing
     * 
     * we should manually guarantee we have enough data in bufferData
     * 
     * @param bufferData contains a 2D array, where the first points to channel, the second points to data
     * @param num_samples total number of available samples per channel
     */
    virtual void processBlock(float *const * bufferData, int num_samples) = 0;

protected:
    int sample_rate_; //!< sample rate
    int num_channels_; //!< number of channels
};
/** \example main.cc
 * This is an example of this class
 */

