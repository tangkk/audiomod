////////////////////////////////////////////////////////////////////////////////
///
/// Classes for easy reading & writing of WAV sound files.
///
/// For big-endian CPU, define BIG_ENDIAN during compile-time to correctly
/// parse the WAV files with such processors.
///
/// Admittingly, more complete WAV reader routines may exist in public domain, but
/// the reason for 'yet another' one is that those generic WAV reader libraries are
/// exhaustingly large and cumbersome! Wanted to have something simpler here, i.e.
/// something that's not already larger than rest of the SoundTouch/SoundStretch program...
///
/// Author        : Copyright (c) Olli Parviainen
/// Author e-mail : oparviai 'at' iki.fi
/// SoundTouch WWW: http://www.surina.net/soundtouch
///
////////////////////////////////////////////////////////////////////////////////
//
// Last changed  : $Date: 2014-10-05 19:20:24 +0300 (Sun, 05 Oct 2014) $
// File revision : $Revision: 4 $
//
// $Id: wavfile.h 200 2014-10-05 16:20:24Z oparviai $
//
////////////////////////////////////////////////////////////////////////////////
//
// License :
//
//  SoundTouch audio processing library
//  Copyright (c) Olli Parviainen
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
////////////////////////////////////////////////////////////////////////////////

#ifndef WAVFILE_H
#define WAVFILE_H

#include <stdio.h>
#include <assert.h>
#include <vector>
#include <string>
#define ST_THROW_RT_ERROR(x)    {assert((const char *)x);}


#ifndef uint
typedef unsigned int uint;
#endif


/// WAV audio file 'riff' section header
typedef struct
{
    char riff_char[4];
    int  package_len;
    char wave[4];
} WavRiff;

/// WAV audio file 'format' section header
typedef struct
{
    char  fmt[4];
    int   format_len;
    short fixed;
    short channel_number;
    int   sample_rate;
    int   byte_rate;
    short byte_per_sample;
    short bits_per_sample;
} WavFormat;

/// WAV audio file 'fact' section header
typedef struct
{
    char  fact_field[4];
    int   fact_len;
    uint  fact_sample_len;
} WavFact;

/// WAV audio file 'data' section header
typedef struct
{
    char  data_field[4];
    uint  data_len;
} WavData;


/// WAV audio file header
typedef struct
{
    WavRiff   riff;
    WavFormat format;
    WavFact   fact;
    WavData   data;
} WavHeader;


/// Base class for processing WAV audio files.
class wavfileBase
{
private:
    /// Conversion working buffer;
    char *convBuff;
    int convBuffSize;
protected:
    wavfileBase();
    virtual ~wavfileBase();

    /// Get pointer to conversion buffer of at min. given size
    void *getConvBuffer(int sizeByte);
};


/// Class for reading WAV audio files.
class WavInFile : protected wavfileBase
{
private:
    /// File Name
    std::string name;

    /// File pointer.
    FILE *fptr;

    /// Position within the audio stream
    long position;

    /// Counter of how many bytes of sample data have been read from the file.
    long dataRead;

    /// WAV header information
    WavHeader header;

    /// Init the WAV file stream
    void init();

    /// Read WAV file headers.
    /// \return zero if all ok, nonzero if file format is invalid.
    int readWavHeaders();

    /// Checks WAV file header tags.
    /// \return zero if all ok, nonzero if file format is invalid.
    int checkCharTags() const;

    /// Reads a single WAV file header block.
    /// \return zero if all ok, nonzero if file format is invalid.
    int readHeaderBlock();

    /// Reads WAV file 'riff' block
    int readRIFFBlock();

    int m_dataStartPosInFile;//数据区开始的地方，随机读取时每次应先seek到数据区开始的地方

public:
    /// Constructor: Opens the given WAV file. If the file can't be opened,
    /// throws 'runtime_error' exception.
    explicit WavInFile(const char *filename);

    explicit WavInFile(FILE *file);

    /// Destructor: Closes the file.
    ~WavInFile();
    
    //fptr is null
    bool isValid(){
        if(fptr == NULL){
            return false;
        }
        else{
            return true;
        }
    }

    /// Rewind to beginning of the file
    void rewind();

    std::string getFileName() const;

    /// Get sample rate.
    uint getSampleRate() const;

    int getWavFormat() const;

    /// Get number of bits per sample, i.e. 8 or 16.
    uint getNumBits() const;

    /// Get sample data size in bytes. Ahem, this should return same information as
    /// 'getBytesPerSample'...
    uint getDataSizeInBytes() const;

    /// Get total number of samples in file.
    uint getNumSamples() const;

    /// Get number of bytes per audio sample (e.g. 16bit stereo = 4 bytes/sample)
    uint getBytesPerSample() const;

    /// Get number of audio channels in the file (1=mono, 2=stereo)
    uint getNumChannels() const;

    /// Get the audio file length in milliseconds
    uint getLengthMS() const;

    /// Returns how many milliseconds of audio have so far been read from the file
    ///
    /// \return elapsed duration in milliseconds
    uint getElapsedMS() const;

    /// Reads audio samples from the WAV file. This routine works only for 8 bit samples.
    /// Reads given number of elements from the file or if end-of-file reached, as many
    /// elements as are left in the file.
    ///
    /// \return Number of 8-bit integers read from the file.
    int read(unsigned char *buffer, int maxElems);

    /// Reads audio samples from the WAV file to 16 bit integer format. Reads given number
    /// of elements from the file or if end-of-file reached, as many elements as are
    /// left in the file.
    ///
    /// \return Number of 16-bit integers read from the file.
    int read(short *buffer,     ///< Pointer to buffer where to read data.
             int maxElems       ///< Size of 'buffer' array (number of array elements).
            );

    /// Reads audio samples from the WAV file to floating point format, converting
    /// sample values to range [-1,1[. Reads given number of elements from the file
    /// or if end-of-file reached, as many elements as are left in the file.
    /// Notice that reading in float format supports 8/16/24/32bit sample formats.
    ///
    /// \return Number of elements read from the file.
    int read(float *buffer,     ///< Pointer to buffer where to read data.
             int maxElems       ///< Size of 'buffer' array (number of array elements).
            );
    /// Reads audio samples form the wave file at time belone to  [startTime endTime] to floating point format
    /// time is in milllion second
    int read(float *buffer,
             long startTime,
             long endTime,
             int maxElems);
    int read(short *buffer, long startTime, long endTime, int maxElems);

    //read audio sampl from startIndex to endIndex, the index is in short cnt;
    int read(short *buffer, long startSampleIndex, long endSampleIndex);

    // read audio into interleave float buffer - tangkk
    // return the actual number of samples read
    int read(float **buffer, int num_samples);

    /// Check end-of-file.
    ///
    /// \return Nonzero if end-of-file reached.
    int eof() const;
};


#include <vector>
#include <string>
class WavOutFile;
extern std::vector<std::string> g_wavfileNameVec;
extern std::vector<WavOutFile*> g_WavOutFileObjVector;

/// Class for writing WAV audio files.
class WavOutFile : protected wavfileBase
{
private:
    /// Pointer to the WAV file
    FILE *fptr;

    /// WAV file header data.
    WavHeader header;

    /// Counter of how many bytes have been written to the file so far.
    int bytesWritten;

    /// Fills in WAV file header information.
    void fillInHeader(const uint sampleRate, const uint bits, const uint channels);

    /// Finishes the WAV file header by supplementing information of amount of
    /// data written to file etc
    void finishHeader();

    /// Writes the WAV file header.
    void writeHeader();

public:
    /// Constructor: Creates a new WAV file. Throws a 'runtime_error' exception
    /// if file creation fails.
    WavOutFile(const char *fileName,    ///< Filename
               int sampleRate,          ///< Sample rate (e.g. 44100 etc)
               int bits,                ///< Bits per sample (8 or 16 bits)
               int channels             ///< Number of channels (1=mono, 2=stereo)
              );

    WavOutFile(FILE *file, int sampleRate, int bits, int channels);

    /// Destructor: Finalizes & closes the WAV file.
    ~WavOutFile();

    /// Write data to WAV file. This function works only with 8bit samples.
    /// Throws a 'runtime_error' exception if writing to file fails.
    void write(const unsigned char *buffer, ///< Pointer to sample data buffer.
               int numElems                 ///< How many array items are to be written to file.
              );

    /// Write data to WAV file. Throws a 'runtime_error' exception if writing to
    /// file fails.
    void write(const short *buffer,     ///< Pointer to sample data buffer.
               int numElems             ///< How many array items are to be written to file.
              );

    /// Write data to WAV file in floating point format, saturating sample values to range
    /// [-1..+1[. Throws a 'runtime_error' exception if writing to file fails.
    void write(const float *buffer,     ///< Pointer to sample data buffer.
               int numElems             ///< How many array items are to be written to file.
              );
    void write(const double *buffer,     ///< Pointer to sample data buffer.
               int numElems             ///< How many array items are to be written to file.
    );

    void write(float *const * buffer, int num_samples);

    //把数组连续写入wav文件,写结束后需要调用EndWriteAllArray或EndWriteArray
    template <typename _T>
    static int WriteArrayTowavfile(const _T *inArray,int size,const char* fpath,float ratio=1,int samplerate=44100,int channel=1,int numBits=16) {
#ifdef DEBUG
        g_wavfileNameVec.size();
        _T *indata = new _T[size];
        for (int i=0; i<size; i++) {
            indata[i] = inArray[i]*ratio;
        }
        bool ffind= false;
        int idx =-1;
        std::string wavFilePath(fpath);
        for (int i=0; i<g_wavfileNameVec.size(); i++) {
            if(g_wavfileNameVec[i].compare(wavFilePath)==0){
                ffind=true;
                idx = i;
                break;
            }
        }
        if (ffind && idx!=-1 &&(g_WavOutFileObjVector.size()> idx)) {
            WavOutFile* wof = g_WavOutFileObjVector[idx];
            wof->write(indata, size);
        }else{
            g_wavfileNameVec.push_back(wavFilePath);
            WavOutFile *wof = new WavOutFile(wavFilePath.c_str(),samplerate,numBits,channel);
            g_WavOutFileObjVector.push_back(wof);
            wof->write(indata,size);
        }
        if(indata){
            delete [] indata;
            indata=NULL;
        }
#endif
        return 1;
    }
    //结束对所有array的写操作，可在程序退出前调用，确保所有arry正常写入wav文件
    static void EndWriteAllArray();
    //传入文件名，结束对某个array的写操作
    static void EndWriteArray(const char* fpath);
};

int audioRead(std::string &inputPath, int samprate, int nchan, std::vector<float> &input);
int audioWrite(std::string &outputPath, int samprate, int nchan, std::vector<float> &output);
#endif
