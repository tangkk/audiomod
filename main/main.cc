#include "wavfile.h"
#include <sys/time.h>
#include <fstream>
#include <cmath>
#include <cerrno>
#include <cstring>
#include "audiomod.h"

using namespace audiomod;

static double getTimeOfDay() {
    struct timeval  tv;
    gettimeofday(&tv, NULL);
    double time_in_mill =
    (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;
    return time_in_mill;
}

static inline bool exists_test (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

int main (int argc, char* argv[]) {   
    
    if (argv[1]==NULL){
        std::cerr << "usage: ./audiomod-exe dafx_name infile outfile <args> (dafx: "
            "constant, "

            // phase vocoder models
            "time_stretch, "
            "normal_pitchshift, "
            "formant_pitchshift, "
            "gender_change, "

            // vocoder models
            "vocoder, "
            "vocoder_chord, "

            // robotic and whispering
            "robotic, "
            "whisper, "

            // delay line models
            "vibrato, "
            "delay, "
            "flanger, "
            "chorus, "

            // modulation models
            "ringmod, "
            "tremolo, "

            // dynamics models
            "compressor, "
            "limiter, "
            "autogain, "

            // reverb models
            "reverb, "

            // filtering models
            "autowah, "
            "phaser, "

            // meters
            "loudnessmeter, "

            // equalizer
            "equalizer, "

            // gain
            "gain, "

            // analyzer
            // "vad, "
            // "makeupgain, "
            "envelope, "
            "pyin, "
            "chromatuning, "
            "hummingest, "
            "keydetection, "
            "chordestimate, "
            "f0tonote,"
            ")" << std::endl;
        return -1;
    }
    if (argv[2]==NULL){
        std::cerr<<"err: input required"<<std::endl;
        return -1;
    }
    if (argv[3]==NULL){
        std::cerr<<"err: output <arg> required"<<std::endl;
        return -1;
    }

    std::string model_name(argv[1]);
    std::string input_file_name(argv[2]);
    std::string output_file_name(argv[3]);

    if (!exists_test(input_file_name) && input_file_name[0] != '-') {
        std::cerr << "input_file doesn't exist" << std::endl;
        return -1;
    }
    
    //create input stream and reader
    WavInFile *input = nullptr;
    std::ifstream txtinput;
    if (input_file_name[0] == '-') {
        // pipe input
        std::freopen(nullptr, "rb", stdin);

        if (std::ferror(stdin))
            throw std::runtime_error(std::strerror(errno));

        std::cerr << "reading from stdin..." << std::endl;
        input = new WavInFile(stdin);
    } else {
        // local file input
        if (model_name == "f0tonote") {
            txtinput.open(input_file_name);
        } else {
            input = new WavInFile(input_file_name.c_str());
        }
    }

    int wav_format = input == nullptr ? 0 : input->getWavFormat();
    std::cerr << "wav_format = " << wav_format << std::endl;

    int bytes_per_sample = input == nullptr ? 0 : input->getBytesPerSample();
    std::cerr << "bytes_per_sample = " << bytes_per_sample << std::endl;

    // int data_len_in_bytes = input->getDataSizeInBytes();
    // std::cerr << "data_len_in_bytes = " << data_len_in_bytes << std::endl;
    
    // num_channels = reader->numChannels;
    int num_channels = input == nullptr ? 1 : input->getNumChannels();
    std::cerr << "num_channels = " << num_channels << std::endl;

    // sample_rate = reader->sampleRate;
    int sample_rate = input == nullptr ? 44100 : input->getSampleRate();
    std::cerr << "sample_rate = " <<  sample_rate << std::endl;

    int file_length = input == nullptr ? 0 : input->getNumSamples();
    if (input_file_name[0] == '-') { // assign file_length from outside (may be different from the actual length)
        float file_dur = atof(&input_file_name[2]);
        file_length = file_dur * sample_rate;
    }
    std::cerr << "file_length = " << file_length << std::endl;

    int output_bits_per_sample = 16;

    int target_LUFS = -15;

    WavOutFile *output = nullptr;
    std::ofstream txtoutput;
    if (model_name == "loudnessmeter" || model_name == "envelope" || model_name == "pyin" || model_name == "f0tonote"
        || model_name == "chromatuning" || model_name == "hummingest" 
        || model_name == "keydetection" || model_name == "chordestimate"
    ) {
        // txtoutput.open(output_file_name);
        // use std::cout as output
    } else {
        output = new WavOutFile(output_file_name.c_str(), sample_rate, output_bits_per_sample, num_channels);
    }

    int block_size = sample_rate / 100 < 480 ? 480 : sample_rate / 100;
    int extra_meta_param = 0; // extra parameter whenever necessary

    
    modbase * m_modbase = nullptr;
    modbase_offline * m_modbase_offline = nullptr;
    modbase_meter * m_modbase_meter = nullptr;
    modbase_analyzer * m_modbase_analyzer = nullptr;

    std::unique_ptr<phasevocoder> m_phasevocoder;

    std::unique_ptr<vibrato> m_vibrato;
    std::unique_ptr<chorus> m_chorus;
    std::unique_ptr<flanger> m_flanger;
    std::unique_ptr<delay> m_delay;

    std::unique_ptr<ringmod> m_ringmod;
    std::unique_ptr<tremolo> m_tremolo;

    std::unique_ptr<compressor> m_compressor;
    std::unique_ptr<limiter> m_limiter;
    
    std::unique_ptr<reverb> m_reverb;
    
    std::unique_ptr<autowah> m_autowah;
    std::unique_ptr<phaser> m_phaser;

    std::unique_ptr<loudnessmeter> m_loudnessmeter;

    std::unique_ptr<equalizer> m_equalizer;

    std::unique_ptr<gain> m_gain;

    std::unique_ptr<envelope> m_envelope;

    std::unique_ptr<pyin> m_pyin;
    std::unique_ptr<pyin> m_pyin_2;
    std::unique_ptr<chromatuning> m_chromatuning;
    std::unique_ptr<keydetection> m_keydetection;
    std::unique_ptr<chordestimate> m_chordestimate;
    std::unique_ptr<f0tonote> m_f0tonote;
    
    if (model_name == "constant") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, 0, CONSTANT));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "time_stretch") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (time_ratio, coremode, fftsize)"<<std::endl;
            return -1;
        }
        float time_ratio = atof(argv[4]);
        int coremode = atoi(argv[5]);
        int fftsize = atoi(argv[6]);
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, time_ratio, 0, NORMAL_STRETCH, coremode, fftsize));
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "normal_pitchshift") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (pitchshift_amount, coremode, fftsize)"<<std::endl;
            return -1;
        }
        float shift_pitch = atof(argv[4]);
        int coremode = atoi(argv[5]);
        int fftsize = atoi(argv[6]);
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, shift_pitch, NORMAL_SHIFT, coremode, fftsize));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "formant_pitchshift") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (pitchshift_amount, coremode, fftsize)"<<std::endl;
            return -1;
        }
        float shift_pitch = atof(argv[4]);
        int coremode = atoi(argv[5]);
        int fftsize = atoi(argv[6]);
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, shift_pitch, FORMANT_PRESERVE, coremode, fftsize));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "gender_change") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (pitchshift_amount, coremode, fftsize)"<<std::endl;
            return -1;
        }
        float shift_pitch = atof(argv[4]);
        int coremode = atoi(argv[5]);
        int fftsize = atoi(argv[6]);
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, shift_pitch, GENDER_CHANGE, coremode, fftsize));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "vocoder") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, 0, VOCODER_ROSENBERG));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "vocoder_chord") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, 0, VOCODER_CHORD));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "robotic") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, 0, ROBOTIC));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "whisper") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phasevocoder = std::unique_ptr<phasevocoder>(new phasevocoder(sample_rate, num_channels, 1, 0, WHISPER));
        m_modbase = m_phasevocoder.get();
        m_modbase_offline = m_phasevocoder.get();
    }
    else if (model_name == "vibrato") {
        if (argc < 6) {
            std::cerr<<"err: not enough para (sweepWidth, frequency)"<<std::endl;
            return -1;
        }
        float width = atof(argv[4]);
        float freq = atof(argv[5]);
        m_vibrato = std::unique_ptr<vibrato>(new vibrato(sample_rate, num_channels, width, freq));
        m_modbase = m_vibrato.get();
    }
    else if (model_name == "ringmod") {
        if (argc < 6) {
            std::cerr<<"err: not enough para (sweepWidth, frequency)"<<std::endl;
            return -1;
        }
        float width = atof(argv[4]);
        float freq = atof(argv[5]);
        m_ringmod = std::unique_ptr<ringmod>(new ringmod(sample_rate, num_channels, width, freq));
        m_modbase = m_ringmod.get();
    }
    else if (model_name == "tremolo") {
        if (argc < 6) {
            std::cerr<<"err: not enough para (frequency, depth)"<<std::endl;
            return -1;
        }
        float freq = atof(argv[4]);
        float depth = atof(argv[5]);
        m_tremolo = std::unique_ptr<tremolo>(new tremolo(sample_rate, num_channels, freq, depth));
        m_modbase = m_tremolo.get();
    }
    else if (model_name == "compressor") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (threshold, ratio, makeup)"<<std::endl;
            return -1;
        }
        float threshold = atof(argv[4]);
        float ratio = atof(argv[5]);
        float makeup = atof(argv[6]);
        m_compressor = std::unique_ptr<compressor>(new compressor(sample_rate, num_channels, threshold, ratio, makeup));
        m_modbase = m_compressor.get();
    }
    else if (model_name == "limiter") {
        if (argc < 6) {
            std::cerr<<"err: not enough para (threshold, makeup)"<<std::endl;
            return -1;
        }
        float threshold = atof(argv[4]);
        float makeup = atof(argv[5]);
        m_limiter = std::unique_ptr<limiter>(new limiter(sample_rate, num_channels, threshold, makeup));
        m_modbase = m_limiter.get();
    }
    else if (model_name == "reverb") {
        if (argc < 9) {
            std::cerr<<"err: not enough para (roomsize, damp, width, dry, wet)"<<std::endl;
            return -1;
        }
        float roomsize = atof(argv[4]);
        float damp = atof(argv[5]);
        float width = atof(argv[6]);
        float dry = atof(argv[7]);
        float wet = atof(argv[8]);
        // default: float roomsize = 0.8f, float damp = 0.9f, float width = 2.f, float dry = 0.9f, float wet = 0.1f (basically 0~1)
        // 0.8 0.9 0.5 0.9 0.1
        m_reverb = std::unique_ptr<reverb>(new reverb(sample_rate, num_channels, roomsize, damp, width, dry, wet));
        m_modbase = m_reverb.get();
    }
    else if (model_name == "autogain") {
        if (argc < 5) {
            std::cerr<<"err: not enough para (target_LUFS)"<<std::endl;
            return -1;
        }
        target_LUFS = atoi(argv[4]);
        std::cerr << "target_LUFS = " << target_LUFS << std::endl;
        m_loudnessmeter = std::unique_ptr<loudnessmeter>(new loudnessmeter(sample_rate, num_channels, block_size));
        m_limiter = std::unique_ptr<limiter>(new limiter(sample_rate, num_channels));
    }
    else if (model_name == "autowah") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_autowah = std::unique_ptr<autowah>(new autowah(sample_rate, num_channels));
        m_modbase = m_autowah.get();
    }
    else if (model_name == "loudnessmeter") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        txtoutput.open(output_file_name);
        m_loudnessmeter = std::unique_ptr<loudnessmeter>(new loudnessmeter(sample_rate, num_channels, block_size));
        m_modbase_meter = m_loudnessmeter.get();
    }
    else if (model_name == "envelope") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        txtoutput.open(output_file_name);
        m_envelope = std::unique_ptr<envelope>(new envelope(sample_rate, num_channels));
        m_modbase_analyzer = m_envelope.get();
    }
    else if (model_name == "pyin") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (blocksize, stepsize, onsetsens)"<<std::endl;
            return -1;
        }
        int blockSize = atoi(argv[4]);
        int stepSize = atoi(argv[5]);
        float pyinOnsetSens = atof(argv[6]);
        txtoutput.open(output_file_name);
        m_pyin = std::unique_ptr<pyin>(new pyin(sample_rate, num_channels, blockSize, stepSize, pyinOnsetSens));
        m_modbase_analyzer = m_pyin.get();
    }
    else if (model_name == "f0tonote") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (sample_rate, stepsize, onsetsens)"<<std::endl;
            return -1;
        }
        int sample_rate = atoi(argv[4]);
        int stepSize = atoi(argv[5]);
        float pyinOnsetSens = atof(argv[6]);
        txtoutput.open(output_file_name);
        m_f0tonote = std::unique_ptr<f0tonote>(new f0tonote(sample_rate, 0, 0, stepSize, pyinOnsetSens));
        m_modbase_analyzer = m_f0tonote.get();
    }
    else if (model_name == "chromatuning") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (blocksize, stepsize, durneeded)"<<std::endl;
            return -1;
        }
        // prefered params 8192, 8192
        int blockSize = atoi(argv[4]); 
        int stepSize = atoi(argv[5]);
        extra_meta_param = atoi(argv[6]);
        txtoutput.open(output_file_name);
        m_chromatuning = std::unique_ptr<chromatuning>(new chromatuning(sample_rate, num_channels, blockSize, stepSize));
        m_modbase_meter = m_chromatuning.get();
        block_size = blockSize;
    }
    else if (model_name == "keydetection") {
        if (argc < 7) {
            std::cerr<<"err: not enough para (blocksize, stepsize, durneeded)"<<std::endl;
            return -1;
        }
        // prefered params 8192, 8192
        int blockSize = atoi(argv[4]); 
        int stepSize = atoi(argv[5]);
        extra_meta_param = atoi(argv[6]);
        txtoutput.open(output_file_name);
        m_keydetection = std::unique_ptr<keydetection>(new keydetection(sample_rate, num_channels, blockSize, stepSize));
        m_modbase_analyzer = m_keydetection.get();
        block_size = blockSize;
    }
    else if (model_name == "chordestimate") {
        if (argc < 6) {
            std::cerr<<"err: not enough para (blocksize, stepsize)"<<std::endl;
            return -1;
        }
        // prefered params 8192, 8192
        int blockSize = atoi(argv[4]); 
        int stepSize = atoi(argv[5]);
        txtoutput.open(output_file_name);
        m_chordestimate = std::unique_ptr<chordestimate>(new chordestimate(sample_rate, num_channels, blockSize, stepSize));
        m_modbase_analyzer = m_chordestimate.get();
        block_size = blockSize;
    }
    else if (model_name == "hummingest") {
        if (argc < 9) {
            std::cerr<<"err: not enough para (blocksizes, stepsizes, onsetsenses, blocksize, stepsize, extrametap)"<<std::endl;
            return -1;
        }

        char* token = std::strtok(argv[4], ",");
        int pyinBlockSizes[2];
        int i = 0;
        while (token != nullptr && i < 2) {
            pyinBlockSizes[i++] = std::atoi(token);
            token = std::strtok(nullptr, ",");
        }

        token = std::strtok(argv[5], ",");
        int pyinStepSizes[2];
        i = 0;
        while (token != nullptr && i < 2) {
            pyinStepSizes[i++] = std::atoi(token);
            token = std::strtok(nullptr, ",");
        }

        token = std::strtok(argv[6], ",");
        float pyinOnsetSenses[2];
        i = 0;
        while (token != nullptr && i < 2) {
            pyinOnsetSenses[i++] = std::atof(token);
            token = std::strtok(nullptr, ",");
        }

        // int pyinBlockSize = atoi(argv[4]); 
        // int pyinStepSize = atoi(argv[5]);
        // float pyinOnsetSens = atof(argv[6]);

        // two different configurations
        m_pyin = std::unique_ptr<pyin>(new pyin(sample_rate, num_channels, pyinBlockSizes[0], pyinStepSizes[0], pyinOnsetSenses[0]));
        m_pyin_2 = std::unique_ptr<pyin>(new pyin(sample_rate, num_channels, pyinBlockSizes[1], pyinStepSizes[1], pyinOnsetSenses[1]));
        // m_modbase_analyzer = m_pyin.get();

        int tuningBlockSize = atoi(argv[7]); 
        int tuningStepSize = atoi(argv[8]);
        m_chromatuning = std::unique_ptr<chromatuning>(new chromatuning(sample_rate, num_channels, tuningBlockSize, tuningStepSize));
        block_size = tuningBlockSize; // use this as the system block size
        // m_modbase_meter = m_chromatuning.get();

        extra_meta_param = atoi(argv[9]); // for pyin1 and pyin2 selection

        txtoutput.open(output_file_name); // file for output

        // autogain modules
        m_loudnessmeter = std::unique_ptr<loudnessmeter>(new loudnessmeter(sample_rate, num_channels, block_size));
        m_limiter = std::unique_ptr<limiter>(new limiter(sample_rate, num_channels));

    }
    else if (model_name == "equalizer") {
        float *paramlist = nullptr;
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }

        // default:
        // 1 200 0.3 1
        // 0 400 0.3 -1.5 
        // 0 1000 0.3 1.5
        // 0 2000 0.3 1.5
        // 0 3000 0.3 1.5
        // 0 4000 0.3 1.5
        // 0 5000 0.3 -1.5
        // 0 6000 0.3 1.0
        // 8 groups (each 4 param): HighPassFilter, LowShelfFilter, Peaking 0-3, HighShelfFilter, LowPassFilter
        if (argc == 36) {
            std::cerr<<"there's a list of 32 params"<<std::endl;
            paramlist = new float[32];
            for (int i=0; i<32; i++) {
                paramlist[i] = atof(argv[i+4]);
                std::cerr << i << ":" << paramlist[i] << " | ";
            }
            std::cerr << std::endl;
        }
        m_equalizer = std::unique_ptr<equalizer>(new equalizer(sample_rate, num_channels, paramlist));
        m_modbase = m_equalizer.get();

        if (paramlist != nullptr) {
            delete[] paramlist;
            paramlist = nullptr;
        }
    }
    else if (model_name == "gain") {
        if (argc < 4) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        float initgain = atof(argv[4]);
        m_gain = std::unique_ptr<gain>(new gain(sample_rate, num_channels, initgain));
        m_modbase = m_gain.get();
    }
    else if (model_name == "chorus") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_chorus = std::unique_ptr<chorus>(new chorus(sample_rate, num_channels));
        m_modbase = m_chorus.get();
    }
    else if (model_name == "flanger") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_flanger = std::unique_ptr<flanger>(new flanger(sample_rate, num_channels, 0.01, 0.6, 0.6));
        m_modbase = m_flanger.get();
    }
    else if (model_name == "delay") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_delay = std::unique_ptr<delay>(new delay(sample_rate, num_channels, 0.3, 0.3, 0.3));
        m_modbase = m_delay.get();
    }
    else if (model_name == "phaser") {
        if (argc < 3) {
            std::cerr<<"err: not enough para ()"<<std::endl;
            return -1;
        }
        m_phaser = std::unique_ptr<phaser>(new phaser(sample_rate, num_channels));
        m_modbase = m_phaser.get();
    }
    else {
        std::cerr << "fx not supported or wrong fx!" << std::endl;
        std::cerr << "or maybe you should put dafx right after ./audiomod" << std::endl;
        return -1;
    }


    // allocate buffers
    std::cerr << "block_size = " << block_size << std::endl;

    float** buff;
    buff = new float* [num_channels];
    for (int i=0; i<num_channels; i++) {
        buff[i] = new float[block_size];
    }

    float** outbuff;
    outbuff = new float* [num_channels];
    for (int i=0; i<num_channels; i++) {
        outbuff[i] = new float[block_size * 4];
    }

    // printf("ready to process...\n");
    if (model_name == "time_stretch") {
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            m_modbase_offline->processInData(buff, num_samples);
            m_modbase_offline->getOutData(outbuff, m_modbase_offline->getOutSamples());
            output->write(outbuff, m_modbase_offline->getOutSamples());
        }
    }
    else if (model_name == "normal_pitchshift" || model_name == "formant_pitchshift" ||
        model_name == "gender_change" || model_name == "vocoder" || model_name == "vocoder_chord" ||
        model_name == "robotic" || model_name == "whisper") {

        int current_output_length = 0;
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            m_modbase_offline->processInData(buff, num_samples);
            // printf("num_samples:%d, num_out_samples:%d\n", num_samples, m_modbase_offline->getOutSamples());
            m_modbase_offline->getOutData(outbuff, m_modbase_offline->getOutSamples());
            output->write(outbuff, m_modbase_offline->getOutSamples());
            current_output_length += m_modbase_offline->getOutSamples();
        }
        printf("process remaining data\n");
        // process remaining data
        for (int i = 0; i < num_channels; i++) {
            memset(buff[i], 0, sizeof(float) * block_size);
        }
        while (current_output_length < file_length) {
            m_modbase_offline->processInData(buff, block_size);
            printf("block_size:%d, num_out_samples:%d\n", block_size, m_modbase_offline->getOutSamples());
            m_modbase_offline->getOutData(outbuff, m_modbase_offline->getOutSamples());
            if (file_length - current_output_length > m_modbase_offline->getOutSamples()) {
                output->write(outbuff, m_modbase_offline->getOutSamples());
                current_output_length += m_modbase_offline->getOutSamples();
            } else {
                int num_to_write = file_length - current_output_length;
                output->write(outbuff, num_to_write);
                current_output_length += num_to_write;
            }
        }
    }
    else if (model_name == "loudnessmeter") {
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            m_modbase_meter->processBlock(buff, num_samples);
        }
        // printf("loudness(dB):%f\n", m_modbase_meter->getScalarMeasurement());
        txtoutput << m_modbase_meter->getScalarMeasurement() << std::endl;
        float dbloudness = m_modbase_meter->getScalarMeasurement();
        std::cerr << "dbloudness(LUFS):" << dbloudness << std::endl;

    }
    else if (model_name == "envelope") {
        int envelope_block_size = sample_rate / 100; // corresponding to 10ms
        std::cerr << "envelope_block_size:" << envelope_block_size << std::endl;
        int envelope_time_step = 10; // in ms unit
        int curtime = 0; // in ms unit
        for (int i=0; i < file_length; i+= envelope_block_size) {
            int num_samples = input->read(buff, envelope_block_size);
            m_modbase_analyzer->processInData(buff, num_samples);
            m_modbase_analyzer->getOutData(outbuff, 1);
            float thisAmp = outbuff[0][0];

            txtoutput << curtime << "\t" << thisAmp << std::endl;
            curtime += envelope_time_step;
        }
        float envelope_mean = m_modbase_analyzer->getScalarMeasurement();
        std::cerr << "envelope_mean:" << envelope_mean << std::endl;
    }
    else if (model_name == "pyin") { // fully offline analyzing, real the whole file and process
        // init a buff of file length
        float** fullLengthbuff;
        fullLengthbuff = new float* [num_channels];
        for (int i=0; i<num_channels; i++) {
            fullLengthbuff[i] = new float[file_length];
        }

        int num_samples = input->read(fullLengthbuff, file_length);
        m_modbase_analyzer->processInData(fullLengthbuff, num_samples);

        for (int i=0; i<num_channels; i++) {
            delete[] fullLengthbuff[i];
        }
        delete[] fullLengthbuff;

        int numNotes = m_modbase_analyzer->getScalarMeasurement();
        
        if (numNotes < block_size * 4) {
            m_modbase_analyzer->getOutData(outbuff, numNotes);
            for (int i=0; i<numNotes*3; i+=3) {
                float st = outbuff[0][i];
                float dur = outbuff[0][i+1];
                float f0 = outbuff[0][i+2];
                if (txtoutput.is_open()) {
                    // std::cerr << "txtoutput is opened" << std::endl;
                    txtoutput << st << "," << dur << "," << f0 << std::endl;
                } else {
                    std::cerr << "txtoutput not opened" << std::endl;
                }
                
            }
        } else {
            std::cerr << "numNotes > block_size * 4, this is strange!" << std::endl;
        }

    }
    else if (model_name == "f0tonote") {
        std::vector<float> f0vec;
        std::vector<float> rmsvec;
        if (txtinput.is_open()) {
            std::string line;
            while (std::getline(txtinput, line)) {
                // Output the line
                // std::cerr << line << std::endl;
                std::stringstream ss(line);

                std::string temp;
                std::getline(ss, temp, ',');
                float rms = std::stof(temp);
                
                std::getline(ss, temp, ',');
                float f0 = std::stof(temp);
                
                rmsvec.push_back(rms);
                f0vec.push_back(f0);
                // std::cerr << "Extracted value: " << value << std::endl;
            }
        }
        // std::cerr << "txtinput done..." << std::endl;

        float* data[2];
        data[0] = f0vec.data();
        data[1] = rmsvec.data();
        m_modbase_analyzer->processInData(data, f0vec.size());
        // std::cerr << "processInData done..." << std::endl;
        int numNotes = m_modbase_analyzer->getScalarMeasurement();
        if (numNotes > 0) {
            m_modbase_analyzer->getOutData(outbuff, numNotes);
            // std::cerr << "getOutData done..." << std::endl;
            for (int i=0; i<numNotes*3; i+=3) {
                float st = outbuff[0][i];
                float dur = outbuff[0][i+1];
                float pitch = outbuff[0][i+2];
                if (txtoutput.is_open()) {
                    // std::cerr << "txtoutput is opened" << std::endl;
                    txtoutput << st << "," << dur << "," << pitch << std::endl;
                } else {
                    std::cerr << "txtoutput not opened" << std::endl;
                }
                
            }
        } else {
            std::cerr << "numNotes == 0, this is strange!" << std::endl;
        }

    }
    else if (model_name == "chromatuning") {
        int targetLen = extra_meta_param == -1 ? file_length : extra_meta_param * sample_rate;
        for (int i = 0; i < targetLen; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            if (num_samples < block_size) {
                // pad remaining buf with zeros
                for (int j=0; j<num_channels; j++) {
                    for (int k=num_samples; k<block_size; k++) {
                        buff[j][k] = 0;
                    }
                }
            }
            m_modbase_meter->processBlock(buff, num_samples);
        }
        float globaltuning = m_modbase_meter->getScalarMeasurement();
        txtoutput << 0 << "," << 0 << "," << globaltuning << std::endl;
    }
    else if (model_name == "keydetection") {
        int targetLen = extra_meta_param == -1 ? file_length : extra_meta_param * sample_rate;
        for (int i = 0; i < targetLen; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            if (num_samples < block_size) {
                // pad remaining buf with zeros
                for (int j=0; j<num_channels; j++) {
                    for (int k=num_samples; k<block_size; k++) {
                        buff[j][k] = 0;
                    }
                }
            }
            m_modbase_analyzer->processInData(buff, num_samples);
        }
        std::vector<std::string> labels;
        m_modbase_analyzer->getOutData(outbuff, 0, &labels);

        int num_keys = m_modbase_analyzer->getScalarMeasurement();
        std::cerr << "num_keys:" << num_keys << std::endl;
        for (int i=0; i<num_keys + 1; i++) { // there's an end label at the end
            txtoutput << labels[4*i] << "," << labels[4*i + 1] << "," << labels[4*i + 2] << "," << labels[4*i + 3] << std::endl;
        }
        // txtoutput << float(file_length) / sample_rate << "," << -1 << std::endl;
        // txtoutput << "------" << std::endl;
        // for (int i=0; i<24; i++) {
        //     txtoutput << outbuff[0][i] << std::endl;
        // }
    }
    else if (model_name == "chordestimate") {
        int targetLen = file_length;
        for (int i = 0; i < targetLen; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            if (num_samples < block_size) {
                // pad remaining buf with zeros
                for (int j=0; j<num_channels; j++) {
                    for (int k=num_samples; k<block_size; k++) {
                        buff[j][k] = 0;
                    }
                }
            }
            m_modbase_analyzer->processInData(buff, num_samples);
        }
        // printf("processInData done.\n");
        std::vector<std::string> labels;
        m_modbase_analyzer->getOutData(outbuff, 0, &labels);

        int num_chords = m_modbase_analyzer->getScalarMeasurement();
        std::cerr << "num_chords:" << num_chords << std::endl;
        for (int i=0; i<num_chords; i++) {
            txtoutput << labels[2*i] << "," << labels[2*i + 1] << std::endl;
        }
        txtoutput << "------" << std::endl;
        int num_chord_notes = int(outbuff[0][0]);
        int j=0;
        for (int i=0; i<num_chord_notes; i++) {
            txtoutput << outbuff[0][j+1] << "," << outbuff[0][j+2] << "," << outbuff[0][j+3] << std::endl;
            j+=3;
        }
    }
    else if (model_name == "hummingest") {
        float** fullLengthbuff;
        fullLengthbuff = new float* [num_channels];
        for (int i=0; i<num_channels; i++) {
            fullLengthbuff[i] = new float[file_length];
        }

        // read the whole file
        int num_samples = input->read(fullLengthbuff, file_length);

        // process loudness, and process tuning
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = std::min(block_size, file_length-i);

            // copy to buff
            for (int c = 0; c < num_channels; c++) {
                int bytesToCopy = num_samples * sizeof(float);
                memcpy(buff[c], &fullLengthbuff[c][i], bytesToCopy);
            }
            if (num_samples < block_size) {
                // pad remaining buf with zeros
                for (int j=0; j<num_channels; j++) {
                    for (int k=num_samples; k<block_size; k++) {
                        buff[j][k] = 0;
                    }
                }
            }
            m_chromatuning->processBlock(buff, num_samples);
            m_loudnessmeter->processBlock(buff, num_samples); // process loudness meter too
        }

        // autogain adjustment
        float dbloudness = m_loudnessmeter->getScalarMeasurement();
        std::cerr << "dbloudness:" << dbloudness << std::endl;
        float target_LUFS = -20;
        float dbMakeUp = target_LUFS - dbloudness;
        std::cerr << "dbMakeUp:" << dbMakeUp << std::endl;

        if (dbMakeUp > 6) {
            m_limiter->setThreshold(-1);
            m_limiter->setMakeUpGain(dbMakeUp);
            for (int i = 0; i < file_length; i+=block_size) {
                int num_samples = std::min(block_size, file_length-i);

                // copy to buff
                for (int c = 0; c < num_channels; c++) {
                    int bytesToCopy = num_samples * sizeof(float);
                    memcpy(buff[c], &fullLengthbuff[c][i], bytesToCopy);
                }
                if (num_samples < block_size) {
                    // pad remaining buf with zeros
                    for (int j=0; j<num_channels; j++) {
                        for (int k=num_samples; k<block_size; k++) {
                            buff[j][k] = 0;
                        }
                    }
                }
                m_limiter->processBlock(buff, num_samples);
                // copy back to original array
                for (int c = 0; c < num_channels; c++) {
                    int bytesToCopy = num_samples * sizeof(float);
                    memcpy(&fullLengthbuff[c][i], buff[c], bytesToCopy);
                }
            }
        }

        // process pYIN
        m_pyin->processInData(fullLengthbuff, num_samples);
        m_pyin_2->processInData(fullLengthbuff, num_samples);

        // write tuning results first
        float globaltuning = m_chromatuning->getScalarMeasurement();
        
        // then pyin results
        int numNotes = m_pyin->getScalarMeasurement();
        int numNotes_2 = m_pyin_2->getScalarMeasurement();
        int pYIN_selection = 0;
        std::cerr << "numNotes:" << numNotes << ", numNotes_2:" << numNotes_2 << std::endl;
        if (numNotes >= numNotes_2 - extra_meta_param) { // take the one that outputs more notes, and prefer pyin_1 always
            m_modbase_analyzer = m_pyin.get();
            pYIN_selection = 1;
            std::cerr << "use pYIN **********:" << std::endl;
        } else {
            m_modbase_analyzer = m_pyin_2.get();
            pYIN_selection = 2;
            numNotes = numNotes_2;
            std::cerr << "use pYIN 2 **********:" << std::endl;
        }

        txtoutput << 0 << "," << 0 << "," << globaltuning << "," << pYIN_selection << "," << dbMakeUp << std::endl;

        if (numNotes < block_size * 4) {
            m_modbase_analyzer->getOutData(outbuff, numNotes);
            for (int i=0; i<numNotes*3; i+=3) {
                float st = outbuff[0][i];
                float dur = outbuff[0][i+1];
                float f0 = outbuff[0][i+2];
                if (txtoutput.is_open()) {
                    // std::cerr << "txtoutput is opened" << std::endl;
                    txtoutput << st << "," << dur << "," << f0 << std::endl;
                } else {
                    std::cerr << "txtoutput not opened" << std::endl;
                }
                
            }
        } else {
            std::cerr << "numNotes > block_size * 4, this is strange!" << std::endl;
        }

        for (int i=0; i<num_channels; i++) {
            delete[] fullLengthbuff[i];
        }
        delete[] fullLengthbuff;
    }
    else if (model_name == "autogain") {
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            m_loudnessmeter->processBlock(buff, num_samples);
        }
        float dbloudness = m_loudnessmeter->getScalarMeasurement();
        std::cerr << "dbloudness:" << dbloudness << std::endl;
        m_limiter->setThreshold(-1);
        // float dbMakeUp = -dbloudness-20 > 0 ? -dbloudness-20 : 0;
        // float dbMakeUp = -dbloudness-15; // 15 here is a fixed number (-15 LUFS target)
        float dbMakeUp = target_LUFS - dbloudness;
        std::cerr << "dbMakeUp:" << dbMakeUp << std::endl;
        m_limiter->setMakeUpGain(dbMakeUp);
        input->rewind();
        for (int i = 0; i < file_length; i+=block_size) {
            int num_samples = input->read(buff, block_size);
            m_limiter->processBlock(buff, num_samples);
            if (m_limiter->outputReady()) {
                output->write(buff, num_samples);
            }
        }
    }
    else { // everything else - real time application
        for (int i = 0; i < file_length; i+=block_size) {
            
            int num_samples = input->read(buff, block_size);
            m_modbase->processBlock(buff, num_samples);
            if (m_modbase->outputReady()) {
                output->write(buff, num_samples);
            } else {
                printf("output not ready block index:%d\n",i / block_size);
            }
        }
    }
    
    m_modbase = nullptr;
    m_modbase_offline = nullptr;
    m_modbase_meter = nullptr;

    for (int i=0; i<num_channels; i++) {
        delete[] buff[i];
        delete[] outbuff[i];
    }
    delete[] buff;
    delete[] outbuff;

    buff = nullptr;
    outbuff = nullptr;

    if (input != nullptr) {
        delete input;
        input = nullptr;
    }

    if (output != nullptr) {
        delete output;
        output = nullptr;
    }

    if (input_file_name == "-") {
        fclose(stdin);
    }

    if (txtoutput.is_open()) {
        // txtoutput.flush();
        txtoutput.close();
        // std::cerr << "txtoutput.close()..." << std::endl;
    }

    if (txtinput.is_open()) {
        txtinput.close();
    }

    // printf("done...\n");
    std::cerr << "done..." << std::endl;
    
    return 0;
}

