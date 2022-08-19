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

    if (!exists_test(input_file_name) && input_file_name != "-") {
        std::cerr << "input_file doesn't exist" << std::endl;
        return -1;
    }
    
    //create input stream and reader
    WavInFile *input = nullptr;
    
    if (input_file_name != "-") {
        input = new WavInFile(input_file_name.c_str());
    } else {
        FILE *fp = std::freopen(nullptr, "rb", stdin);

        if (std::ferror(stdin))
            throw std::runtime_error(std::strerror(errno));

        std::cerr << "reading from stdin..." << std::endl;

        input = new WavInFile(stdin);
    }

    int wav_format = input->getWavFormat();
    std::cerr << "wav_format = " << wav_format << std::endl;

    int bytes_per_sample = input->getBytesPerSample();
    std::cerr << "bytes_per_sample = " << bytes_per_sample << std::endl;

    int data_len_in_bytes = input->getDataSizeInBytes();
    std::cerr << "data_len_in_bytes = " << data_len_in_bytes << std::endl;
    
    // num_channels = reader->numChannels;
    int num_channels = input->getNumChannels();
    std::cerr << "num_channels = " << num_channels << std::endl;

    // sample_rate = reader->sampleRate;
    int sample_rate = input->getSampleRate();
    std::cerr << "sample_rate = " <<  sample_rate << std::endl;

    int file_length = input->getNumSamples();
    std::cerr << "file_length = " << file_length << std::endl;

    int output_bits_per_sample = 16;

    int target_LUFS = -15;

    WavOutFile *output = nullptr;
    std::ofstream txtoutput;
    if (model_name == "loudnessmeter" || model_name == "envelope") {
        // txtoutput.open(output_file_name);
        // use std::cout as output
    } else {
        output = new WavOutFile(output_file_name.c_str(), sample_rate, output_bits_per_sample, num_channels);
    }

    int block_size = sample_rate / 100 < 480 ? 480 : sample_rate / 100;
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
        // txtoutput << m_modbase_meter->getScalarMeasurement() << std::endl;
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

    delete input;

    if (output != nullptr) {
        delete output;
        output = nullptr;
    }

    if (input_file_name == "-") {
        fclose(stdin);
    }

    if (txtoutput.is_open()) {
        txtoutput.close();
    }

    // printf("done...\n");
    // std::cerr << "done..." << std::endl;
    
    return 0;
}

