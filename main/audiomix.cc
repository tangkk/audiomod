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
        std::cerr << "usage: ./audiomix-exe input1 input2 ... (order should be: melody, harmony, bass, drum)" << std::endl;
        return -1;
    }

    int st = getTimeOfDay();

    std::vector<std::string> input_file_name_vec;
    int i = 1;
    while(argv[i] != NULL) {
        std::string input_file_name(argv[i]);
        if (!exists_test(input_file_name) && input_file_name != "-") {
            std::cerr << "input_file doesn't exist" << std::endl;
            return -1;
        }

        input_file_name_vec.push_back(input_file_name);
        i++;
    }

    std::vector<WavInFile*> input_files;
    for(const auto& input_file_name : input_file_name_vec) {
        std::cout << input_file_name << std::endl;    

        WavInFile* input = new WavInFile(input_file_name.c_str());
        input_files.push_back(input);
    }

    std::string output_file_name = "out.wav";
    WavOutFile *output = nullptr;
    // fixed params from the first input
    WavInFile *input_0 = input_files[0];
    const int wav_format = input_0->getWavFormat();
    std::cerr << "wav_format = " << wav_format << std::endl;

    const int bytes_per_sample = input_0->getBytesPerSample();
    std::cerr << "bytes_per_sample = " << bytes_per_sample << std::endl;

    const int data_len_in_bytes = input_0->getDataSizeInBytes();
    std::cerr << "data_len_in_bytes = " << data_len_in_bytes << std::endl;

    const int file_length = input_0->getNumSamples();
    std::cerr << "file_length = " << file_length << std::endl;
    
    const int num_channels = input_0->getNumChannels();
    std::cerr << "num_channels = " << num_channels << std::endl;

    const int sample_rate = input_0->getSampleRate();
    std::cerr << "sample_rate = " <<  sample_rate << std::endl;

    const int block_size = sample_rate / 100 < 480 ? 480 : sample_rate / 100;
    std::cerr << "block_size = " << block_size << std::endl;
    output = new WavOutFile(output_file_name.c_str(), sample_rate, 16, 2);

    std::unique_ptr<limiter> m_limiter = std::unique_ptr<limiter>(new limiter(sample_rate, num_channels, -2, 0));
    // modbase * m_modbase = m_limiter.get();
    std::unique_ptr<loudnessmeter> m_loudnessmeter = std::unique_ptr<loudnessmeter>(new loudnessmeter(sample_rate, num_channels, block_size));
    // modbase_meter * m_modbase_meter = m_loudnessmeter.get();

    // buffer for a block
    float** buff;
    buff = new float* [num_channels];
    for (int i=0; i<num_channels; i++) {
        buff[i] = new float[block_size];
    }

    float** outbuff;
    outbuff = new float* [num_channels];
    for (int i=0; i<num_channels; i++) {
        outbuff[i] = new float[block_size * 4]; // note that this is larger, but it still could only store blocksize elems.
    }

    // get loudness profiles here if necessary

    // do the actual mixing here, go through a limiter if necessary
    std::vector<float> fixed_gain = {
        0.35 * 1.15, // melody
        0.75 * 1.15, // harmony
        0.55 * 1.15, // bass
        0.95 * 1.15 // drum
    };
    for (int i = 0; i < file_length; i+=block_size) { // mix blockwise
        int num_sample_out = 0;
        for (int j = 0; j < input_files.size(); j++) {
            WavInFile *input = input_files[j];
            int num_samples = input->read(buff, block_size);

            float gain = j < 4 ? fixed_gain[j] : 1;

            if (j == 0) {
                for (int c=0; c<num_channels; c++) {
                    for (int k=0; k<num_samples; k++) {
                        outbuff[c][k] = buff[c][k] * gain;
                    }
                }
            } else {
                for (int c=0; c<num_channels; c++) {
                    for (int k=0; k<num_samples; k++) {
                        outbuff[c][k] += buff[c][k] * gain;
                    }
                }
            }

            // process outbuff if input_files all done
            num_sample_out = num_samples;
            if (j == input_files.size() - 1) {
                m_limiter->processBlock(outbuff, num_sample_out);
            }
        }
        output->write(outbuff, num_sample_out);
    }

    for (WavInFile *input : input_files) {
        delete input;
        input = nullptr;
    }
    input_files.clear();
    // m_modbase = nullptr;
    // m_modbase_meter = nullptr;

    if (output != nullptr) {
        delete output;
        output = nullptr;
    }

    for (int i=0; i<num_channels; i++) {
        delete[] buff[i];
        delete[] outbuff[i];
    }
    delete[] buff;
    delete[] outbuff;

    buff = nullptr;
    outbuff = nullptr;

    int et = getTimeOfDay();
    std::cerr << "time est in ms:" << (et - st) << std::endl;
}