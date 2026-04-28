#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
    // Command to convert MP3 to WAV using ffmpeg, outputting to stdout
    const char* command = "ffmpeg -nostdin -i test.mp3 -f wav -";
    // const char* command = "xxd test.wav";

    // Open a pipe to read the output of the ffmpeg command
    FILE* pipe = popen(command, "r");
    if (!pipe) {
        std::cerr << "Failed to open pipe for command: " << command << std::endl;
        return EXIT_FAILURE;
    }

    FILE *csvfile;
    csvfile = fopen("samples.csv", "w");
    

    // Read data from the pipe into a buffer
    std::vector<char> buffer(4096); // Adjust buffer size as needed
    size_t bytesRead = 0;
    while ((bytesRead = fread(buffer.data(), 1, buffer.size(), pipe)) > 0) {
        // Process the data read from the pipe
        // For example, write it to a file or analyze it in memory
        fwrite(buffer.data(), 1, bytesRead, csvfile); // Example of writing to a file
        // fprintf(csvfile, "%d\n", buffer[0]);
    }
    fclose(csvfile);

    // Check for errors in reading from the pipe
    if (ferror(pipe)) {
        std::cerr << "Error reading from pipe" << std::endl;
        pclose(pipe);
        return EXIT_FAILURE;
    }

    // Close the pipe
    if (pclose(pipe) == -1) {
        std::cerr << "Failed to close pipe" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


