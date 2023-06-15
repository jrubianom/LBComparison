#include <iostream>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

void createDir(const std::string& filePath) {
    // Extract the directory path from the file path
    fs::path directoryPath = fs::path(filePath).parent_path();

    // Create the directory if it doesn't exist
    if (!fs::exists(directoryPath)) {
        fs::create_directories(directoryPath);
    }
}
