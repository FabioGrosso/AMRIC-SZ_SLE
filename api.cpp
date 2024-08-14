#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "sz.hpp"

using namespace std;
const double errorBound = 1.0;

const int smallX = 16;
const int bigX = 32;
const int smallY = 16;
const int bigY = 32;
const int smallZ = 24;
const int bigZ = 32;

const int SIZE_X = 512;
const int SIZE_Y = 512;
const int SIZE_Z = 512;

const int BLKSIZE = 16;
const int TOTAL_SIZE = SIZE_X * SIZE_Y * SIZE_Z;

const char* inputFile = "SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";

using DataArray = float*;

// Function to write the box dimensions to a file
void writeBoxFile(const string& fileName, int smallX, int smallY, int smallZ, int bigX, int bigY, int bigZ) {
    ofstream outFile(fileName);
    if (outFile.is_open()) {
        outFile << smallX << " " << smallY << " " << smallZ << " ";
        outFile << bigX << " " << bigY << " " << bigZ;
        outFile.close();
    } else {
        cerr << "Unable to open file for box writing\n";
        exit(1);
    }
}

// Function to read binary data from a file
DataArray readBinaryData(const string& filePath) {
    DataArray data = new float[TOTAL_SIZE];
    ifstream inFile(filePath, ios::binary);

    if (!inFile) {
        cerr << "Unable to open file for reading: " << filePath << endl;
        exit(1);
    }

    inFile.read(reinterpret_cast<char*>(data), TOTAL_SIZE * sizeof(float));
    inFile.close();
    return data;
}

// Function to save binary data to a file
void saveBinaryData(const string& filePath, const DataArray& data) {
    ofstream outFile(filePath, ios::binary);
    if (!outFile.is_open()) {
        cerr << "Failed to open the output file for writing." << endl;
        return;
    }

    outFile.write(reinterpret_cast<const char*>(data), (bigX-smallX)*(bigY-smallY)*(bigZ-smallZ) * BLKSIZE * BLKSIZE * BLKSIZE * sizeof(float));
    outFile.close();
}

// Function to process blocks of data
DataArray processBlocks(DataArray tempData) {
    DataArray oriData = new float[TOTAL_SIZE];
    int index = 0;

    for (size_t k = 0; k < SIZE_Z / BLKSIZE; ++k) {
        for (size_t j = 0; j < SIZE_Y / BLKSIZE; ++j) {
            for (size_t i = 0; i < SIZE_X / BLKSIZE; ++i) {
                for (size_t z = k * BLKSIZE; z < k * BLKSIZE + BLKSIZE; ++z) {
                    for (size_t y = j * BLKSIZE; y < j * BLKSIZE + BLKSIZE; ++y) {
                        for (size_t x = i * BLKSIZE; x < i * BLKSIZE + BLKSIZE; ++x) {
                            oriData[index++] = tempData[x + y * SIZE_X + z * SIZE_X * SIZE_Y];
                        }
                    }
                }
            }
        }
    }

    return oriData;
}

// Function to compress and decompress data
float* AMR_compress(float* oriData, size_t blksize_x, size_t blksize_y, size_t blksize_z, double eb) {
    SZ::Config conf(blksize_z, blksize_y, blksize_x); 
    conf.cmprAlgo = SZ::ALGO_LORENZO_REG;
    conf.lorenzo = true; // only use 1st order lorenzo
    conf.lorenzo2 = false;
    conf.regression = true;
    conf.regression2 = false;
    conf.errorBoundMode = SZ::EB_ABS; 
    conf.absErrorBound = eb; 
    conf.blockSize = 4;
    conf.blkSize = BLKSIZE;
    conf.totalLength = TOTAL_SIZE;

    size_t outSize = 0;
    auto startComp = chrono::high_resolution_clock::now();
    char* compressedData = SZ_compress<float>(conf, oriData, outSize);
    auto endComp = chrono::high_resolution_clock::now();
    double time_takenComp = chrono::duration_cast<chrono::nanoseconds>(endComp - startComp).count() * 1e-9; // Convert nanoseconds to seconds
    cout << "CR: " << (SIZE_X * SIZE_Y * SIZE_Z * 4) / outSize << endl;
    cout << "Time taken by Comp is : "  << fixed << setprecision(5) << time_takenComp << " sec" << endl;

    auto startDecomp = chrono::high_resolution_clock::now();
    float* deData = new float[(bigX-smallX)*(bigY-smallY)*(bigZ-smallZ) * BLKSIZE * BLKSIZE * BLKSIZE];
    SZ_decompress<float>(conf, compressedData, outSize, deData);
    auto endDecomp = chrono::high_resolution_clock::now();
    double time_takenDecomp = chrono::duration_cast<chrono::nanoseconds>(endDecomp - startDecomp).count() * 1e-9; // Convert nanoseconds to seconds
    cout << "Time taken by Decomp is : "  << fixed << setprecision(5) << time_takenDecomp << " sec" << endl;

    delete[] compressedData;
    return deData;
}

// Function to reconstruct the decompressed data
DataArray processBack(DataArray data) {
    DataArray dataBack = new float[(bigX-smallX)*(bigY-smallY)*(bigZ-smallZ) * BLKSIZE * BLKSIZE * BLKSIZE];
    int index = 0;

    for (size_t k = 0; k < bigZ-smallZ; ++k) {
        for (size_t j = 0; j < bigY-smallY; ++j) {
            for (size_t i = 0; i < bigX-smallX; ++i) {
                for (size_t z = k * BLKSIZE; z < k * BLKSIZE + BLKSIZE; ++z) {
                    for (size_t y = j * BLKSIZE; y < j * BLKSIZE + BLKSIZE; ++y) {
                        for (size_t x = i * BLKSIZE; x < i * BLKSIZE + BLKSIZE; ++x) {
                            dataBack[x + y * (bigX-smallX) * BLKSIZE + z * (bigX-smallX) * BLKSIZE * (bigY-smallY) * BLKSIZE] = data[index++];
                        }
                    }
                }
            }
        }
    }

    return dataBack;
}

// High-level function to compress, decompress, and save data
void compressAndSave(const string& inputFilePath, const string& outputFilePath, double errorBound, int smallX, int smallY, int smallZ, int bigX, int bigY, int bigZ) {
    writeBoxFile("box.txt", smallX, smallY, smallZ, bigX, bigY, bigZ);
    DataArray tempData = readBinaryData(inputFilePath);
    DataArray oriData = processBlocks(tempData);

    float* compressedData = AMR_compress(oriData, BLKSIZE, BLKSIZE, TOTAL_SIZE / (BLKSIZE * BLKSIZE), errorBound);
    delete[] oriData;  

    DataArray dataBack = processBack(compressedData);
    saveBinaryData(outputFilePath, dataBack);

    delete[] compressedData;
    delete[] tempData;
    delete[] dataBack;
}

int main() {

    compressAndSave(inputFile, "comp.raw",
    errorBound, smallX, smallY, smallZ, bigX, bigY, bigZ);

    return 0;
}
