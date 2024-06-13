#ifndef SZ_GENERAL_COMPRESSOR_HPP
#define SZ_GENERAL_COMPRESSOR_HPP

#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/frontend/Frontend.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include <cstring>
#include <chrono>

size_t arrSum(size_t *arr, size_t n) {
    size_t sum = 0; 
    for (size_t i = 0; i < n; i++)
        sum += arr[i];
    return sum;
}

namespace SZ {
    template<class T, uint N, class Frontend, class Encoder, class Lossless>
    class SZGeneralCompressor : public concepts::CompressorInterface<T> {
    public:


        SZGeneralCompressor(Frontend frontend, Encoder encoder, Lossless lossless) :
                frontend(frontend), encoder(encoder), lossless(lossless) {
            static_assert(std::is_base_of<concepts::FrontendInterface<T, N>, Frontend>::value,
                          "must implement the frontend interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {
            Timer timer(true);
            size_t blkSize = conf.blkSize;
            size_t szBlk = conf.blockSize;
            size_t estRt = int(std::cbrt(conf.totalLength));
            int* metaTemp = new int[100 * ((estRt - 1) / 6 + 1) * ((estRt - 1) / 6 + 1) * ((estRt - 1) / 6 + 1)];
            size_t regCnt;
            size_t regCntSum = 0;
            std::vector<int> quant_inds;
            std::vector<int> meta_inds;
            size_t bufferSize = 393058300;
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;
            uchar *buffer_temp = buffer;
            size_t blkSize3D = blkSize*blkSize*blkSize;
            size_t blkNum = conf.totalLength/(blkSize3D);

            auto *data_pos = data;
            for (size_t i = 0; i < blkNum; ++i)
            {
                auto quant_blk = frontend.compress(data_pos, blkSize, blkSize, blkSize, szBlk);
                data_pos += blkSize3D;
                quant_inds.insert(quant_inds.end(), quant_blk.begin(), quant_blk.end());
                int* meta_blk = frontend.save(buffer_pos, regCnt, szBlk);
                for (size_t k = 0; k < regCnt; k++)
                {
                    metaTemp[regCntSum] = meta_blk[k];
                    regCntSum++;
                } 
            }
            
            // std::cout << "regCntSum: " << regCntSum << std::endl;

            if (regCntSum) {
                HuffmanEncoder<int> reg_huffman = HuffmanEncoder<int>();
                reg_huffman.preprocess_encode(metaTemp, regCntSum, 0);
                // buffer_temp = buffer_pos;
                reg_huffman.save(buffer_pos);
                // std::cout << "metatree size: " << buffer_pos - buffer_temp << std::endl;
                // buffer_temp = buffer_pos;
                reg_huffman.encode(metaTemp, regCntSum, buffer_pos);
                // std::cout << "meta size: " << buffer_pos - buffer_temp << std::endl;
                reg_huffman.postprocess_encode();
                
            }
            delete [] metaTemp;



            encoder.preprocess_encode(quant_inds, 0);
            // buffer_temp = buffer_pos;
            encoder.save(buffer_pos);
            // std::cout << "datatree size: " << buffer_pos - buffer_temp << std::endl;
            // std::cout << "buffer_pos: " << buffer_pos << std::endl;

            // buffer_temp = buffer_pos;
            encoder.encode(quant_inds, buffer_pos);
            // std::cout << "data size: " << buffer_pos - buffer_temp << std::endl;
            encoder.postprocess_encode();

            assert(buffer_pos - buffer < bufferSize);
            // std::cout << "real bufferSize: " << buffer_pos - buffer << std::endl;

            uchar *lossless_data = lossless.compress(buffer, buffer_pos - buffer, compressed_size);
            lossless.postcompress_data(buffer);

            return lossless_data;

        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;
            size_t init_length = cmpSize;
            size_t regCnt;
            size_t regCntSum = 0;
            //std::cout << "cmpSize: " << cmpSize << std::endl;

            Timer timer(true);
            size_t blkSize = frontend.get_blkSize();
            //std::cout << "blkSize: " << blkSize << std::endl;
            //std::cout << "frontend.get_num_elements()" << frontend.get_num_elements() << std::endl;
            size_t blkSize3D = blkSize*blkSize*blkSize;
            size_t blkNum = frontend.get_num_elements()/(blkSize3D);
            // std::cout << "blkNum: " << blkNum << std::endl;

            auto compressed_data = lossless.decompress(cmpData, remaining_length);
            uchar const *compressed_data_pos = compressed_data;

            uchar const *ptr[blkNum];
            size_t remainArr[blkNum];
            size_t regArrTemp[blkNum];
            size_t count = 0;
            std::vector<int> temp;

            for (size_t j = 0; j < blkNum; j++){
                ptr[count] = compressed_data_pos;
                remainArr[count] = remaining_length;
                frontend.load(compressed_data_pos, remaining_length, 1, 1, 1, regCnt, temp, 0);
                regArrTemp[count] = regCnt;
                regCntSum += regCnt;
                count++;
            }

            std::vector<int> reg_vector ;
            if (regCntSum) {
                HuffmanEncoder<int> reg_huffman = HuffmanEncoder<int>();
                reg_huffman.load(compressed_data_pos, remaining_length);
                reg_vector = reg_huffman.decode(compressed_data_pos, regCntSum);
                reg_huffman.postprocess_decode();
            }

            size_t *regArr = new size_t[count]{0};
            for(int i = 0; i < count; ++i) {
                regArr[i] = regArrTemp[i];
            }

            encoder.load(compressed_data_pos, remaining_length);

            timer.start();
            auto startDequant= std::chrono::high_resolution_clock::now();
            auto quant_inds = encoder.decode(compressed_data_pos, frontend.get_num_elements()/4);
            encoder.postprocess_decode();
            // std::cout << "de quant size: " << quant_inds.size() << std::endl;
            auto endDequant = std::chrono::high_resolution_clock::now();
            double time_takenDequant = std::chrono::duration_cast<std::chrono::nanoseconds>(endDequant - startDequant).count() * 1e-9; // Convert nanoseconds to seconds
            std::cout << "Time taken by Decode is : "  << std::fixed << std::setprecision(5) << time_takenDequant << " sec" << std::endl;
            timer.stop("Decoder");

            count = 0;

            auto *decDataPos = decData;
            size_t tempInd = 0;

            auto startDepre= std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < blkNum/4; ++i) {
                if (regCntSum) {
                    std::vector<int> sub_reg {reg_vector.begin() + arrSum(regArr, count), reg_vector.begin() + arrSum(regArr, count + 1)};
                    frontend.load(ptr[count], remainArr[count], blkSize, blkSize, blkSize, regCnt, sub_reg, 1);
                } else {
                    frontend.load(ptr[count], remainArr[count], blkSize, blkSize, blkSize, regCnt, reg_vector, 1);
                }
                std::vector<int> sub_quant {quant_inds.begin() + tempInd, quant_inds.begin() + tempInd + blkSize3D};
                tempInd += blkSize3D;
                // encoder.postprocess_decode();
                frontend.decompress(sub_quant, decDataPos);
                decDataPos += blkSize3D;
                count++;
            }
            auto endDepre = std::chrono::high_resolution_clock::now();
            double time_takenDepre = std::chrono::duration_cast<std::chrono::nanoseconds>(endDepre - startDepre).count() * 1e-9; // Convert nanoseconds to seconds
            std::cout << "Time taken by Depre is : "  << std::fixed << std::setprecision(5) << time_takenDepre << " sec" << std::endl;

            // std::vector<int> sub_reg, sub_quant(blkSize3D);

            // for (size_t i = 0; i < blkNum; ++i) {
            //     if (regCntSum) {
            //         int start = arrSum(regArr, count);
            //         int end = arrSum(regArr, count + 1);
            //         sub_reg.assign(reg_vector.begin() + start, reg_vector.begin() + end);
            //         frontend.load(ptr[count], remainArr[count], blkSize, blkSize, blkSize, regCnt, sub_reg, 1);
            //     } else {
            //         frontend.load(ptr[count], remainArr[count], blkSize, blkSize, blkSize, regCnt, reg_vector, 1);
            //     }
                
            //     std::copy(quant_inds.begin() + tempInd, quant_inds.begin() + tempInd + blkSize3D, sub_quant.begin());
            //     tempInd += blkSize3D;

            //     frontend.decompress(sub_quant, decDataPos);
            //     decDataPos += blkSize3D;
            //     count++;
            // }

            // std::cout << "down decomp" << std::endl;

            lossless.postdecompress_data(compressed_data);
            return decData;
        }


    private:
        Frontend frontend;
        Encoder encoder;
        Encoder encoder0;
        Lossless lossless;
    };

    template<class T, uint N, class Frontend, class Encoder, class Lossless>
    std::shared_ptr<SZGeneralCompressor<T, N, Frontend, Encoder, Lossless>>
    make_sz_general_compressor(Frontend frontend, Encoder encoder, Lossless lossless) {
        return std::make_shared<SZGeneralCompressor<T, N, Frontend, Encoder, Lossless>>(frontend, encoder, lossless);
    }


}
#endif

