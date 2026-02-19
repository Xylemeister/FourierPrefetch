#ifndef FOURIER_H
#define FOURIER_H


#include <cstdint>
#include <array>


constexpr int  DELTA_BUFF_SIZE = 24;
constexpr int FREQ_BINS_SIZE = 6;
// template<typename T>
// class Buffer {
//     public:
//         Buffer(int size_): size(size_){};
//         int insert(T data);
//         int get()



//     private:
//         int size;
//         std::deque<T> buf;
// };




class Fourier {
    public:

    struct FourierEntry {
        uint64_t confidence;
        uint64_t ip_tag;
        uint64_t last_cl_addr;
        // some array for the deltas
        // the frequency bins?
        std::array<uint64_t, DELTA_BUFF_SIZE>  deltas;
        std::array<uint64_t, FREQ_BINS_SIZE> frequency;
    };

};











#endif
