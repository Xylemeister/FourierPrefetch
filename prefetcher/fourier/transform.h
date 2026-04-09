#ifndef TRANSFORM_H
#define TRANSFORM_H

// This library is the definition of the fourier transforms or any other
// transform that we will do to our address accesses.

#include <cstddef>
#include <cstdint>
#include <array>
#include <functional>
#include <iterator>
#include <memory>
#include <cmath>
#include <iostream>
#include "msl/lru_table.h"



namespace Transform {

    template<std::size_t SIZE>
    class TransformBuf {
        // static constexpr uint8_t DELTA_BUFF_SIZE = 16;
        public:

        using WaveType = std::array<std::array<double, SIZE>, SIZE>;
        static constexpr std::size_t EFFECTIVE_SIZE = SIZE/2 + 1;
        static constexpr double PI = 3.14159265358979323846;

        // this should set all the entries to zeroes
        // when we started sampling we fill out the windows with zeroes for the fourier
        // maybe it is better to just do the FT on valid values? but
        // to my knowledge this is how it's usually done
        TransformBuf(){
            for (std::size_t i = 0; i < SIZE; i++){
                buf_[i] = 0;
            }
        }

        // this push will push to the buffer if buffer is full the oldest gets thrown
        void Insert(int64_t elem){
            buf_[head] = elem;
            head = (head + 1) % SIZE;
        }

        void Transform(){
            // we are using effective size by exploiting the mirror behaviour of dft with real values
            // need to test this
            for (std::size_t i = 0; i < EFFECTIVE_SIZE; i++){
                double re = 0.0, im = 0.0;
                for (std::size_t j = 0; j < SIZE; j++){
                    std::size_t offset = (head + j) % SIZE;
                    re += buf_[offset] * CosTable_[i][j];
                    im -= buf_[offset] * SineTable_[i][j];
                }

                // there are multiple ways of using this value, for now we are implementing it
                // as a power, we could look at other implementations later such as taking sqrt for mag
                // or we can take the real or imaginary directly/separately
                bin_[i] = re*re + im * im;

            }

        }


        const int64_t* viewBuf(){
            return buf_;
        }

        const auto viewBin(){
            return bin_;
        }


        static WaveType ComputeTable(double(*fn)(double)){
            WaveType result;
            for (std::size_t i = 0; i < SIZE; i++){
                for (std::size_t j = 0; j < SIZE; j++){
                    double arg = (double(i) * double(j) * PI * 2.0) / SIZE;
                    result[i][j] = fn(arg);
                }
            }
            return result;
        }

        private:
        // points at the next pos to be written
        size_t head = 0;

        /*
         *  I envision the transform to be something like for loop from head to Size, taking care
         *  of the circular behaviour of this buffer
         */
        int64_t buf_[SIZE];

        std::array<double, EFFECTIVE_SIZE> bin_;
        static inline const WaveType SineTable_ = ComputeTable(std::sin);
        static inline const WaveType CosTable_ = ComputeTable(std::cos);
    };


    class FourierPrefetchV1 {

        // for a start with start with this stupid idea of using the behaviour instead of IP
        // for initiating prefetch, my hypothesis is that different IPs might be using the same
        // data structure with the same pattern thus we can use information learned in other IPs
        // hypothesis is that the pattern/frequency behaviour is unique for different behaviours
        // for now we just use the strongest one, it can be a hash type of thing in the future
        constexpr static size_t WINDOW_SIZE = 16;

        struct tracker_entry {
            uint64_t ip = 0; // the ip tag, to differentiate regions
            uint64_t last_cl_addr = 0; // the last cl addr accessed for calculating deltas
            TransformBuf<WINDOW_SIZE> buf{}; // history of the last 16 deltas, this engine gives us the transform

            auto index() const {
                return ip;
            }

            auto tag() const {
                return ip;
            }
        };


        // struct lookahead_entry {
        //     uint64_t stride;
        // }

        constexpr static std::size_t TRACKER_SETS = 256;
        constexpr static std::size_t TRACKER_WAYS = 4;
        // constexpr static int PREFETCH_DEGREE = 3;

        champsim::msl::lru_table<tracker_entry> table{TRACKER_SETS, TRACKER_WAYS};
        public:
        void debug_lookup(uint64_t cl_addr, uint64_t ip){
            TransformBuf<WINDOW_SIZE> buf;
            auto found= table.check_hit({ip, cl_addr, buf});


            if (found.has_value()){
                auto delta = static_cast<int64_t>(cl_addr) - static_cast<int64_t>(found->last_cl_addr);
                found->last_cl_addr = cl_addr;
                table.fill(*found);

                found->buf.Insert(delta);
                found->buf.Transform();
                table.fill(*found);

                auto view = found->buf.viewBin();

                for (std::size_t i = 0; i < TransformBuf<WINDOW_SIZE>::EFFECTIVE_SIZE; i++){
                    std::cout << " " << view[i] << "";
                }
                std::cout << std::endl;

            }
            else{
                table.fill({ip, cl_addr, buf});
            }

        }


    };

    class FourierPrefetchV2 {
       // this is the second class of fourier prefetch
       // in here, we derive the PF deltas directly from the result of our fourier transform
       // such taht if we have a good strong period then we can prefetch the last "PERIOD" deltas
       // im not sure how this is going to work as the concept of fourier transform doesnt understand time
       // or where we are currently in the wave, we might be able to find some ways to circumvent that with a
       // DIscrete wavelet transform which I still need to understand
       // basically a DWT, instead of using infinite-horizon cosine to sample or find similarity in the wave it uses
       // a finite time one
       // the problem with using the raw  fourier series is that we have infinite time cosine in such a case
    };


    class FourierPrefetchV3{
      // this is another approach where we use the fourier series as an idea to make filtering decisions in our pipeline
    };

} // Transform


#endif
