#ifndef TRANSFORM_H
#define TRANSFORM_H

// This library is the definition of the fourier transforms or any other
// transform that we will do to our address accesses.

#include <cstddef>
#include <cstdint>
#include <array>
#include <functional>
#include <memory>
#include <cmath>



namespace Transform {

    template<std::size_t SIZE>
    class TransformBuf {
        // static constexpr uint8_t DELTA_BUFF_SIZE = 16;
        public:

        using WaveType = std::array<std::array<double, SIZE>, SIZE>;
        static const int EFFECTIVE_SIZE = SIZE/2 + 1;

        // this should set all the entries to zeroes
        // when we started sampling we fill out the windows with zeroes for the fourier
        // maybe it is better to just do the FT on valid values? but
        // to my knowledge this is how it's usually done
        TransformBuf(){
            for (int i = 0; i < SIZE; i++){
                buf_[i] = 0;
            }
        }

        // this push will push to the buffer if buffer is full the oldest gets thrown
        void Insert(uint64_t elem){
            buf_[head] = elem;
            head = (head + 1) % SIZE;
        }

        void Transform(){
            // we are using effective size by exploiting the mirror behaviour of dft with real values
            // need to test this
            for (int i = 0; i < EFFECTIVE_SIZE; i++){
                double re = 0.0, im = 0.0;
                for (int j =0; j < SIZE; j++){
                    int offset = (head + j) % SIZE;
                    re += buf_[offset] * CosTable_[i][j];
                    im -= buf_[offset] * SineTable_[i][j];
                }

                // there are multiple ways of using this value, for now we are implementing it
                // as a power, we could look at other implementations later such as taking sqrt for mag
                // or we can take the real or imaginary directly/separately
                bin_[i] = re*re + im * im;

            }

        }


        const uint64_t* viewBuf(){
            return buf_;
        }

        const double* viewBin(){
            return bin_;
        }


        // i dont think this is even conexpritable buit it is what it is
        constexpr static WaveType ComputeSine(){
            WaveType sin;
            for (int i = 0; i < SIZE; i++){
                for (int j = 0; j < SIZE; j++){
                    double arg = (double(i) * double(j)  * M_PI * 2.0) / SIZE;
                    sin[i][j] = std::sin(arg);

                }
            }

            return sin;
        }

        constexpr static WaveType ComputeCos(){
            WaveType cos;
            for (int i = 0; i < SIZE; i++){
                for (int j = 0; j < SIZE; j++){
                    double arg = (double(i) * double(j)  * M_PI * 2.0) / SIZE;
                    cos[i][j] = std::cos(arg);

                }
            }

            return cos;
        }

        private:
        // points at the next pos to be written
        size_t head = 0;

        /*
         *  I envision the transform to be something like for loop from head to Size, taking care
         *  of the circular behaviour of this buffer
         */
        uint64_t buf_[SIZE];

        std::array<double, EFFECTIVE_SIZE> bin_;
        constexpr static WaveType SineTable_ = ComputeSine();
        constexpr static WaveType CosTable_ = ComputeCos();


    };

} // Transform


#endif
