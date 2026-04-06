#ifndef TRANSFORM_H
#define TRANSFORM_H

// This library is the definition of the fourier transforms or any other
// transform that we will do to our address accesses.

#include <cstddef>
#include <cstdint>
#include <array>


namespace Transform {

    template<std::size_t Size>
    class TransformBuf {
        // static constexpr uint8_t DELTA_BUFF_SIZE = 16;
        public:

        // this should set all the entries to zeroes
        // when we started sampling we fill out the windows with zeroes for the fourier
        // maybe it is better to just do the FT on valid values? but
        // to my knowledge this is how it's usually done
        TransformBuf(){
            for (int i = 0; i < Size; i++){
                buf[i] = 0;
            }
        }

        // this push will push to the buffer if buffer is full the oldest gets thrown
        void Insert(uint64_t elem){
            buf[head] = elem;
            head = (head + 1) % Size;
        }

        const float* Transform();


        const uint64_t* viewBuf(){
            return buf;
        }

        private:
        // points at the next pos to be written
        size_t head = 0;

        /*
         *  I envision the transform to be something like for loop from head to Size, taking care
         *  of the circular behaviour of this buffer
         */
        uint64_t buf[Size];
        float bin[Size];

    };

} // Transform


#endif
