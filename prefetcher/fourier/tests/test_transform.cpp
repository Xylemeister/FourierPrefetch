
#include <cstdint>
#include <vector>
#include <gtest/gtest.h>
#include "../transform.h"


// --- Test the buffer ----
using namespace Transform;


TEST(Transform, BufferTest){
    std::vector<uint64_t> expected = {420, 69, 0, 67};

    auto* fut = new TransformBuf<3>;

    for (auto a: expected){
        fut->Insert(a);
    }

    const uint64_t* it = fut->viewBuf();

    EXPECT_EQ(it[0], 67);
    EXPECT_EQ(it[1], 69);
    EXPECT_EQ(it[2], 0);

    fut->Insert(1);

    EXPECT_EQ(it[0],67);
    EXPECT_EQ(it[1], 1);
    EXPECT_EQ(it[2], 0);


    delete fut;
}
