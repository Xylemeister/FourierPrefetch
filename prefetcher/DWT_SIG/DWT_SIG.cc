#include <cmath>
#include <cstdint>
#include <vector>

class CircularBuffer {
public:
    CircularBuffer(int size) : buffer_(size, 0), size_(size) {}

    void insert(int64_t val) {
        buffer_[head_] = val;
        head_ = (head_ + 1) % size_;
        if (count_ < size_) ++count_;
    }

    int64_t at(int i) const {
        int idx = (head_ - count_ + i + size_) % size_;
        return buffer_[idx];
    }

    int64_t newest() const {
        return buffer_[(head_ - 1 + size_) % size_];
    }

    int64_t oldest() const {
        return buffer_[(head_ - count_ + size_) % size_];
    }

    int size() const { return count_; }
    int capacity() const { return size_; }
    bool full() const { return count_ == size_; }
    void clear() { head_ = 0; count_ = 0; }

private:
    std::vector<int64_t> buffer_;
    int size_ = 0;
    int head_ = 0;
    int count_ = 0;
};


class DWT {
    public:
    DWT(int size, int bsize): buff(size) {
        approx_[0].resize(size / 2);
        detail_[0].resize(size / 2);
        approx_[1].resize(size / 4);
        detail_[1].resize(size / 4);
        approx_[2].resize(size / 8);
        detail_[2].resize(size / 8);
    }

    void insert(int64_t val){
        buff.insert(val);
    }

    const std::vector<double>& approx(int level) const { return approx_[level]; }
    const std::vector<double>& detail(int level) const { return detail_[level]; }

    void transform(){
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        for (int bank = 0; bank < 3; bank++){
            if (bank == 0){
                int n = buff.capacity() / 2;
                for (int i = 0; i < n; i++){
                    double a = static_cast<double>(buff.at(2 * i));
                    double b = static_cast<double>(buff.at(2 * i + 1));
                    approx_[0][i] = (a + b) * inv_sqrt2;
                    detail_[0][i] = (a - b) * inv_sqrt2;
                }
            } else {
                int n = approx_[bank - 1].size() / 2;
                for (int i = 0; i < n; i++){
                    double a = approx_[bank - 1][2 * i];
                    double b = approx_[bank - 1][2 * i + 1];
                    approx_[bank][i] = (a + b) * inv_sqrt2;
                    detail_[bank][i] = (a - b) * inv_sqrt2;
                }
            }
        }
    }


    private:
        CircularBuffer buff;
        std::vector<double> approx_[3];
        std::vector<double> detail_[3];
};
