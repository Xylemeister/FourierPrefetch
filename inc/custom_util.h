#ifndef CUSTOM_H
#define CUSTOM_H


template <class T>
class SetAssociativeCache {
public:
    class Entry {
    public:
        uint64_t key;
        uint64_t index;
        uint64_t tag;
        bool valid;
        T data;
    };

    SetAssociativeCache(int size, int num_ways, int debug_level = 0) :
        size(size), num_ways(num_ways), num_sets(size / num_ways), entries(num_sets, std::vector<Entry>(num_ways)),
        cams(num_sets, std::unordered_map<uint64_t, int>(num_ways)), debug_level(debug_level) {
        // assert(size % num_ways == 0);
        for (int i = 0; i < num_sets; i += 1)
            for (int j = 0; j < num_ways; j += 1)
                entries[i][j].valid = false;
        /* calculate `index_len` (number of bits required to store the index) */
        for (int max_index = num_sets - 1; max_index > 0; max_index >>= 1)
            this->index_len += 1;
    }

    /**
     * Invalidates the entry corresponding to the given key.
     * @return A pointer to the invalidated entry
     */
    Entry* erase(uint64_t key) {
        Entry* entry = this->find(key);
        uint64_t index = key & ((1 << this->index_len) - 1);
        uint64_t tag = key >> this->index_len;
        auto& cam = cams[index];
        int num_erased = cam.erase(tag);
        if (entry)
            entry->valid = false;
        // assert(entry ? num_erased == 1 : num_erased == 0);
        return entry;
    }

    /**
     * @return The old state of the entry that was updated
     */
    Entry insert(uint64_t key, const T& data) {
        Entry* entry = this->find(key);
        if (entry != nullptr) {
            Entry old_entry = *entry;
            entry->data = data;
            return old_entry;
        }
        uint64_t index = key & ((1 << this->index_len) - 1);
        uint64_t tag = key >> this->index_len;
        std::vector<Entry>& set = this->entries[index];
        int victim_way = -1;
        for (int i = 0; i < this->num_ways; i += 1)
            if (!set[i].valid) {
                victim_way = i;
                break;
            }
        if (victim_way == -1) {
            victim_way = this->select_victim(index);
        }
        Entry& victim = set[victim_way];
        Entry old_entry = victim;
        victim = {key, index, tag, true, data};
        auto& cam = cams[index];
        if (old_entry.valid) {
            int num_erased = cam.erase(old_entry.tag);
            // assert(num_erased == 1);
        }
        cam[tag] = victim_way;
        return old_entry;
    }

    Entry* find(uint64_t key) {
        uint64_t index = key & ((1 << this->index_len) - 1);
        uint64_t tag = key >> this->index_len;
        auto& cam = cams[index];
        if (cam.find(tag) == cam.end())
            return nullptr;
        int way = cam[tag];
        Entry& entry = this->entries[index][way];
        // assert(entry.tag == tag && entry.valid);
        if (!entry.valid)
            return nullptr;
        return &entry;
    }

    void flush() {
        for (int i = 0; i < num_sets; i += 1) {
            cams[i].clear();
            for (int j = 0; j < num_ways; j += 1)
                entries[i][j].valid = false;
        }
    }

    /**
     * Creates a table with the given headers and populates the rows by calling `write_data` on all
     * valid entries contained in the cache. This function makes it easy to visualize the contents
     * of a cache.
     * @return The constructed table as a std::string
     */
    std::string log(std::vector<std::string> headers) {
        std::vector<Entry> valid_entries = this->get_valid_entries();
        Table table(headers.size(), valid_entries.size() + 1);
        table.set_row(0, headers);
        for (unsigned i = 0; i < valid_entries.size(); i += 1)
            this->write_data(valid_entries[i], table, i + 1);
        return table.to_string();
    }

    int get_index_len() { return this->index_len; }

    void set_debug_level(int debug_level) { this->debug_level = debug_level; }

protected:
    /* should be overriden in children */
    virtual void write_data(Entry& entry, Table& table, int row) {}

    /**
     * @return The way of the selected victim
     */
    virtual int select_victim(uint64_t index) {
        /* random eviction policy if not overriden */
        return rand() % this->num_ways;
    }

    std::vector<Entry> get_valid_entries() {
        std::vector<Entry> valid_entries;
        for (int i = 0; i < num_sets; i += 1)
            for (int j = 0; j < num_ways; j += 1)
                if (entries[i][j].valid)
                    valid_entries.push_back(entries[i][j]);
        return valid_entries;
    }

    int size;
    int num_ways;
    int num_sets;
    int index_len = 0; /* in bits */
    std::vector<std::vector<Entry>> entries;
    std::vector<std::unordered_map<uint64_t, int>> cams;
    int debug_level = 0;
};

template <class T>
class LRUSetAssociativeCache : public SetAssociativeCache<T> {
    typedef SetAssociativeCache<T> Super;

public:
    LRUSetAssociativeCache(int size, int num_ways, int debug_level = 0) :
        Super(size, num_ways, debug_level), lru(this->num_sets, std::vector<uint64_t>(num_ways)) {}

    void set_mru(uint64_t key) { *this->get_lru(key) = this->t++; }

    void set_lru(uint64_t key) { *this->get_lru(key) = 0; }

    void rp_promote(uint64_t key) { set_mru(key); }

    void rp_insert(uint64_t key) { set_mru(key); }

protected:
    /* @override */
    int select_victim(uint64_t index) {
        std::vector<uint64_t>& lru_set = this->lru[index];
        return min_element(lru_set.begin(), lru_set.end()) - lru_set.begin();
    }

    uint64_t* get_lru(uint64_t key) {
        uint64_t index = key % this->num_sets;
        uint64_t tag = key / this->num_sets;
        // assert(this->cams[index].count(tag) == 1);
        int way = this->cams[index][tag];
        return &this->lru[index][way];
    }

    std::vector<std::vector<uint64_t>> lru;
    uint64_t t = 1;
};

template <class T>
class SRRIPSetAssociativeCache : public SetAssociativeCache<T> {
    typedef SetAssociativeCache<T> Super;

public:
    SRRIPSetAssociativeCache(int size, int num_ways, int debug_level = 0, int max_rrpv = 3) :
        Super(size, num_ways, debug_level), rrpv(this->num_sets, std::vector<uint64_t>(num_ways)),
        max_rrpv(max_rrpv) {}

    void rp_promote(uint64_t key) { *this->get_rrpv(key) = 0; }

    void rp_insert(uint64_t key) { *this->get_rrpv(key) = 2; }

protected:
    /* @override */
    int select_victim(uint64_t index) {
        std::vector<uint64_t>& rrpv_set = this->rrpv[index];
        for (;;) {
            for (int i = 0; i < this->num_ways; i++) {
                if (rrpv_set[i] >= max_rrpv) {
                    return i;
                }
            }
            aging(index);
        }
    }

    uint64_t* get_rrpv(uint64_t key) {
        uint64_t index = key % this->num_sets;
        uint64_t tag = key / this->num_sets;
        // assert(this->cams[index].count(tag) == 1);
        int way = this->cams[index][tag];
        return &this->rrpv[index][way];
    }

private:
    void aging(uint64_t index) {
        std::vector<uint64_t>& rrpv_set = this->rrpv[index];
        for (auto& r : rrpv_set) {
            ADD(r, max_rrpv);
        }
    }

    std::vector<std::vector<uint64_t>> rrpv;
    int max_rrpv;
};






#endif
