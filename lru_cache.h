#ifndef LRU_CACHE_H
#define LRU_CACHE_H
#include <unordered_map>
#include <list>
#include <utility>
#include <mutex>

// Simple thread-safe LRU cache
template <class Key, class Value>
class LRUCache {
public:
    explicit LRUCache(size_t cap) : capacity(cap) {}

    bool get(const Key& k, Value& v) {
        std::lock_guard<std::mutex> lock(mtx);
        auto it = map.find(k);
        if(it==map.end()) return false;
        order.splice(order.begin(), order, it->second.second);
        v = it->second.first;
        return true;
    }

    void put(const Key& k, const Value& v) {
        std::lock_guard<std::mutex> lock(mtx);
        auto it = map.find(k);
        if(it!=map.end()) {
            it->second.first = v;
            order.splice(order.begin(), order, it->second.second);
            return;
        }
        if(map.size()>=capacity) {
            auto last = order.back();
            map.erase(last);
            order.pop_back();
        }
        order.push_front(k);
        map.emplace(k, std::make_pair(v, order.begin()));
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mtx);
        return map.size();
    }

private:
    size_t capacity;
    mutable std::mutex mtx;
    std::list<Key> order;
    std::unordered_map<Key, std::pair<Value, typename std::list<Key>::iterator>> map;
};

#endif
