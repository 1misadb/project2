#pragma once
#include <unordered_map>
#include <list>
#include <mutex>

// Simple thread-safe LRU cache
template<typename Key, typename Value>
class LRUCache {
public:
    explicit LRUCache(size_t cap) : capacity_(cap) {}
    bool get(const Key& k, Value& v) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = map_.find(k);
        if(it == map_.end()) return false;
        items_.splice(items_.begin(), items_, it->second);
        v = it->second->second;
        return true;
    }
    void put(const Key& k, const Value& v) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = map_.find(k);
        if(it != map_.end()) {
            it->second->second = v;
            items_.splice(items_.begin(), items_, it->second);
            return;
        }
        items_.emplace_front(k, v);
        map_[k] = items_.begin();
        if(items_.size() > capacity_) {
            auto last = items_.end();
            --last;
            map_.erase(last->first);
            items_.pop_back();
        }
    }
private:
    size_t capacity_;
    std::list<std::pair<Key, Value>> items_;
    std::unordered_map<Key, typename std::list<std::pair<Key, Value>>::iterator> map_;
    std::mutex mutex_;
};
