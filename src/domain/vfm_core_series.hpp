#pragma once
#include <vector>
#include <cstddef>
using TimeIdx = int;

template <typename FrameT>
class TimeSeries {
public:
    std::size_t nTimes() const { return _t.size(); }
    TimeIdx addTime() { _t.emplace_back(); return (TimeIdx)_t.size() - 1; }
    FrameT&       getTime(TimeIdx t){ return _t[(std::size_t)t]; }
    const FrameT& getTime(TimeIdx t) const { return _t[(std::size_t)t]; }
private:
    std::vector<FrameT> _t;
};