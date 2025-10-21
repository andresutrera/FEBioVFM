// diag/printers/param_table.hpp
#pragma once
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include "diag/felog_bridge.hpp" // provides feLog via GetFEModel()
#include "state/vfm_state.hpp"

namespace diag::printers
{

    namespace detail
    {
        inline void line(size_t wname, int wval, char c = '-')
        {
            const std::string n(wname, c), v(wval, c);
            feLog("+%.*s+%.*s+%.*s+%.*s+\n",
                  (int)wname, n.c_str(), wval, v.c_str(), wval, v.c_str(), wval, v.c_str());
        }
        inline void center(const char *title, int inner)
        {
            if (!title)
                title = "";
            const int len = (int)std::strlen(title);
            const int padL = std::max(0, (inner - len) / 2);
            const int padR = std::max(0, inner - len - padL);
            const std::string L(padL, '='), R(padR, '=');
            feLog(" %.*s%s%.*s \n", padL, L.c_str(), title, padR, R.c_str());
        }
        template <class Range, class NameFn>
        inline size_t name_width(const Range &xs, NameFn name_of)
        {
            size_t w = 4; // "Name"
            for (const auto &x : xs)
                w = std::max(w, name_of(x).size());
            return w;
        }
    } // namespace detail

    // Generic table
    template <class Range, class NameFn, class ValFn, class MinFn, class MaxFn>
    inline void ParameterTableGeneric(const Range &items,
                                      const char *title,
                                      int precision,
                                      NameFn name_of,
                                      ValFn val_of,
                                      MinFn min_of,
                                      MaxFn max_of,
                                      int value_col_width = 20)
    {
        const size_t wname = detail::name_width(items, name_of);
        const int W = value_col_width;
        const int inner = (int)wname + 3 * W + 3;

        detail::center(title, inner);
        detail::line(wname, W, '-');
        feLog("|%-*s|%*s|%*s|%*s|\n",
              (int)wname, "Name", W, "Value", W, "Min", W, "Max");
        detail::line(wname, W, '-');

        for (const auto &it : items)
        {
            const std::string nm = name_of(it);
            const double v = val_of(it);
            const double lo = min_of(it);
            const double hi = max_of(it);
            feLog("|%-*s|%*.*g|%*.*g|%*.*g|\n",
                  (int)wname, nm.c_str(),
                  W, precision, v,
                  W, precision, lo,
                  W, precision, hi);
        }
        detail::line(wname, W, '-');
    }

    // Concrete overload for VFMState::params
    inline void ParameterTable(const std::vector<::VFMParam> &params,
                               const char *title,
                               int precision,
                               int value_col_width = 20)
    {
        auto name_of = [](const ::VFMParam &p) -> std::string
        { return p.spec.name; };
        auto val_of = [](const ::VFMParam &p) -> double
        { return p.value; };
        auto min_of = [](const ::VFMParam &p) -> double
        { return p.spec.lo; };
        auto max_of = [](const ::VFMParam &p) -> double
        { return p.spec.hi; };

        ParameterTableGeneric(params, title, precision, name_of, val_of, min_of, max_of, value_col_width);
    }

} // namespace diag::printers