#pragma once
#include "../Colib/lib/libCo.hpp"

template <class Trace>
Trace& getBaseline(Trace & trace, size_t nb_samples_baseline)
{
  int baseline = 0;
  for (size_t i = 0; i<nb_samples_baseline; ++i) baseline += trace[i];
  return baseline / nb_samples_baseline;
}

template <class Trace>
Trace& removeBaseline(Trace & trace, size_t nb_samples_baseline)
{
  auto baseline = getBaseline(trace, nb_samples_baseline);
  for (auto & sample : trace) sample -= baseline;
}

template <class Trace>
Trace& removeBaseline(Trace & trace, size_t nb_samples_baseline, double baseline)
{
  for (auto & sample : trace) sample -= baseline;
}

template <class Trace>
double getRiseTimeBins(Trace & trace, double start = 0.1, double stop = 0.9, double baseline = 0)
{
  auto [min, max] = minAndMax(trace);
  auto const start_threshold = min-baseline;
  auto const stop_threshold = max-baseline;
  bool start_bin{}, stop_bin{};
  size_t bin_start{}, bin_stop{};
  for (size_t bin = 0; bin<trace.size(); ++bin)
  {
    if (!start && start_threshold < trace[bin]) bin_start = bin;
    else if (!stop && stop_threshold < trace[bin]) {bin_stop = bin; break;}
  }
  return bin_start-bin_stop;
}

template <class Trace>
auto interpolate(const Trace& trace, double fractional_bin) -> typename Trace::value_type 
{
  if (fractional_bin <= 0) return trace.front();
  if (fractional_bin >= trace.size() - 1) return trace.back();

  size_t i = static_cast<size_t>(fractional_bin);
  double fraction = fractional_bin - i;

  auto const& y0 = trace[i];
  auto const& y1 = trace[i + 1];

  return static_cast<typename Trace::value_type>(y0 + (y1 - y0) * fraction);
}