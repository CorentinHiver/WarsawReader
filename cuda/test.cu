#include "../AnalysisLib/CFD.hpp"
#include <cuda_runtime.h>
#include <memory>
#include "../LibCo/print.hpp"
#include "../LibCo/Timer.hpp"
#include "../LibCo/randomCo.hpp"

/*
So, overall it works. However the process is mainly sped up if the resources are correcly managed.
*/

__global__
void batch_cfd_kernel(const double* traces, double* cfds, int shift, double fraction, int n_samples, int output_samples) {
    int trace_idx = blockIdx.x;
    int thread_idx = threadIdx.x;
    int stride = blockDim.x;

    const double* trace = &traces[trace_idx * n_samples];
    double* cfd = &cfds[trace_idx * output_samples];

    for (int i = thread_idx; i < output_samples; i += stride) {
        int bin = i + 2 * shift;
        cfd[i] = fraction * trace[bin] - trace[bin - shift];
    }
}


class BatchCFD {
public:
  BatchCFD() noexcept = default;

  BatchCFD(const std::vector<std::vector<double>>& traceVec, int shift, double fraction)
      : n_traces(traceVec.size()), shift(shift), fraction(fraction)
  {
    n_samples = traceVec[0].size();  // Assume all traces are same length
    output_samples = n_samples - 2 * shift;

    size_t trace_bytes = sizeof(double) * n_traces * n_samples;
    size_t cfd_bytes = sizeof(double) * n_traces * output_samples;

    std::vector<double> trace_flat(n_traces * n_samples);
    for (size_t i = 0; i < n_traces; ++i)
        std::copy(traceVec[i].begin(), traceVec[i].end(), trace_flat.begin() + i * n_samples);

    // Allocate and copy to GPU
    cudaMalloc(&d_traces, trace_bytes);
    cudaMalloc(&d_cfds, cfd_bytes);
    cudaMemcpy(d_traces, trace_flat.data(), trace_bytes, cudaMemcpyHostToDevice);

    // Launch kernel
    int threadsPerBlock = 256;
    batch_cfd_kernel<<<n_traces, threadsPerBlock>>>(d_traces, d_cfds, shift, fraction, n_samples, output_samples);
    cudaDeviceSynchronize();

    // Copy results back
    cfds.resize(n_traces * output_samples);
    cudaMemcpy(cfds.data(), d_cfds, cfd_bytes, cudaMemcpyDeviceToHost);

    // Free memory
    cudaFree(d_traces);
    cudaFree(d_cfds);
  }

  std::vector<std::vector<double>> getResults() const 
  {
    std::vector<std::vector<double>> result(n_traces, std::vector<double>(output_samples));
    for (size_t i = 0; i < n_traces; ++i) 
      std::copy(cfds.begin() + i * output_samples, cfds.begin() + (i + 1) * output_samples, result[i].begin());
    return result;
  }

  std::vector<CFD*> cfds;

private:

  int n_traces, n_samples, shift, output_samples;
  double fraction;

  double* d_traces = nullptr;
  double* d_cfds = nullptr;
  std::vector<double> cfds;
};


int main() {
  int n_traces = 10000;
  int n_samples = 512;
  int shift = 2;
  double fraction = 0.5;

  std::vector<CFD*> traces(n_traces);
  BatchCFD cfds;
  for (int i = 0; i < n_traces; ++i)
  {
    traces[i] = new CFD();
    traces[i]->cfd.resize(n_samples);
    for (int j = 0; j < n_samples; ++j) traces[i]->cfd[j] = std::sin(j / 10.0 + i * 0.001);
    cfds.cfds.push_back(traces[i]);
  }
  
  {
    BatchCFD cfds(traces, shift, fraction);
    Timer timer;
    auto results = cfds.getResults();
    std::cout << "Time: " << timer() << " seconds\n";
  }

  for (auto const & trace : traces) cfd_o.push_back(trace);

  Timer timer;
  for (auto & cfd : cfd_o) cfd.calculate(shift, fraction);
  std::cout << "Time: " << timer() << " seconds\n";
}

//nvcc -Xptxas -O3,-v -std=c++17 -o cfd_gpu test.cu 