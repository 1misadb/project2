#pragma once
#include "geometry.h"
#include <vector>

bool overlapBatchCUDA(const std::vector<Paths64>& A,
                      const std::vector<Paths64>& B,
                      std::vector<bool>& out);
bool cudaAvailable();
