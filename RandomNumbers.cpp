//  Class for random number generation - functions implementation

#include "RandomNumbers.h"

//Class MTRNG
//  Constructor with an int as input
MTRNG::MTRNG(std::uint_fast32_t seed) : m_MT{seed} {}
//  () Operator
double MTRNG::operator()() {
    return m_MT() / static_cast<double>(m_MT.max());
}

