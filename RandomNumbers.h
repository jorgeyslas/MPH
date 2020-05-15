//  Class for random number generation

#ifndef RandomNumbers_h
#define RandomNumbers_h

#include <random>

// Abstrac class - This is the base for any uniform random number generator
class AbsRNG {
public:
  virtual double operator()() = 0;
};

//  MTRNG uses the MT algorithm for number generation
class MTRNG: public AbsRNG {
    std::mt19937 m_MT;
public:
    // Constructor with an int as input
    MTRNG(std::uint_fast32_t seed = 1);
    // Operator ()  generates the random numbers
    double operator()();
};


#endif /* RandomNumbers_h */
