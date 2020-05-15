//Class for time measurament

#ifndef Timer_h
#define Timer_h

#include <chrono>

class Timer {
private:
    // Type aliases to make accessing nested type easier
    using clock_t = std::chrono::high_resolution_clock;
    using second_t = std::chrono::duration<double, std::ratio<1>>;
    std::chrono::time_point<clock_t> m_beg;
 
public:
    Timer() : m_beg(clock_t::now()) {} //Initiates the timer
    void reset() { // Resets timer
        m_beg = clock_t::now();
    }
    double elapsed() const { // Time elapsed
        return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

#endif /* Timer_h */
