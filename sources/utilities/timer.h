/*
    This file is part of the Tetrahedral Trees library.

    Author(s): Song Ho Ahn (song.ahn@gmail.com) [original author]
               Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Tetrahedral Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Tetrahedral Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Tetrahedral Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <iostream>

#ifndef TIMER_H_DEF
#define TIMER_H_DEF

#ifdef WIN32   // Windows system specific
#include <windows.h>
#else          // Unix based system specific
#include <sys/time.h>
#endif

/**
 * @brief A class representing an high resolution timer.
 * This timer is able to measure the elapsed time with 1 micro-second accuracy in both Windows, Linux and Unix system
 */
class Timer
{
public:
    /**
     * @brief A constructor
     */
    Timer();
    /**
     * @brief A destructor
     */
    ~Timer() {}

    /**
     * @brief A public procedure that starts a timer
     */
    void start();
    /**
     * @brief A public procedure that stops a timer
     */
    void stop();
    /**
     * @brief A public procedure that returns the elapsed time in second
     *
     * @return double
     */
    inline double get_elapsed_time() { return this->get_elapsed_time_in_sec(); }
    /**
     * @brief A public procedure that returns the elapsed time in second
     *
     * @return double
     */
    inline double get_elapsed_time_in_sec() { return this->get_elapsed_time_in_microsec() * 0.000001; }
    /**
     * @brief A public procedure that returns the elapsed time in milli-second
     *
     * @return double
     */
    inline double get_elapsed_time_in_millisec() { return this->get_elapsed_time_in_microsec() * 0.001; }
    /**
     * @brief A public procedure that returns the elapsed time in micro-second
     *
     * @return double
     */
    double get_elapsed_time_in_microsec();
    /**
     * @brief A public procuderes that prints the elapsed time (in seconds) with a user-defined caption
     *
     * @param str, a string referring to the caption
     */
    inline void print_elapsed_time(std::string str) { std::cerr<<str<<get_elapsed_time_in_sec()<<std::endl; }

protected:

private:
    /// A private variable representing starting time in micro-second
    double startTimeInMicroSec;
    /// A private variable representing ending time in micro-second
    double endTimeInMicroSec;
    /// A private variable representing a stop flag
    int    stopped;
#ifdef WIN32
    LARGE_INTEGER frequency;                    // ticks per second
    LARGE_INTEGER startCount;                   //
    LARGE_INTEGER endCount;                     //
#else
    /// A private variable representing the starting tick of the timer
    timeval startCount;
    /// A private variable representing the ending tick of the timer
    timeval endCount;
#endif
};

#endif // TIMER_H_DEF
