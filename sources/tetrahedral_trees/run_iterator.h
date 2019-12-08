/*
    This file is part of the Tetrahedral Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Kenneth Weiss (weiss27@llnl.gov)

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

#ifndef RUN_ITERATOR_H
#define RUN_ITERATOR_H

#include <vector>
#include <list>
#include <utility>
#include <iostream>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits.hpp>
#include <boost/integer.hpp>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
/////  Iterator to go through *runs* of integers
/////	positive integers indicate a run of length 1
/////	negative integers indicate a larger run, whose length is given by the following number
/////  We assume a 1-based array (since 0 does not have a negative)
/////  We also assume that the length of a run includes the final index
/////		e.g.  the following run (-5,3) corresponds to the following sequence (5,6,7,8)
/////		-- which includes 3 numbers (6,7,8) after the original number (5)
/////  R. Fellegara - K. Weiss --  March 2013
///////////////////////////////////////////////////////////////////////////////////////

// Predeclare types
template<typename Incrementable> class run_iterator;

typedef std::vector<int> int_vect;
typedef int_vect::iterator int_vect_iter;
typedef int_vect::const_iterator int_vect_const_iter;
typedef run_iterator<int> RunIterator;
typedef std::pair<RunIterator, RunIterator> RunIteratorPair;

/**
 * @brief A class that implements the expansion of compressed arrays
 *
 */
template<typename Incrementable>
class run_iterator : public boost::iterator_facade<
        run_iterator <Incrementable>
        , Incrementable
        , boost::forward_traversal_tag
        , Incrementable
        >
{
public:
    /**
     * @brief
     *
     */
    typedef typename boost::iterator_facade< run_iterator<Incrementable>, Incrementable, boost::forward_traversal_tag> super_t;
    /**
     * @brief
     *
     */
    typedef std::vector<Incrementable>				RunContainer;		// note -- hardcoded for vectors -- if other types are needed, we can change this
    /**
     * @brief
     *
     */
    typedef typename RunContainer::iterator			RunContainerIt;
    /**
     * @brief
     *
     */
    typedef typename RunContainer::const_iterator	RunContainerCIt;

private:
    /**
     * @brief
     *
     */
    enum RunMode { NO_RUN, IN_RUN };

public:

    /**
     * @brief A constructor method
     *
     * @param _runIter
     */
    run_iterator(RunContainerIt const& _runIter = RunContainerIt())
        : runIter(_runIter)
        , runIterEnd(_runIter)
        , currentElement(Incrementable())
        , runElementsRemaining(0)
        , runMode(NO_RUN)
    {}
    /**
     * @brief A constructor method
     *
     * @param _runIter
     */
    run_iterator(RunContainerCIt const& _runIter)
        : runIter(_runIter)
        , runIterEnd(_runIter)
        , currentElement(Incrementable())
        , runElementsRemaining(0)
        , runMode(NO_RUN)
    {}

    /**
     * @brief  A constructor method
     *
     * @param _runIter refers to the beginning of the vector
     * @param _runIterEnd refers to the end of the vector
     */
    run_iterator(RunContainerIt const& _runIter, RunContainerIt const& _runIterEnd)
        : runIter(_runIter)
        , runIterEnd(_runIterEnd)
        , currentElement(0)
        , runElementsRemaining(0)
        , runMode(NO_RUN)
    {
        updateElement();
    }

    /**
     * @brief  A constructor method
     *
     * @param _runIter refers to the beginning of the vector
     * @param _runIterEnd refers to the end of the vector
     */
    run_iterator(RunContainerCIt const& _runIter, RunContainerCIt const& _runIterEnd)
        : runIter(_runIter)
        , runIterEnd(_runIterEnd)
        , currentElement(0)
        , runElementsRemaining(0)
        , runMode(NO_RUN)
    {
        updateElement();
    }

    /**
     * @brief A copy-constructor method
     *
     * @param other
     */
    template<typename OtherIncremental> run_iterator(run_iterator<OtherIncremental> const& other)
        : runIter(other.runIter)
        , runIterEnd(other.runIterEnd)
        , currentElement(other.currentElement)
        , runElementsRemaining(other.runElementsRemaining)
        , runMode(other.runMode)
    {}

    /**
     * @brief A public static method that return a pair containing the begin and end iterator of the vector
     *
     * @param intvec the input vector
     * @return RunIteratorPair the pair containing the begin and end iterators
     */
    static inline RunIteratorPair make_run_iterator_pair(std::vector<Incrementable> const& intvec)
    {
        return std::make_pair( RunIterator(intvec.begin(), intvec.end()), RunIterator(intvec.end()) );
    }

    /**
     * @brief A public method that counts the real number of elements contained in a compressed vector
     *
     * The procedure simply expand each element at once
     *
     * @param intvec
     * @return int
     */
    int elementCount(std::vector<Incrementable> const& intvec);
    /**
     * @brief A public method that counts the real number of elements contained in a compressed vector
     *
     * The procedure visit each entry once without expanding the runs, as it increment the counter inline when encountering a run
     *
     * @param intvec the input vector
     * @return int
     */
    int elementCountFast(std::vector<Incrementable> const& intvec);

private:
    friend class boost::iterator_core_access;

    /**
     * @brief A private method that increment the iterator
     *
     */
    void increment()
    {
        // go to next element
        switch( runMode)
        {
        case NO_RUN:
            ++runIter;
            break;
        case IN_RUN:
            if( runElementsRemaining == 0)
            {
                runMode = NO_RUN;
                ++runIter;
            }
            else
            {
                --runElementsRemaining;
            }
            break;
        }

        updateElement();
    }

    /**
     * @brief A private method that updates the currentElement variable
     *
     */
    void updateElement()
    {
        //  This is meant to prevent dereferences iterators pointing to the end of the vector
        if( runIter == runIterEnd)
            return;

        // Either get the current value or the run value
        switch( runMode)
        {
        case NO_RUN:
            currentElement = *runIter;
            if(currentElement < 0)		// negative indexes indicate the start of a run, the next element gives the run length
            {
                runMode = IN_RUN;
                currentElement = -currentElement;
                runElementsRemaining = *(++runIter);
            }
            break;
        case IN_RUN:
            currentElement = *runIter  - runElementsRemaining - *(runIter-1) ;
            break;
        }
    }

    /**
     * @brief A private method that returns the currentElement variable
     *
     * @return const Incrementable
     */
    Incrementable const& dereference() const { return currentElement; }

    /**
     * @brief A private method that compares two different run_iterator variables
     *
     * Nota: we are not comparing against the currentElement
     *
     * @param other
     * @return bool
     */
    template<typename OtherIncremental> bool equal(run_iterator<OtherIncremental> const& other) const
    {
        switch(runMode)
        {
        case IN_RUN:
            return (runIter == other.runIter)
                    && (IN_RUN == other.runMode)
                    && (runElementsRemaining == other.runElementsRemaining)
                    ;
        case NO_RUN:
        default:
            return (runIter == other.runIter)
                    && (NO_RUN == other.runMode)
                    ;
        }
    }

private:
    RunContainerCIt	runIter; /**< A private variable representing a constant iterator to the beginning of the vector */
    RunContainerCIt	runIterEnd; /**< A private variable representing a constant iterator to the end of the vector */

    Incrementable currentElement; /**< A private variable representing the current element of the run */
    Incrementable runElementsRemaining; /**< A private variable representing a decreasing counter of the elements remaining in a counter */

    RunMode runMode; /**< A private variable representing the current run mode. The mode can be IN_RUN or NO_RUN */
};

template<typename Incrementable> int run_iterator<Incrementable>::elementCount(std::vector<Incrementable> const& intvec)
{
    int count = 0;
    for(RunIteratorPair itPair = make_run_iterator_pair(intvec); itPair.first != itPair.second; ++itPair.first)
        ++count;
    return count;
}

template<typename Incrementable> int run_iterator<Incrementable>::elementCountFast(std::vector<Incrementable> const& intvec)
{
    int count = 0;
    for(typename run_iterator<Incrementable>::RunContainerCIt it = intvec.begin(); it != intvec.end(); ++it)
        count += (*it > 0) ? 1 : 1 + *(++it) ;
    return count;
}


#endif // RUN_ITERATOR_H
