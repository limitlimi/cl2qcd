/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <stdint.h>

/**
 * A utility class to measure times
 * using the monotonic timer at a
 * micro-second accuracy (if the system permits)
 */
class Timer
{
	public:
		/**
		 * Default and only constructor.
		 * Imidiateley starts the measurement.
		 */
		Timer();

		/**
		 * Resets the timer to start at 0 again.
		 */
		void reset();

		/**
		 * Recieves the time passed since timer start/reset.
		 * The value returned is in microseconds.
		 */
		uint64_t getTime();

		/**
		 * Recieves the time passed sind timer start/reset
		 * and resets the timer.
		 * The value returned is in microseconds.
		 */
		uint64_t getTimeAndReset();
		
	private:
		/**
		 * Time of the timer start/reset.
		 * This is the offset that needs to be substracted from
		 * later measurements to get the time difference.
		 */
		uint64_t start;

		/**
		 * Calculates the difference in microseconds between the two events
		 */
		uint64_t getDifference( uint64_t start, uint64_t end ) const;
	

		/**
		 * Retrieves the current Timestamp in an OS-specific manner
		 */
		uint64_t getTimestamp() const;
};

#endif // TIMER_H
