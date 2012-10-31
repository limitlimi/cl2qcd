/** @file
 * PRNG PRNG unit declaration
 */

#ifndef _PHYSICS_PRNG_
#define _PHYSICS_PRNG_

#include "../hardware/system.hpp"
#include "../hardware/buffers/prng_buffer.hpp"

/**
 * This package contains the actual "business" logic of the library.
 */
namespace physics {

	namespace lattices {
		class Gaugefield;
	}

	/**
	 * The PRNG PRNG for host and device.
	 *
	 * WARNING: Must only be used as a singleton!
	 */
	class PRNG {
		friend physics::lattices::Gaugefield;

		public:
			/**
			 * Initialize the ranlux instance
			 */
			PRNG(const hardware::System& system);

			~PRNG();

			/*
			 * non-copyable
			 */
			PRNG& operator=(const PRNG&) = delete;
			PRNG(const PRNG&) = delete;
			PRNG() = delete;

			/**
			 * Get random numbers
			 */
			double get_double() noexcept;

		private:
			/**
			 * Reference to the PRNG Buffers used on each device
			 */
			std::vector<const hardware::buffers::PRNGBuffer*> buffers;

			/**
			 * Reference to the system this PRNG is for.
			 */
			const hardware::System& system;
	};

}

#endif /* _PHYSICS_PRNG_ */