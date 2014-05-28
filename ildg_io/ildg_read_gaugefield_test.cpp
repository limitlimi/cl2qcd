/** @file
 * Testcases for the ildg I/O class
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ildg_read_gaugefield
#include <boost/test/unit_test.hpp>

#include "ildg_read_gaugefield.h"

int expectedPrecision = 64;
std::string nameOfExistingGaugefieldFile = std::string(SOURCEDIR) + "/ildg_io/conf.00200";

void checkDefaults(sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.ly_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lz_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lt_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.prec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.flavours_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.beta_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver_source, 0);
}

void checkMetadataOfSpecificGaugefieldFile(sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.ly_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lz_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lt_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.prec_source, 64);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries_source, 18432);
  BOOST_REQUIRE_EQUAL(toCheck.flavours_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr_source, 200);
  /**
   * Some values (time!) are not properly written to file, or their values are not correctly gathered.
   * TODO: repair!
   */
  BOOST_REQUIRE_EQUAL(toCheck.time_source, -619635472);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver_source, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver_source, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue_source, 0.571077);
  BOOST_REQUIRE_EQUAL(toCheck.beta_source, 5.69);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_source, 0.125);
  BOOST_REQUIRE_EQUAL(toCheck.mu_source, 0.006);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar_source, 0);
}

BOOST_AUTO_TEST_CASE(defaults)
{
  sourcefileparameters srcFileParams;
  checkDefaults(srcFileParams);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithFileException)
{
  sourcefileparameters srcFileParams;
  std::string nameOfNonexistingGaugefieldFile = "thisfileshouldnotbethere";
  char * bufferToStoreGaugefield;
  BOOST_CHECK_THROW(srcFileParams.readsourcefile(nameOfNonexistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield), File_Exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithWrongPrecision)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  int wrongPrecision = 27;
  BOOST_CHECK_THROW(srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), wrongPrecision, &bufferToStoreGaugefield), std::exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldSuccess)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  BOOST_REQUIRE_NO_THROW(srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield));
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckMetadata)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  checkMetadataOfSpecificGaugefieldFile(srcFileParams);
}

//todo: this test probably will never work!
BOOST_AUTO_TEST_CASE(readInGaugefieldCheckBufferSize)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
	srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  size_t expectedSizeOfBuffer = srcFileParams.num_entries_source * sizeof(hmc_float);
  size_t actualSizeOfBuffer = sizeof(bufferToStoreGaugefield);
	BOOST_REQUIRE_EQUAL(expectedSizeOfBuffer, actualSizeOfBuffer);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckChecksum)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  uint32_t referenceChecksumA = 171641288;
  uint32_t referenceChecksumB = 3618036129;
  srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.checksum, true);
}




