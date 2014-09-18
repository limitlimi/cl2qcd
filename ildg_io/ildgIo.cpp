/** @file
 * ildg IO utilities
 *
 * Copyright 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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
#include "ildgIo.hpp"

#include "ildgIo_gaugefield.hpp"
#include "../meta/util.hpp"
#include "checksum.h"
#include <cassert>
#include "../meta/version.hpp"
#include "../executables/exceptions.h"

static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const meta::Inputparameters& parameters);
static void copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters);
static Checksum calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters);

Matrixsu3 * ildgIo::readGaugefieldFromSourcefile(std::string ildgfile, const meta::Inputparameters * parameters, int & trajectoryNumberAtInit, double & plaq)
{
	Matrixsu3 * gf_host;
	char * gf_ildg;

	IldgIoReader_gaugefield reader(ildgfile.c_str(), parameters->get_precision(), &gf_ildg);
	
	Checksum checksum = calculate_ildg_checksum(gf_ildg, reader.parameters.getSizeInBytes(), *parameters);

	//todo: this should not be that explicit here!	
	gf_host = new Matrixsu3[meta::get_vol4d(*parameters) * 4];
	copy_gaugefield_from_ildg_format(gf_host, gf_ildg, reader.parameters.num_entries, *parameters);
	delete[] gf_ildg;

	trajectoryNumberAtInit = reader.parameters.trajectorynr;
	plaq = reader.parameters.plaquettevalue;

	//todo: move this to destructor or so...
	reader.parameters.checkAgainstInputparameters(parameters);
	reader.parameters.checkAgainstChecksum(checksum, parameters->get_ignore_checksum_errors(), ildgfile);

	return gf_host;
}

static size_t getBufferSize_gaugefield(const meta::Inputparameters * parameters) noexcept
{
	return 2 * NC * NC * NDIM * meta::get_volspace(*parameters) * parameters->get_ntime() * sizeof(hmc_float);
}

void ildgIo::writeGaugefieldToFile(std::string outputfile, Matrixsu3 * host_buf, const meta::Inputparameters * parameters, int number, double plaq)
{
	const size_t gaugefield_buf_size = getBufferSize_gaugefield(parameters);
	char * gaugefield_buf = new char[gaugefield_buf_size];

	copy_gaugefield_to_ildg_format(gaugefield_buf, host_buf, *parameters);

	const Checksum checksum = calculate_ildg_checksum(gaugefield_buf, gaugefield_buf_size, *parameters);

	sourcefileparameters_values srcFileParameters(parameters, number, plaq, checksum, version);
	
	IldgIoWriter_gaugefield writer(gaugefield_buf, gaugefield_buf_size, srcFileParameters, outputfile);

	delete[] gaugefield_buf;
}

static hmc_float make_float_from_big_endian(const char* in)
{
	union {
		char b[sizeof(hmc_float)];
		hmc_float f;
	} val;

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
		val.b[i] = in[sizeof(hmc_float) - 1 - i];
	}
	return val.f;
}

static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const meta::Inputparameters& parameters)
{
	//little check if arrays are big enough
	if ((int) (meta::get_vol4d(parameters) *NDIM * NC * NC * 2) != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t NSPACE = parameters.get_nspace();
	int cter = 0;
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in a complex array
						hmc_complex tmp [NC][NC];
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								tmp[m][n].re = make_float_from_big_endian(&gaugefield_tmp[pos * sizeof(hmc_float)]);
								tmp[m][n].im = make_float_from_big_endian(&gaugefield_tmp[(pos + 1) * sizeof(hmc_float)]);
								cter++;
							}
						}
						//store su3matrix tmp in our format
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);

						//copy hmc_su3matrix to Matrixsu3 format
						Matrixsu3 destElem;
						destElem.e00 = tmp[0][0];
						destElem.e01 = tmp[0][1];
						destElem.e02 = tmp[0][2];
						destElem.e10 = tmp[1][0];
						destElem.e11 = tmp[1][1];
						destElem.e12 = tmp[1][2];
						destElem.e20 = tmp[2][0];
						destElem.e21 = tmp[2][1];
						destElem.e22 = tmp[2][2];

						gaugefield[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)] = destElem;
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}
}

static Checksum calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters)
{
	const size_t elem_size = 4 * sizeof(Matrixsu3);

	const size_t NT = inputparameters.get_ntime();
	const size_t NS = inputparameters.get_nspace();

	if(nbytes != (NT * NS * NS * NS * elem_size)) {
		logger.error() << "Buffer does not contain a gaugefield!";
		throw Invalid_Parameters("Buffer size not match possible gaugefield size", (NT * NS * NS * NS * elem_size), nbytes);
	}

	Checksum checksum;

	size_t offset = 0;
	for(uint32_t t = 0; t < NT; ++t) {
		for(uint32_t z = 0; z < NS; ++z) {
			for(uint32_t y = 0; y < NS; ++y) {
				for(uint32_t x = 0; x < NS; ++x) {
					assert(offset < nbytes);
					uint32_t rank = ((t * NS + z) * NS + y) * NS + x;
					checksum.accumulate(&buf[offset], elem_size, rank);
					offset += elem_size;
				}
			}
		}
	}

	logger.debug() << "Calculated Checksum: " << checksum;
	return checksum;
}

static void make_big_endian_from_float(char* out, const hmc_float in)
{
	union {
		char b[sizeof(hmc_float)];
		hmc_float f;
	} val;

	val.f = in;

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
		out[i] = val.b[sizeof(hmc_float) - 1 - i];
	}
}

static void copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters)
{
	const size_t NSPACE = parameters.get_nspace();
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);
						hmc_complex destElem [NC][NC];

						Matrixsu3 srcElem = source_in[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)];
						destElem[0][0] = srcElem.e00;
						destElem[0][1] = srcElem.e01;
						destElem[0][2] = srcElem.e02;
						destElem[1][0] = srcElem.e10;
						destElem[1][1] = srcElem.e11;
						destElem[1][2] = srcElem.e12;
						destElem[2][0] = srcElem.e20;
						destElem[2][1] = srcElem.e21;
						destElem[2][2] = srcElem.e22;

						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								make_big_endian_from_float(&dest[pos * sizeof(hmc_float)], destElem[m][n].re);
								make_big_endian_from_float(&dest[(pos + 1) * sizeof(hmc_float)], destElem[m][n].im);
							}
						}
					}
				}
			}
		}
	}
}
