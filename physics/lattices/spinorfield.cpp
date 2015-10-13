/** @file
 * Implementation of the physics::lattices::Spinorfield class
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "spinorfield.hpp"
#include <cassert>
#include <stdexcept>
#include "../../hardware/code/spinors.hpp"
#include "../../hardware/code/fermions.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "../../host_functionality/host_geometry.h"

static std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(const hardware::System& system, const bool place_on_host);

physics::lattices::Spinorfield::Spinorfield(const hardware::System& system, const bool place_on_host)
	: system(system), spinorfieldParametersInterface(new SpinorfieldParametersImplementation(system.get_inputparameters())),
	  buffers(allocate_buffers(system, place_on_host)), place_on_host(place_on_host)
{
}

static  std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(const hardware::System& system, const bool place_on_host)
{
	using hardware::buffers::Plain;

	std::vector<const Plain<spinor>*> buffers;
	for(auto device: system.get_devices()) {
		buffers.push_back(new Plain<spinor>(hardware::code::get_spinorfieldsize(device->get_mem_lattice_size()), device, place_on_host));
	}
	return buffers;
}

physics::lattices::Spinorfield::~Spinorfield()
{
    //TODO: Remove this delete
    delete spinorfieldParametersInterface;
	clear_buffers();
}


void physics::lattices::Spinorfield::clear_buffers()
{
for(auto buffer: buffers) {
		delete buffer;
	}
	buffers.clear();
}

void physics::lattices::Spinorfield::fill_buffers()
{
	if(buffers.size() != 0) {
		return;
	}

	buffers = allocate_buffers(system, place_on_host);
}

std::vector<physics::lattices::Spinorfield *> physics::lattices::create_spinorfields(const hardware::System& system, const size_t n, const bool place_on_host)
{
	std::vector<Spinorfield *> fields;
	fields.reserve(n);

	for(size_t i = 0; i < n; ++i) {
		fields.push_back(new Spinorfield(system, place_on_host));
	}

	return fields;
}

void physics::lattices::release_spinorfields(const std::vector<physics::lattices::Spinorfield *> fields)
{
for(auto field: fields) {
		delete field;
	}
}

const std::vector<const hardware::buffers::Plain<spinor> *> physics::lattices::Spinorfield::get_buffers() const noexcept
{
	return buffers;
}

void physics::lattices::Spinorfield::gamma5() const
{
for(auto buffer: buffers) {
		auto fermion_code = buffer->get_device()->getFermionCode();
		fermion_code->gamma5_device(buffer);
	}
}

hmc_complex physics::lattices::scalar_product(const Spinorfield& left, const Spinorfield& right)
{
	const Scalar<hmc_complex> res(left.system);
	scalar_product(&res, left, right);
	return res.get();
}

void physics::lattices::scalar_product(const Scalar<hmc_complex>* res, const Spinorfield& left, const Spinorfield& right)
{
	auto res_buffers = res->get_buffers();
	auto left_buffers = left.get_buffers();
	auto right_buffers = right.get_buffers();
	size_t num_buffers = res_buffers.size();

	if(num_buffers != left_buffers.size() || num_buffers != right_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto res_buf = res_buffers[i];
		auto left_buf = left_buffers[i];
		auto right_buf = right_buffers[i];
		auto device = res_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_complex_to_scalar_product_device(left_buf, right_buf, res_buf);
	}

	res->sum();
}

hmc_float physics::lattices::squarenorm(const Spinorfield& field)
{
	const Scalar<hmc_float> res(field.system);
	squarenorm(&res, field);
	return res.get();
}

void physics::lattices::squarenorm(const Scalar<hmc_float>* res, const Spinorfield& field)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	// TODO implemente for more than one device
	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not sue the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_float_to_global_squarenorm_device(field_buf, res_buf);
	}

	res->sum();
}

void physics::lattices::Spinorfield::zero() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->getSpinorCode();
		spinor_code->set_zero_spinorfield_device(buffer);
	}
}

void physics::lattices::Spinorfield::cold() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->getSpinorCode();
		spinor_code->set_spinorfield_cold_device(buffer);
	}
}

void physics::lattices::Spinorfield::gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(buffers.size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto spin_buf = buffers[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->getSpinorCode()->generate_gaussian_spinorfield_device(spin_buf, prng_buf);
	}

	update_halo();
}

void physics::lattices::saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxpy_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Spinorfield& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxpy_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x)
{
	const Scalar<hmc_complex> alpha_buf(out->system);
	alpha_buf.store(alpha);
	sax(out, alpha_buf, x);
}

void physics::lattices::sax(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->sax_device(x_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z)
{
	const Scalar<hmc_complex> alpha_buf(out->system);
	const Scalar<hmc_complex> beta_buf(out->system);
	alpha_buf.store(alpha);
	beta_buf.store(beta);
	saxsbypz(out, alpha_buf, x, beta_buf, y, z);
}

void physics::lattices::saxsbypz(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Scalar<hmc_complex>& beta, const Spinorfield& y, const Spinorfield& z)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	auto z_bufs = z.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size() || out_bufs.size() != z_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxsbypz_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Spinorfield& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Spinorfield::update_halo() const
{
	hardware::buffers::update_halo<spinor>(buffers, system);
}

void physics::lattices::Spinorfield::import(const spinor * const host) const
{
	logger.trace() << "importing spinorfield";
	if(buffers.size() == 1) {
		buffers[0]->load(host);
	} else {
		auto const _device = buffers.at(0)->get_device();
		auto const local_size = _device->get_local_lattice_size();
		size_4 const halo_size(local_size.x, local_size.y, local_size.z, _device->get_halo_size());
		auto const grid_size = _device->getGridSize();
		if(grid_size.x != 1 || grid_size.y != 1 || grid_size.z != 1) {
			throw Print_Error_Message("Not implemented!", __FILE__, __LINE__);
		}
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();

			size_4 offset(0, 0, 0, device->getGridPos().t * local_size.t);
			logger.debug() << offset;
			const size_t local_volume = get_vol4d(local_size);
			buffer->load(&host[get_global_pos(offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs())], local_volume);

			const size_t halo_volume = get_vol4d(halo_size);
			size_4 halo_offset(0, 0, 0, (offset.t + local_size.t) % spinorfieldParametersInterface->getNt());
			logger.debug() << halo_offset;
			logger.trace() << get_global_pos(halo_offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs());
			logger.trace() << halo_volume;
			logger.trace() << get_elements();
			assert(get_global_pos(halo_offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs()) + halo_volume <= get_elements());
			buffer->load(&host[get_global_pos(halo_offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs())], halo_volume, local_volume);

			halo_offset = size_4(0, 0, 0, (offset.t + spinorfieldParametersInterface->getNt() - halo_size.t) % spinorfieldParametersInterface->getNt());
			logger.debug() << halo_offset;
			assert(get_global_pos(halo_offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs()) + halo_volume <= get_elements());
			buffer->load(&host[get_global_pos(halo_offset, spinorfieldParametersInterface->getNt(), spinorfieldParametersInterface->getNs())], halo_volume, local_volume + halo_volume);

		}
	}
	logger.trace() << "import complete";
}

unsigned physics::lattices::Spinorfield::get_elements() const noexcept
{
    return spinorfieldParametersInterface->getNumberOfElements();
}

void physics::lattices::fill_window(const physics::lattices::Spinorfield* out, const physics::lattices::Spinorfield& src, const size_t idx)
{
	auto src_buf = src.get_buffers()[idx];
	spinor * const tmp = new spinor[src_buf->get_elements()];

	src_buf->dump(tmp);

	for(auto out_buf: out->get_buffers()) {
		out_buf->load(tmp);
	}

	delete[] tmp;
}
