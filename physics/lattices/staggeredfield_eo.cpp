/** @file
 * Implementation of the physics::lattices::Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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


#include "staggeredfield_eo.hpp"
#include "../../meta/util.hpp"
#include <cassert>
#include "../../hardware/code/spinors_staggered.hpp"
//#include "../../hardware/code/fermions_staggered.hpp"
#include "../../meta/type_ops.hpp"
#include "../../hardware/buffers/halo_update.hpp"
//For hardware::code::get_eoprec_spinorfieldsize()
#include "../../hardware/code/spinors.hpp"

static std::vector<const hardware::buffers::SU3vec *> allocate_buffers(const hardware::System& system);

physics::lattices::Staggeredfield_eo::Staggeredfield_eo(const hardware::System& system)
	: system(system), staggaredfieldEoParametersInterface(new StaggaredfieldEoParametersImplementation(system.get_inputparameters())),
	  buffers(allocate_buffers(system))
{
}

static std::vector<const hardware::buffers::SU3vec *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::SU3vec;

	auto devices = system.get_devices();
	std::vector<const SU3vec*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new SU3vec(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device));
	}
	return buffers;
}

physics::lattices::Staggeredfield_eo::~Staggeredfield_eo()
{
    for(auto buffer: buffers) {
		delete buffer;
	}
    //TODO: remove the following delete
    delete staggaredfieldEoParametersInterface;
}

const std::vector<const hardware::buffers::SU3vec *> physics::lattices::Staggeredfield_eo::get_buffers() const noexcept
{
	return buffers;
}

hmc_complex physics::lattices::scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right)
{
	const Scalar<hmc_complex> res(left.system);
	scalar_product(&res, left, right);
	return res.get();
}

hmc_float physics::lattices::scalar_product_real_part(const Staggeredfield_eo& left, const Staggeredfield_eo& right)
{
	const Scalar<hmc_float> res(left.system);
	scalar_product_real_part(&res, left, right);
	return res.get();
}

void physics::lattices::scalar_product(const Scalar<hmc_complex>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right)
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
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->set_complex_to_scalar_product_eoprec_device(left_buf, right_buf, res_buf);
	}
	res->sum();
}

void physics::lattices::scalar_product_real_part(const Scalar<hmc_float>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right)
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
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->set_float_to_scalar_product_real_part_eoprec_device(left_buf, right_buf, res_buf);
	}
	res->sum();
}

hmc_float physics::lattices::squarenorm(const Staggeredfield_eo& field)
{
	const Scalar<hmc_float> res(field.system);
	squarenorm(&res, field);
	return res.get();
}

void physics::lattices::squarenorm(const Scalar<hmc_float>* res, const Staggeredfield_eo& field)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->set_float_to_global_squarenorm_eoprec_device(field_buf, res_buf);
	}
	res->sum();
}

void physics::lattices::Staggeredfield_eo::set_zero() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_staggered_code();
		spinor_code->set_zero_spinorfield_eoprec_device(buffer);
	}
}

void physics::lattices::Staggeredfield_eo::set_cold() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_staggered_code();
		spinor_code->set_cold_spinorfield_eoprec_device(buffer);
	}
}

void physics::lattices::Staggeredfield_eo::set_gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(buffers.size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto spin_buf = buffers[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->get_spinor_staggered_code()->set_gaussian_spinorfield_eoprec_device(spin_buf, prng_buf);
	}
	if((system.get_devices()).size()!=1)
		update_halo();
}

void physics::lattices::sax(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::sax(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x)
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
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::sax(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::sax(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x)
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
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::sax(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x)
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
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], index_alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], index_alpha, out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	
	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const hmc_float beta, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_float>& beta, const Staggeredfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	
	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x, const Vector<hmc_float>& beta, const int index_beta, const Staggeredfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	
	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], beta_bufs[i], index_alpha, index_beta, out_buf);
	}
}

void physics::lattices::saxpbypz(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();
	auto z_bufs = z.get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size() || out_bufs.size() != z_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpbypz(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z)
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
		device->get_spinor_staggered_code()->saxpbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::sax_vec_and_squarenorm(const Vector<hmc_float>* res, const Vector<hmc_float>& alpha, const Staggeredfield_eo& x)
{
	auto x_buffers = x.get_buffers();
	auto alpha_buffers = alpha.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = x_buffers.size();

	if(num_buffers != res_buffers.size() || num_buffers != alpha_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto x_buf = x_buffers[i];
		auto alpha_buf = alpha_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = x_buf->get_device();
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->sax_vectorized_and_squarenorm_eoprec_device(x_buf, alpha_buf, alpha.get_vector_size(), res_buf);
	}
	//res->sum();
}


template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::scalar_product>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("scalar_product_staggered_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::scalar_product_real_part>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("scalar_product_real_part_staggered_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("global_squarenorm_staggered_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("sax_cplx_staggered_eoprec"); //the same as sax_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("sax_real_staggered_eoprec"); //the same as sax_real_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpy_cplx_staggered_eoprec"); //the same as saxpy_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpy_real_staggered_eoprec"); //the same as saxpy_real_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpby>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpby_cplx_staggered_eoprec"); //the same as saxpby_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpby>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpby_real_staggered_eoprec"); //the same as saxpby_real_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::saxpbypz>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpbypz_cplx_staggered_eoprec");//the same as saxpbypz_cplx_arg_stagg_eoprec
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Staggeredfield_eo& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Staggeredfield_eo::update_halo() const
{
	//So far, this is simple a throw
	throw Print_Error_Message("Halo update is not implemented, yet.", __FILE__, __LINE__);
}

//This code is useful for using the template pseudo_randomize defined in lattices/util.hpp
void physics::lattices::Staggeredfield_eo::import(const su3vec * const host) const
{
	logger.trace() << "importing staggeredfield_eo";
	if(buffers.size() == 1) {
		buffers[0]->load(host);
	} else {
		throw Print_Error_Message("Import not implemented for multi device staggeredfield_eo!", __FILE__, __LINE__);
	}
	logger.trace() << "import complete";
}

unsigned physics::lattices::Staggeredfield_eo::get_elements() const noexcept
{
	return meta::get_vol4d(staggaredfieldEoParametersInterface->getNt(), staggaredfieldEoParametersInterface->getNs());
}
