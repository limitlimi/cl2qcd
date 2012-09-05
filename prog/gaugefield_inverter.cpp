#include "gaugefield_inverter.h"

#include "meta/util.hpp"


Opencl_Module_Fermions* Gaugefield_inverter::get_task_solver()
{
	return (Opencl_Module_Fermions*)opencl_modules[task_solver];
}

Opencl_Module_Correlator* Gaugefield_inverter::get_task_correlator()
{
	return (Opencl_Module_Correlator*)opencl_modules[task_correlator];
}

void Gaugefield_inverter::init_tasks()
{
	//allocate host-memory for solution- and source-buffer
	int num_sources;
	if(get_parameters().get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters().get_num_sources();

	size_t bufsize = num_sources * meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize / 1024. / 1024. / 1024. << " GByte";
	solution_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];
	source_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];

	task_solver = 0;
	task_correlator = 1;

	opencl_modules = new Opencl_Module* [get_num_tasks()];


	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_solver] = new Opencl_Module_Fermions(get_parameters());
	get_task_solver()->init(queue[task_solver], get_max_compute_units(task_solver), get_double_ext(task_solver), task_solver);

	opencl_modules[task_correlator] = new Opencl_Module_Correlator(get_parameters());
	get_task_correlator()->init(queue[task_correlator], get_max_compute_units(task_correlator), get_double_ext(task_correlator), task_correlator);


	int spinorfield_size = sizeof(spinor) * meta::get_spinorfieldsize(get_parameters());

	clmem_corr = get_task_correlator()->create_rw_buffer(spinorfield_size * num_sources);
	clmem_source = get_task_correlator()->create_rw_buffer(spinorfield_size);
}

void Gaugefield_inverter::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
}

void Gaugefield_inverter::finalize_opencl()
{
	/// @todo this must be generalized if more than one device is used for one task
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		opencl_modules[ntask]->finalize();
	}
	Gaugefield_hybrid::finalize_opencl();

	cl_int clerr = clReleaseMemObject(clmem_corr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	clerr = clReleaseMemObject(clmem_source);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);


	logger.debug() << "free solution buffer";
	delete [] solution_buffer;
	logger.debug() << "free source buffer";
	delete [] source_buffer;
}

void Gaugefield_inverter::sync_solution_buffer()
{
	size_t sfsize = 12 * meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	get_task_correlator()->copy_buffer_to_device(solution_buffer, get_clmem_corr(), sfsize);
}

void Gaugefield_inverter::perform_inversion(usetimer* solver_timer)
{
	int use_eo = get_parameters().get_use_eo();

	//decide on type of sources
	int num_sources;
	if(get_parameters().get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters().get_num_sources();

	Opencl_Module_Fermions * solver = get_task_solver();

	//allocate host-memory for tmp-buffer
	size_t sfsize = meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	spinor* sftmp = new spinor [meta::get_spinorfieldsize(get_parameters())];

	int spinorfield_size = sizeof(spinor) * meta::get_spinorfieldsize(get_parameters());
	cl_mem clmem_res = solver->create_rw_buffer(spinorfield_size);

	//apply stout smearing if wanted
	if(get_parameters().get_use_smearing() == true) {
		solver->smear_gaugefield(solver->get_gaugefield(), NULL);
	}

	//for CG, one needs a hermitian matrix...
	if(get_parameters().get_use_cg() == true) {
		logger.fatal() << "CG usage requires a hermitian matrix. This is not implemented yet...";
		//the call shoul be like this
		//::QplusQminus_eo f_eo(solver);
	}
	//::Aee f_eo(solver);
	::QplusQminus_eo f_eo(solver);
	::M f_neo(solver);
	Matrix_Function & f = (use_eo) ? static_cast<Matrix_Function &>(f_eo) : static_cast<Matrix_Function &>(f_neo);

	for(int k = 0; k < num_sources; k++) {
		//copy source from to device
		//NOTE: this is a blocking call!
		logger.debug() << "copy pointsource between devices";
		solver->copy_buffer_to_device(&source_buffer[k * meta::get_vol4d(get_parameters())], get_clmem_source(), sfsize);

		logger.debug() << "calling solver..";
		solver->solver(f, clmem_res, get_clmem_source(), solver->get_gaugefield(), solver_timer);

		//add solution to solution-buffer
		//NOTE: this is a blocking call!
		logger.debug() << "add solution...";
		solver->get_buffer_from_device(clmem_res, &solution_buffer[k * meta::get_vol4d(get_parameters())], sfsize);
	}

	if(get_parameters().get_use_smearing() == true) {
		solver->unsmear_gaugefield(solver->get_gaugefield());
	}

	delete [] sftmp;
	cl_int clerr = clReleaseMemObject(clmem_res);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
}

void Gaugefield_inverter::flavour_doublet_correlators(std::string corr_fn)
{
	using namespace std;

	//for now, make sure clmem_corr is properly filled; maybe later we can increase performance a bit by playing with this...
	sync_solution_buffer();

	//suppose that the buffer on the device has been filled with the prior calculated solutions of the solver
	logger.debug() << "start calculating correlators...";

	ofstream of;
	of.open(corr_fn.c_str(), ios_base::app);
	if( !of.is_open() ) throw File_Exception(corr_fn);
	of << "# flavour doublet correlators" << endl;
	if(get_parameters().get_corr_dir() == 3) {
		of << "# format: J P z real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), z spatial distance, value (aggregate x y z)" << endl;
	} else {
		of << "# format: J P t real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), t timelike distance, value (aggregate x y z)" << endl;
	}


	int num_corr_entries =  0;
	switch (get_parameters().get_corr_dir()) {
		case 0 :
			num_corr_entries = get_parameters().get_ntime();
			break;
		case 3 :
			num_corr_entries = get_parameters().get_nspace();
			break;
		default :
			stringstream errmsg;
			errmsg << "Correlator direction " << get_parameters().get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}

	size_t buffersize = sizeof(hmc_float) * num_corr_entries;

	hmc_float* host_result_ps = new hmc_float [num_corr_entries];
	hmc_float* host_result_sc = new hmc_float [num_corr_entries];
	hmc_float* host_result_vx = new hmc_float [num_corr_entries];
	hmc_float* host_result_vy = new hmc_float [num_corr_entries];
	hmc_float* host_result_vz = new hmc_float [num_corr_entries];
	hmc_float* host_result_ax = new hmc_float [num_corr_entries];
	hmc_float* host_result_ay = new hmc_float [num_corr_entries];
	hmc_float* host_result_az = new hmc_float [num_corr_entries];

	cl_mem result_ps;
	cl_mem result_sc;
	cl_mem result_vx;
	cl_mem result_vy;
	cl_mem result_vz;
	cl_mem result_ax;
	cl_mem result_ay;
	cl_mem result_az;
	//LZ usually the correlator is calculated on CPU, then we don't need to copy all those buffers...
	bool needcopy = true;
	if( get_task_correlator()->get_device_type() == CL_DEVICE_TYPE_CPU ) {
		needcopy = false;
		result_ps = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ps);
		result_sc = get_task_correlator()->create_uhp_buffer(buffersize, host_result_sc);
		result_vx = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vx);
		result_vy = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vy);
		result_vz = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vz);
		result_ax = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ax);
		result_ay = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ay);
		result_az = get_task_correlator()->create_uhp_buffer(buffersize, host_result_az);
	} else {
		result_ps = get_task_correlator()->create_rw_buffer(buffersize);
		result_sc = get_task_correlator()->create_rw_buffer(buffersize);
		result_vx = get_task_correlator()->create_rw_buffer(buffersize);
		result_vy = get_task_correlator()->create_rw_buffer(buffersize);
		result_vz = get_task_correlator()->create_rw_buffer(buffersize);
		result_ax = get_task_correlator()->create_rw_buffer(buffersize);
		result_ay = get_task_correlator()->create_rw_buffer(buffersize);
		result_az = get_task_correlator()->create_rw_buffer(buffersize);
	}

	logger.info() << "calculate correlators..." ;
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ps"), get_clmem_corr(), result_ps);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("sc"), get_clmem_corr(), result_sc);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vx"), get_clmem_corr(), result_vx);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vy"), get_clmem_corr(), result_vy);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vz"), get_clmem_corr(), result_vz);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ax"), get_clmem_corr(), result_ax);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ay"), get_clmem_corr(), result_ay);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("az"), get_clmem_corr(), result_az);

	//the correlator_device calcultions are all non-blocking, hence we need a barrier here:
	if( needcopy == false )
		clFinish(queue[task_correlator]);

	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	if( needcopy == true )
		get_task_correlator()->get_buffer_from_device(result_ps, host_result_ps, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << host_result_ps[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << host_result_ps[j] << endl;
	}


	//the scalar (J=0, P=0)
	if( needcopy == true )
		get_task_correlator()->get_buffer_from_device(result_sc, host_result_sc, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << host_result_sc[j] << endl;
	}


	//the vector (J=1, P=1)
	if( needcopy == true ) {
		get_task_correlator()->get_buffer_from_device(result_vx, host_result_vx, buffersize);
		get_task_correlator()->get_buffer_from_device(result_vy, host_result_vy, buffersize);
		get_task_correlator()->get_buffer_from_device(result_vz, host_result_vz, buffersize);
	}
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (host_result_vx[j] + host_result_vy[j] + host_result_vz[j]) / 3. << "\t" << host_result_vx[j] << "\t" << host_result_vy[j] << "\t" << host_result_vz[j] << endl;
	}


	//the axial vector (J=1, P=0)
	if( needcopy == true ) {
		get_task_correlator()->get_buffer_from_device(result_ax, host_result_ax, buffersize);
		get_task_correlator()->get_buffer_from_device(result_ay, host_result_ay, buffersize);
		get_task_correlator()->get_buffer_from_device(result_az, host_result_az, buffersize);
	}
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (host_result_ax[j] + host_result_ay[j] + host_result_az[j]) / 3. << "\t" << host_result_ax[j] << "\t" << host_result_ay[j] << "\t" << host_result_az[j] << endl;
	}


	of << endl;
	of.close();
	delete [] host_result_ps;
	delete [] host_result_sc;
	delete [] host_result_vx;
	delete [] host_result_vy;
	delete [] host_result_vz;
	delete [] host_result_ax;
	delete [] host_result_ay;
	delete [] host_result_az;

	cl_int clerr = CL_SUCCESS;
	clerr = clReleaseMemObject(result_ps);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_sc);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vx);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_ax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_ay);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_az);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
}

void Gaugefield_inverter::create_sources()
{
	//create sources on the correlator-device and save them on the host
	size_t sfsize = meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	if(get_parameters().get_use_pointsource() == true) {
		logger.debug() << "start creating point-sources...";
		for(int k = 0; k < 12; k++) {
			get_task_correlator()->create_point_source_device(get_clmem_source(), k, meta::get_source_pos_spatial(get_parameters()), get_parameters().get_pointsource_t());
			logger.debug() << "copy pointsource to host";
			get_task_correlator()->get_buffer_from_device(get_clmem_source(), &source_buffer[k * meta::get_vol4d(get_parameters())], sfsize);
		}
	} else {
		logger.debug() << "start creating stochastic-sources...";
		int num_sources = get_parameters().get_num_sources();
		for(int k = 0; k < num_sources; k++) {
			get_task_correlator()->create_stochastic_source_device(get_clmem_source());
			logger.debug() << "copy stochastic-source to host";
			get_task_correlator()->get_buffer_from_device(get_clmem_source(), &source_buffer[k * meta::get_vol4d(get_parameters())], sfsize);
		}
	}
}

cl_mem Gaugefield_inverter::get_clmem_corr()
{
	return clmem_corr;
}


cl_mem Gaugefield_inverter::get_clmem_source()
{
	return clmem_source;
}
