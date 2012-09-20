/** @file
 * Defines utility classes for kernel compilation for OpenCL devices.
 */

#ifndef _OPENCL_COMPILER_H_
#define _OPENCL_COMPILER_H_

#include <vector>
#include <string>

#include "exceptions.h"

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

/**
 * A collection of OpenCL source files for repeated passing to the OpenCL compiler.
 */
class ClSourcePackage {

public:
	/**
	 * Create an empty package with the given options
	 */
	ClSourcePackage(std::string options = std::string()) : options(options) {};

	/**
	 * Copy constructor
	 */
	ClSourcePackage(const ClSourcePackage& src) : files(src.files), options(src.options) {};

	/**
	 * Create a package based from an array of filenames and potentially
	 * another package
	 *
	 * @param names An array of filenames.
	 * @param num How many filenames.
	 * @param base Optionally a package to base this package on.
	 */
	ClSourcePackage(const std::vector<const char *>& files, std::string options) : files(files), options(options) {};

	/**
	 * Add another source file to the list of sources.
	 */
	ClSourcePackage operator <<(const char *file);

	/**
	 * Add a predefined group of packages.
	 */
	ClSourcePackage operator <<(const ClSourcePackage& package);

	/**
	 * Allow assignment of packages
	 */
	ClSourcePackage operator =(const ClSourcePackage& package);

	/**
	 * Get the list of files within the package.
	 */
	const std::vector<const char *> getFiles() const;

	/**
	 * Get the list of files within the package.
	 */
	const std::string getOptions() const;

private:
	/**
	 * Collection of all the filenames part of the package
	 */
	std::vector<const char *> files;
	/**
	 * Collection of all build options set for this package
	 */
	std::string options;
};

/**
 * Parameter collection object for OpenCL compile.
 */
class TmpClKernel {

public:

	/**
	 * All purpose constructor.
	 */
	TmpClKernel(const std::string kernel_name, const std::string build_options,
	            const cl_context context, cl_device_id device,
	            const std::vector<const char *> files = std::vector<const char *>())
		: kernel_name(kernel_name), build_options(build_options), context(context),
		  device(device), files(files) { };

	/**
	 * Conversion operator to a real OpenCL kernel object, will trigger build.
	 */
	operator cl_kernel() const;

	/**
	 * Add another source file to the list of sources.
	 */
	TmpClKernel operator <<(const char *file) const;

	/**
	 * Add a predefined group of packages.
	 */
	TmpClKernel operator <<(const ClSourcePackage& package) const;

private:
	/**
	 * The name of the kernel to compile.
	 */
	const std::string kernel_name;

	/**
	 * The build options to use for the kernel.
	 */
	const std::string build_options;

	/**
	 * The OpenCL context we are working in.
	 */
	const cl_context context;

	/**
	 * The device to compile for.
	 */
	cl_device_id device;

	/**
	 * The files required for compilation.
	 */
	const std::vector<const char *> files;

	/**
	 * Print resource requirements of a kernel object.
	 *
	 * All information is dumped to the trace.
	 *
	 * @param kernel The kernel of which to query the information.
	 */
	void printResourceRequirements(const cl_kernel kernel) const;

	/**
	 * Generate an MD5 string uniquely identifying the OpenCL program binary.
	 */
	std::string generateMD5() const;

	void dumpBinary(cl_program program, std::string md5) const;

	cl_program loadBinary(std::string md5) const;
};

#endif /* _OPENCL_COMPILER_H_ */
