#include <string>
#include <vector>
#include <functional>

struct Result
{
	Result(int status_)
		: status(status_)
	{}

	Result(double leq_m_, double leq_nw_)
		: status(0)
		, leq_m(leq_m_)
		, leq_nw(leq_nw_)
	{}

	/** 0 on success, or
	 *
	 * -100: Either channel_corrections contained a different number of
	 *  calibrations than; number of channels in the file or it was empty and  the
	 *  program cannot infer one from the number of channels. Please specify a
	 *  values in channel_corrections.
	 *
	 * -101: Failed to open the sound file.
	 *
	 * -102: buffer_size_ms is not an integer number of samples at the sound file's rate.
	 */
	int status;

	double leq_m;
	double leq_nw;
};


Result calculate_file(
	std::string sound_filename,
	std::vector<double> channel_corrections,
	int buffer_size_ms,
	int number_of_filter_interpolation_points,
	int num_cpu
	);


Result calculate_function(
	std::function<int64_t (double*, int64_t)> get_audio_data,
	int channels,
	int sample_rate,
	int bits_per_sample,
	std::vector<double> channel_corrections,
	int buffer_size_ms,
	int number_of_filter_interpolation_points,
	int num_cpu
	);


double convert_log_to_linear_single(double in);
