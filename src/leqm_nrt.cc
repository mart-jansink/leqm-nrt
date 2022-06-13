/*
    leqm-nrt is a  non-real-time implementation
    of Leq(M) measurement according to ISO 21727:2004(E)
    "Cinematography -- Method of measurement of perceived
    loudness of motion-picture audio material"

    Copyright (C) 2011-2013, 2017-2018 Luca Trisciani

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "leqm_nrt.h"

#include <stdio.h>
#include <math.h>
#ifdef LEQM_NRT_WITH_LIBSNDFILE
#include <sndfile.h>
#endif
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <iso646.h>

#include <vector>
#include <memory>
#include <iostream>
#include <thread>
#include <mutex>
#include <functional>

using namespace leqm_nrt;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Version 0.0.18 (C) Luca Trisciani 2011-2013, 2017-2018
// Tool from the DCP-Werkstatt Software Bundle

namespace leqm_nrt {

class Worker
{
public:
	Worker(std::vector<double> buffer, int nsamples, int nch, std::vector<double> const& ir, Sum* sum, std::vector<double> chconf)
		: _buffer(buffer)
		, _nsamples(nsamples)
		, _nch(nch)
		, _ir(ir)
		, _sum(sum)
		, _chconf(chconf)
	{
		_thread = std::thread(&Worker::process, this);
	}

	Worker(Worker& other) = delete;
	bool operator=(Worker&) = delete;
	Worker(Worker&& other) = delete;
	bool operator=(Worker&&) = delete;

	~Worker()
	{
		try {
			_thread.join();
		} catch (...)
		{
		}
	}


private:
	double sum_and_short_term_avrg(std::vector<double>const & channel_accumulator, int nsamples) const
	{
		double stsum = 0.0;
		for (auto i = 0; i < nsamples; i++) {
			stsum += channel_accumulator[i];

		}
		return stsum / nsamples;
	}

	void accumulate_ch(std::vector<double>& ch_accumulator, std::vector<double> const& input_channel, int nsamples) const
	{
		for (auto i = 0; i < nsamples; i++) {
			ch_accumulator[i] += input_channel[i];
		}
	}

	std::vector<double> rectify(std::vector<double> const& input) const
	{
		std::vector<double> squared;
		for (auto i: input) {
			squared.push_back(pow(i, 2));
		}
		return squared;
	}

	std::vector<double> convolve(std::vector<double> const& signal, std::vector<double> const& ir) const
	{
		std::vector<double> result(signal.size());
		double sum = 0.0;
		for (auto i = 0U; i < signal.size(); i++) {
			auto m = i;
			for (int l = ir.size()- 1; l >= 0; l--, m++) {
				if (m >= signal.size()) {
					m -= signal.size();
				}
				sum += signal[m] * ir[l];
			}
			result[i] = sum;
			sum = 0.0;
		}
		return result;
	}

	void process()
	{
		int const frames = _nsamples / _nch;

		std::vector<double> ch_sum_accumulator_norm(frames);
		std::vector<double> ch_sum_accumulator_conv(frames);

		for (int ch = 0; ch < _nch; ch++) {

			std::vector<double> normalized_buffer(frames);

			for (int n = ch, m = 0; n < _nsamples; n += _nch, m++) {
				// use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
				//so not normalized but calibrated
				normalized_buffer[m] = _buffer[n] * _chconf[ch]; //this scale amplitude according to specified calibration
			}

			auto convolved_buffer = convolve(normalized_buffer, _ir);
			accumulate_ch(ch_sum_accumulator_norm, rectify(normalized_buffer), frames);
			accumulate_ch(ch_sum_accumulator_conv, rectify(convolved_buffer), frames);
		}

		_sum->sum_samples(ch_sum_accumulator_norm, ch_sum_accumulator_conv, frames);
	}

	std::vector<double> _buffer;
	int _nsamples;
	int _nch;
	std::vector<double> const& _ir;
	Sum* _sum;
	std::vector<double> _chconf;

	std::thread _thread;
};

}


//the following is different from version 1 because interpolate between db and not linear. Conversion from db to lin must be done after.
//it is also different for the way it interpolates between DC and 31 Hz
// Pay attention that also arguments to the functions are changed
static std::vector<double> equalinterval2(int points, int samplingfreq, int bitdepthsoundfile)
{
	//ISO 21727:2004(E)
	// M Weighting
	double const frequencies[] = {31, 63, 100, 200, 400, 800, 1000, 2000, 3150, 4000, 5000, 6300, 7100, 8000, 9000, 10000, 12500, 14000, 16000, 20000, 31500};
	double const db[] = {-35.5, -29.5, -25.4, -19.4, -13.4, -7.5, -5.6, 0.0, 3.4, 4.9, 6.1, 6.6, 6.4, 5.8, 4.5, 2.5, -5.6, -10.9, -17.3, -27.8, -48.3};
	int constexpr points_in_standard_ccir_filter = 21;

	std::vector<double> freq_response(points);
	//calculate miminum attenuation depending on the bitdeph (minus one), that is −6.020599913 dB per bit in eccess to sign
	double const dcatt = ((double) (bitdepthsoundfile - 1))*(-6.020599913) + 20.00; //in dB
	//double dcatt = -90.3;
	double const pass = ((double) (samplingfreq >> 1)) / ((double) points);
	for (int ieq = 0, i = 0; ieq < points; ieq++) {
		double const freq = ieq * pass;

		if (freq == 0.0) {
			freq_response[ieq] = dcatt;
		} else if (freq < frequencies[0]) { // this has a lot of influence on final Leq(M) value
			freq_response[ieq] = ((db[0] - dcatt) / (frequencies[0] - 0)) * freq + dcatt;
			continue;
		} else {

			if (freq >= frequencies[i] && freq < frequencies[i+1]) {
				freq_response[ieq] = ( (db[i+1] - db[i]) / (frequencies[i+1] - frequencies[i])) * (freq - frequencies[i]) + db[i];
			} else if (freq >= frequencies[i+1]) {
				while (freq >= frequencies[i+1]) {
					i++;
					if ((i + 1) >= points_in_standard_ccir_filter) {
						break;
					}
				}
				if ((i+1) < points_in_standard_ccir_filter) {
					freq_response[ieq] = ( (db[i+1] - db[i]) / (frequencies[i+1] - frequencies[i])) * (freq - frequencies[i]) + db[i];
				} else {
					freq_response[ieq] = ( (1 - db[i]) / ((samplingfreq >> 1) - frequencies[i])) * (freq - frequencies[i]) + db[i];
				}
			}
		}
	}
	return freq_response;
}


static std::vector<double> convert_log_to_linear(std::vector<double> const& in)
{
	std::vector<double> out(in.size());
	for (auto i = 0U; i < in.size(); i++) {
		out[i] = powf(10, in[i] / 20.0);
	}
	return out;
}


static std::vector<double> inverse_fft(std::vector<double> const& freq_response)
{
	int const npoints = freq_response.size();
	std::vector<double> ir(npoints * 2);
	for (int n = 0; n < npoints; n++) {
		double parsum = 0.0;
		double partial = 0.0;

		for (int m = 1; m <= npoints - 1; m++) {
			partial = cos(2.0 * M_PI * m * ((n - (npoints * 2.0 - 1) / 2) / (npoints * 2.0)));
			parsum = parsum + freq_response[m] * partial;
		}
		ir[n] = (freq_response[0] + 2.0 * parsum) / (npoints * 2.0);
	}

	for (int n = 0; n < npoints; n++) {
		ir[npoints + n] = ir[npoints - (n + 1)];
	}

	return ir;
}


#ifdef LEQM_NRT_WITH_LIBSNDFILE
Result leqm_nrt::calculate_file(
	std::string sound_filename,
	std::vector<double> channel_corrections,
	int buffer_size_ms,
	int number_of_filter_interpolation_points,
	int num_cpu
	)
{
	SF_INFO sf_info;
	sf_info.format = 0;
	SNDFILE* file = sf_open(sound_filename.c_str(), SFM_READ, &sf_info);
	if (!file) {
		return {-101};
	}

	int bitdepth = 0;
	switch(sf_info.format & SF_FORMAT_SUBMASK) {
		// all signed bitdepth
		case 0x0001:
			bitdepth = 8;
			break;
		case 0x0002:
			bitdepth = 16;
			break;
		case 0x0003:
			bitdepth = 24;
			break;
		case 0x0004:
			bitdepth = 32;
			break;
		default:
			printf("No known bitdepth! Exiting ...\n");
			return -1;
	}

	Calculator calculator(
		sf_info.channels,
		sf_info.samplerate,
		bitdepth,
		channel_corrections,
		buffer_size_ms,
		number_of_filter_interpolation_points,
		num_cpu
		);

	while (true) {
		std::vector<double> buffer(4096);
		auto read = sf_read_double(file, buffer.data(), buffer.size());
		if (read <= 0) {
			break;
		}

		buffer.resize(read);
		calculator.add(buffer);
	}

	sf_close(file);

	return {calculator.leq_m(), calculator.leq_nw()};
}
#endif


std::vector<double> calculate_ir(double number_of_filter_interpolation_points, int sample_rate, int bits_per_sample)
{
	return inverse_fft(convert_log_to_linear(equalinterval2(number_of_filter_interpolation_points, sample_rate, bits_per_sample)));
}


/** @return List of default channel corrections for the given number of channels,
 *  or an empty vector if we can't offer a default.
 */
std::vector<double> default_channel_corrections(int channels)
{
	std::vector<double> corr;

	if (channels == 6) {
		double conf51[] = {0, 0, 0, 0, -3, -3};
		for (auto cind = 0; cind < channels; cind++) {
			corr.push_back(convert_log_to_linear_single(conf51[cind]));
		}
	}

	return corr;
}


double leqm_nrt::convert_log_to_linear_single(double in)
{
	return powf(10, in / 20.0f);
}


Calculator::Calculator(
		int channels,
		int sample_rate,
		int bits_per_sample,
		std::vector<double> channel_corrections,
		int buffer_size_ms,
		int number_of_filter_interpolation_points,
		int num_cpu
	)
	: _channels(channels)
	, _channel_corrections(channel_corrections)
	, _num_cpu(num_cpu)
{
	if ((sample_rate * buffer_size_ms) % 1000) {
		throw BadBufferSizeError();
	}

	if (static_cast<int>(channel_corrections.size()) != channels) {
		channel_corrections = default_channel_corrections(channels);
	}
	if (static_cast<int>(channel_corrections.size()) != channels) {
		throw BadChannelCorrectionsError();
	}

	_ir = calculate_ir(number_of_filter_interpolation_points, sample_rate, bits_per_sample);
	_buffer.resize((sample_rate * channels * buffer_size_ms) / 1000);
}


void Calculator::process_buffer()
{
	if (_buffer_free_offset == 0) {
		return;
	}

	_workers.push_back(
			std::make_shared<Worker>(
				_buffer,
				_buffer_free_offset,
				_channels,
				_ir,
				&_sum,
				_channel_corrections
				)
			);

	if (static_cast<int>(_workers.size()) == _num_cpu) {
		_workers.clear();
	}

	_buffer_free_offset = 0;
}


void Calculator::add(std::vector<double> samples)
{
	size_t samples_offset = 0;

	while (samples_offset < samples.size()) {
		/* Copy some data into our buffer */
		auto to_copy = std::min(samples.size() - samples_offset, _buffer.size() - _buffer_free_offset);
		memcpy(_buffer.data() + _buffer_free_offset, samples.data() + samples_offset, to_copy * sizeof(double));
		samples_offset += to_copy;
		_buffer_free_offset += to_copy;

		/* Process the buffer if it's full */
		if (_buffer_free_offset == _buffer.size()) {
			process_buffer();
		}
	}
}


double Calculator::leq_m()
{
	process_buffer();
	_workers.clear();

	return _sum.leqm();
}


double Calculator::leq_nw()
{
	process_buffer();
	_workers.clear();

	return _sum.rms();
}
