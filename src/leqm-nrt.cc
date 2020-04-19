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

#include "leqm-nrt.h"

#include <stdio.h>
#include <math.h>
#include <sndfile.h>
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

// Version 0.0.18 (C) Luca Trisciani 2011-2013, 2017-2018
// Tool from the DCP-Werkstatt Software Bundle

class Sum
{
public:
	void sum_samples(std::vector<double> const& input_samples, std::vector<double> const& c_input_samples, int nsamples)
	{
		_mutex.lock();
		_nsamples += nsamples;
		for (auto i = 0; i < nsamples; i++) {
			_sum  += input_samples[i];
			_csum += c_input_samples[i];
		}
		_mutex.unlock();
	}

	int nsamples() const
	{
		return _nsamples;
	}

	/*
	   How the final offset is calculated without reference to a test tone:
	   P0 is the SPL reference 20 uPa
	   Reference SPL is RMS ! So 85 SPL over 20 uPa is 10^4.25 x 0.000020 = 0.355655882 Pa (RMS),
	   but Peak value is 0.355655882 x sqr(2) = 0.502973372 that is 20 x log ( 0.502973372 / 0.000020) = 88.010299957
	   To that one has to add the 20 dB offset of the reference -20dBFS: 88.010299957 + 20.00 = 108.010299957

	   But ISO 21727:2004(E) ask for a reference level "measured using an average responding meter". So reference level is not 0.707, but 0.637 = 2/pi
	*/

	double mean() const
	{
		return pow(_sum / _nsamples, 0.500);
	}

	double cmean() const
	{
		return pow(_csum / _nsamples, 0.500);
	}

	double rms() const
	{
		return 20 * log10(mean()) + 108.010299957;
	}

	double leqm() const
	{
		return 20 * log10(cmean()) + 108.010299957;
	}

private:
	double _csum = 0.0; // convolved sum
	double _sum = 0.0; // flat sum
	int _nsamples = 0;
	std::mutex _mutex;
};

class Worker
{
public:
	Worker(std::vector<double> buffer, int nsamples, int nch, int npoints, std::vector<double> const& ir, Sum* sum, std::vector<double> chconf)
		: _buffer(buffer)
		, _nsamples(nsamples)
		, _nch(nch)
		, _npoints(npoints)
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

	int accumulate_ch(std::vector<double>& ch_accumulator, std::vector<double> const& input_channel, int nsamples) const
	{
		for (auto i = 0; i < nsamples; i++) {
			ch_accumulator[i] += input_channel[i];
		}
		return 0;
	}

	//rectify, square and sum
	int rectify(std::vector<double>& squared, std::vector<double> const& input_samples, int nsamples) const
	{
		for (auto i = 0; i < nsamples; i++) {
			squared[i] = powf(input_samples[i], 2);
		}
		return 0;

	}

	void convolv_buff(std::vector<double> const& sig_in, std::vector<double>& sig_out, std::vector<double> const& impresp, int sigin_dim, int impresp_dim) const
	{
		double sum = 0.0;
		for (int i = 0; i < sigin_dim; i++) {
			int m = i;
			for (int l = impresp_dim - 1; l >= 0; l--, m++) {
				if (m >= sigin_dim) {
					m -= sigin_dim;
				}
				sum += sig_in[m] * impresp[l];
			}
			sig_out[i] = sum;
			sum = 0.0;
		}
	}

	void process()
	{
		int const frames = _nsamples / _nch;

		std::vector<double> sum_and_square_buffer(frames);
		std::vector<double> c_sum_and_square_buffer(frames);
		std::vector<double> ch_sum_accumulator_norm(frames);
		std::vector<double> ch_sum_accumulator_conv(frames);

		for (int ch = 0; ch < _nch; ch++) {

			std::vector<double> normalized_buffer(frames);
			std::vector<double> convolved_buffer(frames);

			for (int n = ch, m = 0; n < _nsamples; n += _nch, m++) {
				// use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
				//so not normalized but calibrated
				normalized_buffer[m] = _buffer[n] * _chconf[ch]; //this scale amplitude according to specified calibration
			}

			//convolution
			convolv_buff(normalized_buffer, convolved_buffer, _ir, frames, _npoints * 2);
			//rectify, square und sum
			rectify(c_sum_and_square_buffer, convolved_buffer, frames);
			rectify(sum_and_square_buffer, normalized_buffer, frames);

			accumulate_ch(ch_sum_accumulator_norm, sum_and_square_buffer, frames);
			accumulate_ch(ch_sum_accumulator_conv, c_sum_and_square_buffer, frames);

		}

		_sum->sum_samples(ch_sum_accumulator_norm, ch_sum_accumulator_conv, frames);
	}

	std::vector<double> _buffer;
	int _nsamples;
	int _nch;
	int _npoints;
	std::vector<double> const& _ir;
	Sum* _sum;
	std::vector<double> _chconf;

	std::thread _thread;
};


void equalinterval2(double const freqsamples[], double const * freqresp, std::vector<double>& eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile);
std::vector<double> convert_log_to_linear(std::vector<double> const& in);
double convert_log_to_linear_single(double in);
double inputcalib (double dbdiffch);
std::vector<double> inversefft2(std::vector<double> const& eqfreqresp, int npoints);


Result calculate_file(
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

	auto result = calculate_function(
		[file](double* buffer, int64_t samples) -> int64_t {
			return sf_read_double(file, buffer, samples);
		},
		sf_info.channels,
		sf_info.samplerate,
		bitdepth,
		channel_corrections,
		buffer_size_ms,
		number_of_filter_interpolation_points,
		num_cpu
		);

	sf_close(file);

	return result;
}


Result calculate_function(
	std::function<int64_t (double*, int64_t)> get_audio_data,
	int channels,
	int sample_rate,
	int bits_per_sample,
	std::vector<double> channel_corrections,
	int buffer_size_ms,
	int number_of_filter_interpolation_points,
	int num_cpu
	)
{
	int constexpr origpoints = 21; //number of points in the standard CCIR filter

	std::vector<double> channel_conf_cal;

	//postprocessing parameters
	if (static_cast<int>(channel_corrections.size()) == channels) {
		for (auto i: channel_corrections) {
			channel_conf_cal.push_back(convert_log_to_linear_single(i));

		}
	} else if (channel_corrections.empty() && channels == 6) {
		double conf51[] = {0, 0, 0, 0, -3, -3};
		for (auto cind = 0; cind < channels; cind++) {
			channel_conf_cal.push_back(convert_log_to_linear_single(conf51[cind]));
		}
	} else {
		return {-100};
	}

	if ((sample_rate * buffer_size_ms) % 1000) {
		return -102;
	}

	int buffer_size_samples = (sample_rate * channels * buffer_size_ms) / 1000;
	std::vector<double> buffer(buffer_size_samples);


	//ISO 21727:2004(E)
	// M Weighting
	double const freqsamples[] = {31, 63, 100, 200, 400, 800, 1000, 2000, 3150, 4000, 5000, 6300, 7100, 8000, 9000, 10000, 12500, 14000, 16000, 20000, 31500};
	double const freqresp_db[] = {-35.5, -29.5, -25.4, -19.4, -13.4, -7.5, -5.6, 0.0, 3.4, 4.9, 6.1, 6.6, 6.4, 5.8, 4.5, 2.5, -5.6, -10.9, -17.3, -27.8, -48.3};

	std::vector<double> eqfreqresp_db(number_of_filter_interpolation_points);

	equalinterval2(freqsamples, freqresp_db, eqfreqresp_db, number_of_filter_interpolation_points, sample_rate, origpoints, bits_per_sample);
	auto eqfreqresp = convert_log_to_linear(eqfreqresp_db);

	auto ir = inversefft2(eqfreqresp, number_of_filter_interpolation_points);

	// read through the entire file

	Sum totsum;
	Result result;
	sf_count_t samples_read = 0;

	// Main loop through audio file

	int worker_id = 0;
	std::vector<std::shared_ptr<Worker>> worker_args;

	while ((samples_read = get_audio_data(buffer.data(), buffer_size_samples)) > 0) {
		worker_args.push_back(
			std::make_shared<Worker>(
				buffer,
				samples_read,
				channels,
				number_of_filter_interpolation_points,
				ir,
				&totsum,
				channel_conf_cal
				)
			);

		worker_id++;

		if (worker_id == num_cpu) {
			worker_id = 0;
			worker_args.clear();
		}
	}

	if (worker_id != 0) {
		for (int idxcpu = 0; idxcpu < worker_id; idxcpu++) { //worker_id is at this point one unit more than threads launched
			worker_args.clear();
		}
	}

	result.leq_nw = totsum.rms();
	result.leq_m = totsum.leqm();

	return result;
}


//the following is different from version 1 because interpolate between db and not linear. Conversion from db to lin must be done after.
//it is also different for the way it interpolates between DC and 31 Hz
// Pay attention that also arguments to the functions are changed
void equalinterval2(double const freqsamples[], double const freqresp_db[], std::vector<double>& eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile) {
	double freq;


	//calculate miminum attenuation depending on the bitdeph (minus one), that is −6.020599913 dB per bit in eccess to sign
	double dcatt = ((double) (bitdepthsoundfile - 1))*(-6.020599913) + 20.00; //in dB
	//double dcatt = -90.3;
	double pass = ((double) (samplingfreq >> 1)) / ((double) points);
	for (int ieq = 0, i = 0; ieq < points; ieq++) {
		freq = ieq*pass;

		if (freq == 0.0) {
			eqfreqresp[ieq] = dcatt;
		} else if (freq < freqsamples[0]) { // this has a lot of influence on final Leq(M) value
			eqfreqresp[ieq] = ((freqresp_db[0] - dcatt) / (freqsamples[0] - 0)) * freq + dcatt;
			//eqfreqresp[ieq] = freqresp_db[0]; // Is this meaningful? Shouldn't I interpolate between 0 Hz and 31 Hz? Otherwise for DC I have -35.5 dB
			continue;
		} else {

			if ((freq >= freqsamples[i]) && (freq < freqsamples[i+1])) {
				eqfreqresp[ieq] = ((freqresp_db[i+1] - freqresp_db[i])/(freqsamples[i+1] - freqsamples[i]))*(freq - freqsamples[i]) + freqresp_db[i];
			} else if (freq >=freqsamples[i+1]) {
				while(freq >= freqsamples[i+1]) {
					i++;
					if ((i + 1) >= origpoints) {
						break;
					}
				}
				if ((i+1) < origpoints) {
					eqfreqresp[ieq] = ((freqresp_db[i+1] - freqresp_db[i])/(freqsamples[i+1] - freqsamples[i]))*(freq- freqsamples[i]) + freqresp_db[i];
				} else {
					eqfreqresp[ieq] = ((1 - freqresp_db[i])/(((double) (samplingfreq >> 1)) - freqsamples[i]))*(freq- freqsamples[i]) + freqresp_db[i];
				}
			}
		}
	}
}


std::vector<double> convert_log_to_linear(std::vector<double> const& in)
{
	std::vector<double> out(in.size());
	for (auto i = 0U; i < in.size(); i++) {
		out[i] = powf(10, in[i] / 20.0);
	}
	return out;
}


double convert_log_to_linear_single(double in)
{
	return powf(10, in / 20.0f);
}


std::vector<double> inversefft2(std::vector<double> const& eqfreqresp, int npoints)
{
	std::vector<double> ir(npoints * 2);
	for (int n = 0; n < npoints; n++) {
		double parsum = 0.0;
		double partial = 0.0;

		for (int m = 1; m <= npoints - 1; m++) {
			partial = cos(2.0 * M_PI * m * ((n - (npoints * 2.0 - 1) / 2) / (npoints * 2.0)));
			parsum = parsum + eqfreqresp[m] * partial;
		}
		ir[n] = (eqfreqresp[0] + 2.0 * parsum) / (npoints * 2.0);
	}

	for (int n = 0; n < npoints; n++) {
		ir[npoints + n] = ir[npoints - (n + 1)];
	}

	return ir;
}

// scale input according to required calibration
// this could be different for certain digital cinema formats
double inputcalib(double dbdiffch)
{
	return pow(10, dbdiffch / 20);
}
