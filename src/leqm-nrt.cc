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

#ifdef _WIN32
#include <windows.h>
#elif __APPLE__
#include <sys/param.h>
#include <sys/sysctl.h
#endif

#include <vector>
#include <memory>
#include <iostream>
#include <thread>
#include <mutex>

// Version 0.0.18 (C) Luca Trisciani 2011-2013, 2017-2018
// Tool from the DCP-Werkstatt Software Bundle



// COMPILATION
// compile for DEBUG with gcc -g -DEBUG -lsndfile -lfftw3 -lm -lpthread -lrt -o leqm-nrt leqm-nrt.cpp
//otherwise  gcc -lsndfile -lm -lpthread -lrt -o leqm-nrt leqm-nrt.c

//#define DEBUG


class Sum
{
public:
	Sum()
	{
		csum = 0.0;
		sum = 0.0;
		nsamples = 0;
		cmean = 0.0;
		mean = 0.0; // Do I write anything here?
		leqm = 0.0;
		rms = 0.0;
	}

	void sumsamples(double * inputsamples, double * cinputsamples, int nsamples_)
	{
		_mutex.lock();
		nsamples += nsamples_;
		for (auto i = 0; i < nsamples_; i++) {
			sum  += inputsamples[i];
			csum += cinputsamples[i];
		}
		_mutex.unlock();
	}

	double csum; // convolved sum
	double sum; // flat sum
	int nsamples;
	double cmean; //convolved mean
	double mean;
	double leqm;
	double rms;

private:
	std::mutex _mutex;
};

class Worker
{
public:
	Worker(double* buffer, int buffer_size_samples, int nsamples, int nch, int npoints, double* ir, Sum* sum, double* chconf, int shorttermindex, double* shorttermarray, int leqm10flag)
		: _arg_buffer(new double[buffer_size_samples])
		, _nsamples(nsamples)
		, _nch(nch)
		, _npoints(npoints)
		, _ir(ir)
		, _sum(sum)
		, _chconf(chconf)
		, _shorttermindex(shorttermindex)
		, _shorttermarray(shorttermarray)
		, _leqm10flag(leqm10flag)
	{
		memcpy(_arg_buffer, buffer, nsamples * sizeof(double));
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
		delete[] _arg_buffer;
	}


private:
	double sumandshorttermavrg(double * channelaccumulator, int nsamples) const
	{
		double stsum = 0.0;
		for (int i=0; i < nsamples; i++) {
			stsum += channelaccumulator[i];

		}
		return stsum / (double) nsamples;
	}

	int accumulatech(double * chaccumulator, double * inputchannel, int nsamples) const
	{
		for (int i = 0; i < nsamples; i++) {
			chaccumulator[i] += inputchannel[i];
		}
		return 0;
	}

	//rectify, square and sum
	int rectify(double * squared, double * inputsamples, int nsamples) const
	{
		for (int i = 0; i < nsamples; i++) {
			squared[i] = (double) powf(inputsamples[i], 2);
		}
		return 0;

	}

	int convolv_buff(double * sigin, double * sigout, double * impresp, int sigin_dim, int impresp_dim) const
	{
		double sum = 0.0;
		for (int i = 0; i < sigin_dim; i++) {

			int m = i;
			for (int l = impresp_dim - 1; l >=0; l--,m++) {
				if (m >= sigin_dim) {
					m -= sigin_dim;
				}
				sum += sigin[m]*impresp[l];
			}
			sigout[i] = sum;
			sum=0.0;
		}
		return 0;
	}

	void process()
	{
		int const frames = _nsamples / _nch;

		auto sum_and_square_buffer = new double[frames];
		auto c_sum_and_square_buffer = new double[frames];
		auto ch_sum_accumulator_norm = new double[frames];
		auto ch_sum_accumulator_conv = new double[frames];

		for (int i = 0; i < frames; i++) {
			sum_and_square_buffer[i] = 0.0;
			c_sum_and_square_buffer[i] = 0.0;
			ch_sum_accumulator_norm[i] = 0.0;
			ch_sum_accumulator_conv[i] = 0.0;
		}

		for (int ch = 0; ch < _nch; ch++) {

			auto normalized_buffer = new double[frames];
			auto convolved_buffer = new double[frames];

			for (int n=ch, m= 0; n < _nsamples; n += _nch, m++) {
				// use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
				//so not normalized but calibrated
				normalized_buffer[m] = _arg_buffer[n] * _chconf[ch]; //this scale amplitude according to specified calibration
			}

			//convolution
			convolv_buff(normalized_buffer, convolved_buffer, _ir, frames, _npoints * 2);
			//rectify, square und sum
			rectify(c_sum_and_square_buffer, convolved_buffer, frames);
			rectify(sum_and_square_buffer, normalized_buffer, frames);

			accumulatech(ch_sum_accumulator_norm, sum_and_square_buffer, frames);
			accumulatech(ch_sum_accumulator_conv, c_sum_and_square_buffer, frames);

			delete[] normalized_buffer;
			delete[] convolved_buffer;

		} // loop through channels

		//Create a function for this also a tag so that the worker know if he has to do this or not

		if (_leqm10flag) {
			_shorttermarray[_shorttermindex] = sumandshorttermavrg(ch_sum_accumulator_conv, frames);
#ifdef DEBUG
			printf("%d: %.6f\n", _shorttermindex, _shorttermarray[_shorttermindex]);
#endif
		}

		_sum->sumsamples(ch_sum_accumulator_norm, ch_sum_accumulator_conv, frames);

		delete[] sum_and_square_buffer;
		delete[] c_sum_and_square_buffer;
		delete[] ch_sum_accumulator_norm;
		delete[] ch_sum_accumulator_conv;
	}

	double* _arg_buffer;
	int _nsamples;
	int _nch;
	int _npoints;
	double* _ir;
	Sum* _sum;
	double* _chconf;
	int _shorttermindex;
	double* _shorttermarray;
	int _leqm10flag;

	std::thread _thread;
};

int equalinterval( double * freqsamples, double * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints);
int equalinterval2( double freqsamples[], double * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile);
int convloglin(double * in, double * out, int points);
double convlinlog_single(double in);
double convloglin_single(double in);
double inputcalib (double dbdiffch);
int meanoverduration(Sum * oldsum);
void  inversefft1(double * eqfreqresp, double * ir, int npoints);
void  inversefft2(double * eqfreqresp, double * ir, int npoints);
void * worker_function(void * argfunc);
void logleqm(FILE * filehandle, double featuretimesec, Sum * oldsum);
void logleqm10(FILE * filehandle, double featuretimesec, double longaverage);

int calculate(int numcalread, SNDFILE* file, SF_INFO& sfinfo, char* soundfilename, double* channelconfcalvector, double* tempchcal, int leqmlog, FILE* leqmlogfile, int leqm10, FILE* leqm10logfile, struct timespec& starttime, int timing, double* shorttermaveragedarray, int bitdepth, int numbershortperiods, int buffersizems, int buffer_size_samples, int numCPU, int samplingfreq, int npoints, int origpoints, int leqnw);

int main(int argc, const char ** argv)
{
	int npoints = 64; // This value is low for precision. Calibration is done with 32768 point.
	int origpoints = 21; //number of points in the standard CCIR filter
	int samplingfreq; // this and the next is defined later taking it from sound file
	int bitdepth;
	// double normalizer;
	int timing = 0;
	struct timespec starttime;
	int fileopenstate = 0;
	int leqm10 = 0;
	int leqmlog = 0;
#if defined __unix__ || defined  __APPLE__
	int numCPU = sysconf(_SC_NPROCESSORS_ONLN) - 1;
#elif defined _WIN64 || defined _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	int numCPU = sysinfo.dwNumberOfProcessors - 1;
#endif

	double * channelconfcalvector = nullptr;
	printf("leqm-nrt  Copyright (C) 2011-2013, 2017-2018 Luca Trisciani\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it\nunder the GPL v3 licence.\nProgram will use 1 + %d slave threads.\n", numCPU);
	//SndfileHandle file;
	SNDFILE *file;
	file=NULL;
	SF_INFO sfinfo;
	FILE *leqm10logfile;
	leqm10logfile = NULL;
	FILE *leqmlogfile;
	leqmlogfile = NULL;
	int buffersizems = 850; //ISO 21727:2004 do not contain any indication, TASA seems to indicate 1000, p. 8
	int buffer_size_samples;
	double tempchcal[128];
	int numcalread = 0;
	double * shorttermaveragedarray = nullptr;
	int numbershortperiods = 0;
	int parameterstate = 0;
	int leqnw = 0;

	char soundfilename[1024];
	// This is a requirement of sndfile library, do not forget it.

	memset(&sfinfo, 0, sizeof(sfinfo));


	if (argc == 1)
	{ const char helptext[] = "Order of parameters is free.\nPossible parameters are:\n-convpoints <integer number> \tNumber of interpolation points for the filter. Default 64\n-numcpus <integer number> \tNumber of slave threads to speed up operation.\n-timing \t\t\tFor benchmarking speed.\n-leqnw\t Print out Leq without M Weighting\n-chconfcal <db correction> <db correction> <etc. so many times as channels>\n-logleqm10\n-logleqm\n-buffersize <milliseconds>\n";
		printf(helptext);
		printf("Please indicate a sound file to be processed.\n");
		return 0;
	}


	for (int in = 1; in < argc;) {

		if (!(strncmp(argv[in], "-", 1) == 0)) {
			if (fileopenstate == 0) {
				if(! (file = sf_open(argv[in], SFM_READ, &sfinfo))) {
					printf("Error while opening audio file, could not open  %s\n.", argv[in]);
					puts(sf_strerror(NULL));
					return 1;
				}

				strcpy(soundfilename, argv[in]);
				fileopenstate = 1;
				printf("Opened file: %s\n", argv[in]);
				printf("Sample rate: %d\n", sfinfo.samplerate);
				printf("Channels: %d\n", sfinfo.channels);
				printf("Format: %d\n", sfinfo.format);
				printf("Frames: %d\n", (int) sfinfo.frames);
				channelconfcalvector = (double *) malloc(sizeof(double) * sfinfo.channels);
				in++;
				continue;
			} else {
				free(channelconfcalvector);
				channelconfcalvector = NULL;
				return 0;
			}
		}


		if (strcmp(argv[in], "-chconfcal") == 0) {
			/* as the order of parameter is free I have to postpone
			   the check for consistency with the number of channels.
			   So first create a temporary array, whose number of element will be checked after
			   the parsing of the command line parameters is finished.
			   The calibration will be expressed in dB on the command line and converted to multiplier
			   here so that it can be stored as a factor in the channelconfcalvector.
			   */

			in++;
			for (;;)  {
				if (in < argc) {
					//if (!(strncmp(argv[in], "-", 1) == 0)) { //changed this to allow negative numbers
					if (!(strncmp(argv[in], "-", 1) == 0) || isdigit(argv[in][1])) {
						tempchcal[numcalread++]=atof(argv[in++]);
					} else break;
				} else break;

			} //for
			continue;
		}

		if (strcmp(argv[in], "-convpoints") == 0) {
			npoints = atoi(argv[in + 1]);
			in+=2;
			printf("Convolution points sets to %d.\n", npoints);
			continue;

		}


		if (strcmp(argv[in], "-version") == 0) {
			in++;
			printf("leqm-nrt version 0.18\n");
			continue;

		}
		if (strcmp(argv[in], "-numcpus") == 0) {
			numCPU= atoi(argv[in + 1]);
			in+=2;
			printf("Number of threads manually set to %d. Default is number of cores in the system minus one.\n", numCPU);
			continue;

		}
		if (strcmp(argv[in], "-timing") == 0) {
			timing = 1;
			in++;
			printf("Execution time will be measured.\n");
			continue;

		}

		if (strcmp(argv[in], "-logleqm10") == 0) {
			leqm10 = 1;
			in++;
			printf("Leq(M)10 data will be logged to the file leqm10.txt\n");
			continue;

		}
		if (strcmp(argv[in], "-logleqm") == 0) {
			leqmlog = 1;
			in++;
			printf("Leq(M) data will be logged to the file leqmlog.txt\n");
			continue;

		}

		if (strcmp(argv[in], "-leqnw") == 0) {
			leqnw = 1;
			in++;
			printf("Leq(nW) - unweighted -  will be outputted.\n");
			continue;

		}

		if (strcmp(argv[in], "-buffersize") == 0) {
			buffersizems = atoi(argv[in + 1]);
			in+=2;
			printf("Buffersize will be set to %d milliseconds.\n", buffersizems);
			continue;

		}

		if (parameterstate==0) {
			break;
		}
	}

	return calculate(numcalread, file, sfinfo, soundfilename, channelconfcalvector, tempchcal, leqmlog, leqmlogfile, leqm10, leqm10logfile, starttime, timing, shorttermaveragedarray, bitdepth, numbershortperiods, buffersizems, buffer_size_samples, numCPU, samplingfreq, npoints, origpoints, leqnw);
}

int calculate(int numcalread, SNDFILE* file, SF_INFO& sfinfo, char* soundfilename, double* channelconfcalvector, double* tempchcal, int leqmlog, FILE* leqmlogfile, int leqm10, FILE* leqm10logfile, struct timespec& starttime, int timing, double* shorttermaveragedarray, int bitdepth, int numbershortperiods, int buffersizems, int buffer_size_samples, int numCPU, int samplingfreq, int npoints, int origpoints, int leqnw)
{
	// Open audio file

	//postprocessing parameters
	if (numcalread == sfinfo.channels) {
		for (int cind = 0; cind < sfinfo.channels; cind++) {
			channelconfcalvector[cind] = convloglin_single(tempchcal[cind]);

		}
	}
	else if ((numcalread == 0) && (sfinfo.channels == 6)) {
		double conf51[6] = {0, 0, 0, 0, -3, -3};
		for (int cind = 0; cind < sfinfo.channels; cind++) {
			channelconfcalvector[cind] = convloglin_single(conf51[cind]);
		}

	} else {

		printf("Either you specified a different number of calibration than number of channels in the file or you do not indicate any calibration and the program cannot infer one from the number of channels. Please specify a channel calibration on the command line.\n");

		free(channelconfcalvector);
		channelconfcalvector = NULL;
		return 0;
	}




	if (leqm10) {
		char tempstring[1536];
		strcpy(tempstring, soundfilename);
		strcat(tempstring, ".leqm10.txt");
		leqm10logfile = fopen(tempstring, "w");
		if (leqm10logfile == NULL) {
			printf("Could not open file to write log leqm10 data!\n");
		}
	}




	if (leqmlog) {
		char tempstring[1536];
		strcpy(tempstring, soundfilename);
		strcat(tempstring, ".leqmlog.txt");
		leqmlogfile = fopen(tempstring, "w");
		if (leqmlogfile == NULL) {
			printf("Could not open file to write log leqm data!\n");
		}
	}


	if (timing) {
		clock_gettime(CLOCK_MONOTONIC, &starttime);
	}

	// reading to a double or float buffer with sndfile take care of normalization
	/*
   	static double  buffer[BUFFER_LEN]; // it seems this must be static. I don't know why
   	*/
	double * buffer;
	// buffer = new double [BUFFER_LEN];
	//buffer_size_samples = (sfinfo.samplerate*sfinfo.channels*buffersizems)/1000;
	if ((sfinfo.samplerate*buffersizems)%1000) {
		printf("Please fine tune the buffersize according to the sample rate\n");
		//close file
		// free memory
		// write a function to do that
		return 1;
	}

	buffer_size_samples = (sfinfo.samplerate*sfinfo.channels*buffersizems)/1000;
	buffer = new double[buffer_size_samples];

	samplingfreq = sfinfo.samplerate;

	if(leqm10) {

		//if duration < 10 mm exit

		double featdursec = sfinfo.frames / sfinfo.samplerate;
		if ((featdursec/60.0) < 10.0) {
			printf("The audio file is too short to measure Leq(m10).\n");
			return 0;
		}

		//how many short periods in overall duration
		int remainder = sfinfo.frames % (sfinfo.samplerate*buffersizems/1000);
		if (remainder == 0)  numbershortperiods = sfinfo.frames/(sfinfo.samplerate*buffersizems/1000);
		else  numbershortperiods = sfinfo.frames/(sfinfo.samplerate*buffersizems/1000) + 1;

		//allocate array
		shorttermaveragedarray = (double *) malloc(sizeof(*shorttermaveragedarray)*numbershortperiods);
	}


	//End opening audio file

	//ISO 21727:2004(E)
	// M Weighting
	double freqsamples[] = {31, 63, 100, 200, 400, 800, 1000, 2000, 3150, 4000, 5000, 6300, 7100, 8000, 9000, 10000, 12500, 14000, 16000, 20000, 31500};
	double freqresp_db[] = {-35.5, -29.5, -25.4, -19.4, -13.4, -7.5, -5.6, 0.0, 3.4, 4.9, 6.1, 6.6, 6.4, 5.8, 4.5, 2.5, -5.6, -10.9, -17.3, -27.8, -48.3};

	double * eqfreqresp_db;
	eqfreqresp_db = (double *) malloc(sizeof(*eqfreqresp_db)*npoints);

	double * eqfreqsamples;
	eqfreqsamples = (double *) malloc(sizeof(*eqfreqsamples)*npoints);
	double * eqfreqresp;
	eqfreqresp = (double *) malloc(sizeof(*eqfreqresp)*npoints);
	double * ir;
	ir = (double *) malloc(sizeof(*ir)*npoints*2);


	// And what to do for floating point sample coding?

	switch(sfinfo.format & SF_FORMAT_SUBMASK) {
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




	equalinterval2(freqsamples, freqresp_db, eqfreqsamples, eqfreqresp_db, npoints, samplingfreq, origpoints, bitdepth);
	convloglin(eqfreqresp_db, eqfreqresp, npoints);

#ifdef DEBUG
	for (int i=0; i < npoints; i++) {
		printf("%d\t%.2f\t%.2f\t%.6f\n", i, eqfreqsamples[i], eqfreqresp_db[i], eqfreqresp[i]);
	}
#endif

	inversefft2(eqfreqresp, ir, npoints);

	// read through the entire file

	auto totsum = new Sum();
	sf_count_t samples_read = 0;

	// Main loop through audio file

	int worker_id = 0;
	std::vector<std::shared_ptr<Worker>> worker_args;
	int staindex = 0; //shorttermarrayindex


	while((samples_read = sf_read_double(file, buffer, buffer_size_samples)) > 0) {
		worker_args.push_back(std::make_shared<Worker>(
					buffer, buffer_size_samples, samples_read, sfinfo.channels, npoints, ir, totsum, channelconfcalvector, leqm10 ? staindex++ : 0, leqm10 ? shorttermaveragedarray : 0, leqm10 ? 1 : 0)
				);
		worker_id++;

		if (worker_id == numCPU) {
			worker_id = 0;
			worker_args.clear();
			//simply log here your measurement it will be a multiple of your threads and your buffer
			if (leqmlog) {
				meanoverduration(totsum); //update leq(m) until now and log it
				logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) sfinfo.samplerate), totsum );
			} //endlog
		}


		//end while worker_id
		/// End looping cores
	} // main loop through file

	//here I should wait for rest worker (< numcpu)
	//but I need to dispose of thread id.
	if (worker_id != 0) { // worker_id = 0 means the number of samples was divisible through the number of cpus
		for (int idxcpu = 0; idxcpu < worker_id; idxcpu++) { //worker_id is at this point one unit more than threads launched
			worker_args.clear();
		}
		//also log here for a last value
		if (leqmlog) {
			meanoverduration(totsum); //update leq(m) until now and log it
			logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) sfinfo.samplerate), totsum );
		} //endlog
	}
	// mean of scalar sum over duration

	meanoverduration(totsum);
	if (leqnw) {
		printf("Leq(nW): %.4f\n", totsum->rms); // Leq(no Weighting)
	}
	printf("Leq(M): %.4f\n", totsum->leqm);

	if(timing) {
		struct timespec stoptime;
		long stoptimenanoseconds;
		long executionnanoseconds;
		clock_gettime(CLOCK_MONOTONIC, &stoptime);

		if (stoptime.tv_nsec < starttime.tv_nsec) {
			stoptimenanoseconds = 1000000000 + stoptime.tv_nsec;
		} else {
			stoptimenanoseconds = stoptime.tv_nsec;
		}
		executionnanoseconds = stoptimenanoseconds - starttime.tv_nsec;
		printf("Total execution time is %.6f seconds\n", ((double) stoptime.tv_sec) - ((double) starttime.tv_sec) + ((double) executionnanoseconds / 1000000000.00));
	}


	if (leqm10) {

		//Take the array with the short term accumulators
		double interval = 10.0;
		//create a rolling average according to rolling interval
		int rollint; // in short 10*60 = 600 sec 600/0.850


		//how many element of the array to consider for the rollint?
		//that is how many buffersizems in the interval - interval could be parameterized(?)
		double tempint = 60.0 * interval / (((double) buffersizems) /1000.0);
		rollint = (int) tempint;
		//dispose of the rest
		if (tempint - ((double) rollint) > 0) {
			rollint += 1;
		}
		//two loops
		//external loop
		int indexlong = 0;
		while(indexlong < (numbershortperiods - rollint)) {

			double accumulator = 0;
			//internal loop
			double averagedaccumulator = 0;
			for (int indexshort = 0; indexshort < rollint; indexshort++) {

				accumulator += shorttermaveragedarray[indexshort+indexlong];
			} //end internal loop
			averagedaccumulator = accumulator/((double) rollint);
			logleqm10(leqm10logfile, ((double) (indexlong+rollint)) * ((double) buffersizems / 1000.0), averagedaccumulator);
			indexlong++;
		} //end external loop

		fclose(leqm10logfile);
		free(shorttermaveragedarray);
		shorttermaveragedarray = NULL;
	}


	if (leqmlog) {

		fclose(leqmlogfile);
	}

	sf_close(file);

	free(eqfreqsamples);
	eqfreqsamples = NULL;
	free(eqfreqresp_db);
	eqfreqresp_db=NULL;
	free(eqfreqresp);
	eqfreqresp = NULL;
	free(ir);
	ir = NULL;
	free(channelconfcalvector);
	channelconfcalvector = NULL;

	free(totsum);
	totsum = NULL;
	delete[] buffer;
	buffer = nullptr;
}





	//to get impulse response frequency response at equally spaced intervals is needed

	int equalinterval( double * freqsamples, double  * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints) {
		double freq;
		// int findex = 0;
		// int rindex = 0;
		double pass = ((double) (samplingfreq >> 1)) / ((double) points);
		for (int ieq = 0, i = 0; ieq < points; ieq++) {
			freq = ieq*pass;
			eqfreqsamples[ieq] = freq;

			if ((freq == 0.0) || (freq < freqsamples[1])) {
				eqfreqresp[ieq] = freqresp[0];
				continue;
			} else {

				if ((freq >= freqsamples[i]) && (freq < freqsamples[i+1])) {
					eqfreqresp[ieq] = ((freqresp[i+1] - freqresp[i])/(freqsamples[i+1] - freqsamples[i]))*(freq - freqsamples[i]) + freqresp[i];
				} else if (freq >=freqsamples[i+1]) {
					while(freq >= freqsamples[i+1]) {
						i++;
						if ((i + 1) >= origpoints) {
							break;
						}
					}
					if ((i+1) < origpoints) {
						eqfreqresp[ieq] = ((freqresp[i+1] - freqresp[i])/(freqsamples[i+1] - freqsamples[i]))*(freq- freqsamples[i]) + freqresp[i];
					} else {
						eqfreqresp[ieq] = ((1 - freqresp[i])/(((double) (samplingfreq >> 1)) - freqsamples[i]))*(freq- freqsamples[i]) + freqresp[i];
					}
				}
			}
		}
		return 0;
	}





	//the following is different from version 1 because interpolate between db and not linear. Conversion from db to lin must be done after.
	//it is also different for the way it interpolates between DC and 31 Hz
	// Pay attention that also arguments to the functions are changed
	int equalinterval2( double freqsamples[], double  freqresp_db[], double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile) {
		double freq;


		//calculate miminum attenuation depending on the bitdeph (minus one), that is −6.020599913 dB per bit in eccess to sign
		double dcatt = ((double) (bitdepthsoundfile - 1))*(-6.020599913) + 20.00; //in dB
		//double dcatt = -90.3;
		double pass = ((double) (samplingfreq >> 1)) / ((double) points);
		for (int ieq = 0, i = 0; ieq < points; ieq++) {
			freq = ieq*pass;
			eqfreqsamples[ieq] = freq;

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
		return 0;
	}






	int convloglin(double * in, double * out, int points) {
		for (int i = 0; i < points; i++) {
			out[i] = powf(10, (in[i]/20.0));
		}

		return 0;
	}


	double convlinlog_single(double in) {
		double out;
		out = log(in)*20.0f;
		return out;
	}


	double convloglin_single(double in) {
		double out;
		out = powf(10, in/20.0f);
		return out;
	}

	// convolution


	void  inversefft2(double * eqfreqresp, double * ir, int npoints) {
		for (int n = 0; n < npoints; n++) {
			double parsum = 0.0;
			double partial = 0.0;

			for (int m = 1; m <= npoints -1; m++) {
				partial = cos(2.0*M_PI*((double) m)*( ( ((double) n) - ( ((double) npoints) * 2.0 -1 ) / 2 ) / ( ((double) npoints) * 2.0) ));
				parsum = parsum + eqfreqresp[m]*partial;
			}
			ir[n] = (eqfreqresp[0] + 2.0 * parsum)/((double) npoints * 2.0);
#ifdef DEBUG
			printf("%.4f\n", ir[n]);
#endif
		}
		for (int n = 0; n < npoints; n++) {
			ir[npoints+n] = ir[npoints-(n + 1)];
#ifdef DEBUG
			printf("%.4f\n", ir[npoints+n]);
#endif
		}


	}

	// scale input according to required calibration
	// this could be different for certain digital cinema formats
	double inputcalib (double dbdiffch) {

		double coeff = pow(10, dbdiffch/20);
		return coeff;

	}

	int initbuffer(double * buffertoinit, int nsamples) {
		for (int i = 0; i < nsamples; i++) {
			buffertoinit[i] = 0.0;

		}
		return 0;
	}

	int meanoverduration(Sum * oldsum) {
		oldsum->mean = pow(oldsum->sum / ((double) oldsum->nsamples), 0.500);
		oldsum->cmean = pow(oldsum->csum / ((double) oldsum->nsamples), 0.500);
		oldsum->rms = 20*log10(oldsum->mean) + 108.010299957;
		oldsum->leqm = 20*log10(oldsum->cmean) +  108.010299957;//

		/*
		   How the final offset is calculated without reference to a test tone:
		   P0 is the SPL reference 20 uPa
		   Reference SPL is RMS ! So 85 SPL over 20 uPa is 10^4.25 x 0.000020 = 0.355655882 Pa (RMS),
		   but Peak value is 0.355655882 x sqr(2) = 0.502973372 that is 20 x log ( 0.502973372 / 0.000020) = 88.010299957
		   To that one has to add the 20 dB offset of the reference -20dBFS: 88.010299957 + 20.00 = 108.010299957
		   */
		/*But ISO 21727:2004(E) ask for a reference level "measured using an average responding meter". So reference level is not 0.707, but 0.637 = 2/pi
		*/
		return 0;
	}

	void logleqm(FILE * filehandle, double featuretimesec, Sum * oldsum) {

		fprintf(filehandle, "%.4f", featuretimesec);
		fprintf(filehandle, "\t");
		fprintf(filehandle, "%.4f\n", oldsum->leqm);


	}

	void logleqm10(FILE * filehandle, double featuretimesec, double longaverage) {
		double leqm10 = 20*log10(pow(longaverage, 0.500)) +  108.010299957;
		fprintf(filehandle, "%.4f", featuretimesec);
		fprintf(filehandle, "\t");
		fprintf(filehandle, "%.4f\n", leqm10);

	}
