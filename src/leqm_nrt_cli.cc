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



// COMPILATION
// compile for DEBUG with gcc -g -DEBUG -lsndfile -lfftw3 -lm -lpthread -lrt -o leqm-nrt leqm-nrt.cpp
//otherwise  gcc -lsndfile -lm -lpthread -lrt -o leqm-nrt leqm-nrt.c

//#define DEBUG

int main(int argc, const char ** argv)
{
	int number_of_filter_interpolation_points = 64; // This value is low for precision. Calibration is done with 32768 point.
	int num_cpu = std::thread::hardware_concurrency() - 1;

	printf("leqm-nrt  Copyright (C) 2011-2013, 2017-2018 Luca Trisciani\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it\nunder the GPL v3 licence.\nProgram will use 1 + %d slave threads.\n", num_cpu);
	int buffer_size_ms = 850; //ISO 21727:2004 do not contain any indication, TASA seems to indicate 1000, p. 8
	std::vector<double> channel_corrections;
	int parameterstate = 0;
	bool display_leqnw = false;
	std::string sound_filename;

	if (argc == 1)
	{ const char helptext[] = "Order of parameters is free.\nPossible parameters are:\n-convpoints <integer number> \tNumber of interpolation points for the filter. Default 64\n-numcpus <integer number> \tNumber of slave threads to speed up operation.\n-timing \t\t\tFor benchmarking speed.\n-leqnw\t Print out Leq without M Weighting\n-chconfcal <db correction> <db correction> <etc. so many times as channels>\n-logleqm10\n-logleqm\n-buffersize <milliseconds>\n";
		printf(helptext);
		printf("Please indicate a sound file to be processed.\n");
		return 0;
	}

	for (int in = 1; in < argc;) {

		if (!(strncmp(argv[in], "-", 1) == 0)) {
			sound_filename = argv[in];
			in++;
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
					if (!(strncmp(argv[in], "-", 1) == 0) || isdigit(argv[in][1])) {
						channel_corrections.push_back(leqm_nrt::convert_log_to_linear_single(atof(argv[in++])));
					} else break;
				} else break;

			} //for
			continue;
		}

		if (strcmp(argv[in], "-convpoints") == 0) {
			number_of_filter_interpolation_points = atoi(argv[in + 1]);
			in+=2;
			printf("Convolution points sets to %d.\n", number_of_filter_interpolation_points);
			continue;

		}


		if (strcmp(argv[in], "-version") == 0) {
			in++;
			printf("leqm-nrt version 0.18\n");
			continue;

		}
		if (strcmp(argv[in], "-numcpus") == 0) {
			num_cpu= atoi(argv[in + 1]);
			in+=2;
			printf("Number of threads manually set to %d. Default is number of cores in the system minus one.\n", num_cpu);
			continue;

		}

		if (strcmp(argv[in], "-leqnw") == 0) {
			display_leqnw = true;
			in++;
			printf("Leq(nW) - unweighted -  will be outputted.\n");
			continue;

		}

		if (strcmp(argv[in], "-buffersize") == 0) {
			buffer_size_ms = atoi(argv[in + 1]);
			in+=2;
			printf("Buffersize will be set to %d milliseconds.\n", buffer_size_ms);
			continue;

		}

		if (parameterstate==0) {
			break;
		}
	}

	auto result = leqm_nrt::calculate_file(sound_filename, channel_corrections, buffer_size_ms, number_of_filter_interpolation_points, num_cpu);

	if (display_leqnw) {
		printf("Leq(nW): %.4f\n", result.leq_nw); // Leq(no Weighting)
	}
	printf("Leq(M): %.4f\n", result.leq_m);

	return result.status;
}

