/* FM-index of alignment with gaps
    Copyright (C) 2015-2019  Hyunjoon Kim

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see https://www.gnu.org/licenses/.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_helper.hpp>
#include "csaa.hpp"
#include "gapinfo.hpp"

using namespace sdsl;
using namespace std;

void PrintUsage(int argc, char *argv[]) {
    std::cout << "usage : " << argv[0] << " saa_directory fma_directory num_strings pattern_file" << std::endl;
}

int main(int argc, char* argv[]) {
	if(argc < 5) {
		PrintUsage(argc, argv);
		return 0;
	}

	const int saa_dens = 32, isaa_dens = 32;
	string info_dir = string(argv[1]) + "/";
	string fma_dir = string(argv[2]) + "/";
	uint32_t num_string = atoi(argv[3]);
	string pattern_file = string(argv[4]);

	string gap_file = info_dir + "gapinfo.vec";
    string index_suffix = ".fma.gap.no_isaa." + to_string(saa_dens) + "." + to_string(isaa_dens);
    string index_file   = fma_dir + to_string(num_string)+index_suffix;

	int_vector<> gap_vector;
    load_from_file(gap_vector, gap_file);
    GapInfo gap_info(num_string, gap_vector);

    csaa<rrr_vector<127>, saa_dens, isaa_dens> fma(num_string, info_dir, gap_info.GetAlignLength());
    if (!load_from_file(fma, index_file)) {
        cout << "No index "<< index_file << " located. Building index now." << endl;
        construct(fma, info_dir, 1); // generate index
        store_to_file(fma, index_file); // save it
    }
    else{
        cout << "Load " << index_file << endl;
    }

	int patternNum = 0;
	clock_t begin, end;
	double count_elapsed = 0, locate_elapsed = 0;

	ifstream in(pattern_file);
	string query;
	while(getline(in, query)) {
		++patternNum;
		std::cout << "query " << patternNum << ": " << query << std::endl;
		begin = clock();
		// Counts the number of occurrences of a pattern in FMA.
		// \param fma           The FMA object.
		// \param query.begin() Iterator to the begin of the pattern(inclusive).
		// \param query.end()   Iterator to the end of the pattern(exclusive).
		// \return              The number of occurrences of the pattern in FMA.
		size_t occs = count(fma, query.begin(), query.end());
		end = clock();
		count_elapsed += double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "# of occurrences: " << occs << std::endl;
		if(occs > 0) {
			begin = clock();
			// Calculates all occurrences of a pattern in FMA.
			// \param fma           The FMA object.
			// \param query.begin() Iterator to the begin of the pattern(inclusive).
			// \param query.end()   Iterator to the end of the pattern(exclusive).
			// \return              A vector containing the occrrences of the pattern in FMA.
			auto locations = locate(fma, query.begin(), query.end());
			end = clock();
			locate_elapsed += double(end - begin) / CLOCKS_PER_SEC;

			for(int i = 0; i < occs; ++i) {
				int seqID = locations[2 * i] - 1;
                size_t indPos = gap_info.ConvertToInd(seqID, locations[2*i+1]);
				std::cout << "(" << seqID + 1 << ", " << indPos << "): ";

                size_t extract_begin = gap_info.ConvertToAlign(seqID, indPos);
                size_t extract_end = gap_info.ConvertToAlign(seqID, indPos + query.size() - 1);
                auto s   = extract(fma, locations[2*i], extract_begin, extract_end, query.size(), gap_info);
				std::cout << s << std::endl;
			}
		}
	}
	in.close();

	std::cout << "Pattern number: " << patternNum << std::endl;
	std::cout << "Count elapsed: " << count_elapsed << " sec" << std::endl;
	std::cout << "Count average: " << count_elapsed / patternNum << " sec" << std::endl;
	std::cout << "Locate elapsed: " << locate_elapsed << " sec" << std::endl;
	std::cout << "Locate average: " << locate_elapsed / patternNum << " sec" << std::endl;
	
	return 0;
}