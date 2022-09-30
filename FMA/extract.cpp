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
    std::cout << "usage : " << argv[0] << " saa_directory fma_directory num_strings pattern_file pattern_length" << std::endl;
}

int main(int argc, char* argv[]) {
	if(argc < 6) {
		PrintUsage(argc, argv);
		return 0;
	}

	const int saa_dens = 32, isaa_dens = 32;
	string info_dir = string(argv[1]) + "/";
	string fma_dir = string(argv[2]) + "/";
	uint32_t num_string = atoi(argv[3]);
	string pattern_file = string(argv[4]);
	uint32_t pattern_size = atoi(argv[5]);

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

	int patternNum = 0, seqNum, offset;
	clock_t begin, end;
	double extract_elapsed = 0;

	ifstream in(pattern_file);
	while(in >> seqNum >> offset) {
		++patternNum;
		size_t extract_begin = gap_info.ConvertToAlign(seqNum - 1, offset);
		size_t extract_end = gap_info.ConvertToAlign(seqNum - 1, offset + pattern_size - 1);

		begin = clock();
		// Retrieve the substring T_i[begin..end] of the original i-th text T_i to text[0..end-begin+1].
		// \param fma           The FMA object.
		// \param seqNum        ID of the sequence which should be extracted (seqNum-th sequence).
		// \param extract_begin Position of the first character which should be extracted (inclusive).
		// \param extract_end   Position of the last character which should be extracted (inclusive).
		// \param pattern_size  The length of the pattern that should be extracted .
		// \param gap_info      The GapInfo object.
		// \return s            The object holding the extracted text. 
		auto s = extract(fma, seqNum, extract_begin, extract_end, pattern_size, gap_info);
		end = clock();
		extract_elapsed += double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "query " << patternNum << " (" << seqNum << ", " << offset << "): " << s << std::endl;
	}
	in.close();

	std::cout << "Pattern number: " << patternNum << std::endl;
	std::cout << "Extract elapsed: " << extract_elapsed << " sec" << std::endl;
	std::cout << "Extract average: " << extract_elapsed / patternNum << " sec" << std::endl;
	
	return 0;
}