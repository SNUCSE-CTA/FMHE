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

#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <ctime>
#include <fstream>
#include "csaa.hpp"
#include "gapinfo.hpp"
using namespace sdsl;
using namespace std;

void PrintUsage(int argc, char *argv[])
{
    std::cout << "usage : " << argv[0] << " saa_directory num_strings output_directory" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc < 4) 
    {
        PrintUsage(argc, argv);
        return 0;
    }

    const int saa_dens = 32, isaa_dens = 32;
    string info_dir = string(argv[1]) + "/";
    uint32_t num_string = atoi(argv[2]);
    string output_dir = string(argv[3]) + "/";
    
    string gap_file = info_dir + "gapinfo.vec";
    string index_suffix = ".fma.gap.no_isaa." + to_string(saa_dens) + "." + to_string(isaa_dens);
    string index_file   = output_dir + to_string(num_string)+index_suffix;
    
    int_vector<> gap_vector;
    load_from_file(gap_vector, gap_file);
    GapInfo gap_info(num_string, gap_vector);
    
    cout << "total align length: " << gap_info.GetAlignLength() << endl;
    csaa<rrr_vector<127>, saa_dens, isaa_dens> fma(num_string, info_dir, gap_info.GetAlignLength());
    if (!load_from_file(fma, index_file)) {
        cout << "No index "<< index_file << " located. Building index now." << endl;
        construct(fma, info_dir, 1); // generate index
        store_to_file(fma, index_file); // save it
    }
    else{
        cout << "Load " << index_file << endl;
    }
    
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fma) << " MiB." << endl;
    double saa_sample_size      = size_in_mega_bytes(fma.saa_sample);
    double marked_size          = size_in_mega_bytes(fma.saa_sample.marked);
    double rank_marked_size     = size_in_mega_bytes(fma.saa_sample.rank_marked);
    double strv_size            = size_in_mega_bytes(fma.saa_sample.strv);
    double ns_isaa_sample_size  = size_in_mega_bytes(fma.ns_isaa_sample);
    double index_table_size     = size_in_mega_bytes(fma.index_table);
    cout << "saa_sample requires " << saa_sample_size << " MiB." << endl;
    cout << "In saa_sample, saa_sample.marked requires " << marked_size << " MiB." << endl;
    cout << "In saa_sample, saa_sample.rank_marked requires " << rank_marked_size << " MiB." << endl; 
    cout << "In saa_sample, saa_sample.strv requires " << strv_size << " MiB." << endl; 
    cout << "In saa_sample, saa_sample.sample requires " << saa_sample_size - marked_size - rank_marked_size - strv_size << " MiB." << endl;
    
    cout << "isaa_sample requires " << ns_isaa_sample_size + index_table_size << " MiB." << endl;
    cout << "gap_vector requires " << size_in_mega_bytes(fma.gap_vector) << " MiB." << endl; 
    cout << "saa sampling count is " << fma.saa_sample.size() << endl; 
    fma.print();
    cout << "the number of vectors in strv: " << fma.saa_sample.strv.size()/num_string << endl;
    cout << "fma size: " << fma.size()<< endl;
    cout << "alignment length: " << fma.text_size << endl;
}


