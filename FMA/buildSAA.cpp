/* FM-index of alignment with gaps
    Copyright (C) 2015-2019  Seunghwan Min

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

#include "vcf.hpp"
#include "gapinfo.hpp"
#include "saa.hpp"

void PrintUsage(int argc, char *argv[])
{
    std::cout << "usage : " << argv[0] << " sampling_rate reference_file vcf_file output_directory [num_threads]" << std::endl;
}

int main(int argc, char *argv[])
{
    if(argc < 5)
    {
        PrintUsage(argc, argv);
        return 0;
    }

    SAAConfig config;
    config.samplingDensity = atoi(argv[1]);
    std::string refPath(argv[2]);
    std::string vcfPath(argv[3]);
    std::string outDir = std::string(argv[4]) + "/";
    if(argc >= 6) config.numThread = atoi(argv[5]);
    else config.numThread = 1;

    VCF vcf;
    vcf.Parsing(vcfPath, refPath, outDir);

    int numSeq = vcf.GetNumSeq();
    SAA saa(numSeq, config, outDir, outDir);
    saa.BuildSAA();
    

    return 0;
}