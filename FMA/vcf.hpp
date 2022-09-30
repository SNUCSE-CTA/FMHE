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
/*! \file vcf.hpp
    \brief vcf.hpp contains a class for parsing the vcf file.
    \author Seunghwan Min
    \date 2017.1
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "gapinfo.hpp"
#include "io.hpp"

#pragma once

const unsigned int BUFFERSIZE = 1024 * 1024;

class VCF
{
private:
    struct VCFInfo
    {
        int Chrom;
        int Pos;
        std::string ID;
        std::string Ref;
        std::string Alt;
        std::vector<std::string> Alts;
        std::string Qual;
        std::string Filter;
        std::string Info;
        std::string Format;
        std::vector<std::string> AllelInfos;
    };

    int numChrom;
    int numSeq;

    char delim;
    std::ifstream vcf;
    std::string line;
    VCFInfo curInfo;


public:
    VCF() : numChrom(0), numSeq(0) {}
    ~VCF()
    {
        if(vcf.is_open()) vcf.close();
    }

    inline int GetNumSeq() {return numSeq;}

    void OpenFile(const std::string fname);
    void ReadHeader();
    void ReadNext(bool first);
    void Parsing(const std::string& vcfPath, const std::string& refPath, const std::string& outDir);

    

};