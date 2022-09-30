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
/*! \file io.hpp
    \brief io.hpp contains an implementation of load and store functions.
    \author Seunghwan Min
    \date 2017.1
*/

#include <fstream>
#include <string>
#include <map>

#include <sdsl/io.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/construct.hpp>

#pragma once

using namespace sdsl;

int GetLength(std::string fname);

void LoadSequence(std::string fname, char *s, int length);

void LoadSequence(std::string fname, std::string& s);

void StoreNonCommon(std::string fname, std::vector<int>& ncPos);

void StoreNonShared(std::string fname, 
            std::vector<int>& NSPos, 
            int alignLength, 
            int numSeq);
