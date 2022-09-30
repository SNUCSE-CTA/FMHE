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

#include "io.hpp"

int GetLength(std::string fname)
{
    std::ifstream in(fname, std::ifstream::ate | std::ifstream::binary);

    // File open error
    if(!in.is_open()) return -1;

    int ret = in.tellg();
    in.close();
    return ret;
}

void LoadSequence(std::string fname, char *s, int length)
{
    std::ifstream in(fname, std::ifstream::binary);

    // File open error
    if(!in.is_open()) return;

    in.read(s, length);
    in.close();
}

void LoadSequence(std::string fname, std::string& s)
{
    std::ifstream in(fname, std::ifstream::in);
    
    // File open error
    if(!in.is_open()) return;

    std::getline(in, s);
    in.close();
}

void StoreNonCommon(std::string fname, std::vector<int>& ncPos)
{
    std::ofstream out(fname, std::ofstream::out);

    // File open error
    if(!out.is_open()) return;

    for(int i = 0; i < ncPos.size() / 2; ++i) 
        out << ncPos[2 * i] << " " << ncPos[2 * i + 1] << "\n";
    out.close();
}

void StoreNonShared(std::string fname, 
            std::vector<int>& NSPos, 
            int alignLength, 
            int numSeq)
{
    int_vector<> NS(NSPos.size(), 0, bits::hi(alignLength) + 1);
    for(int i = 0; i < NSPos.size(); ++i) NS[i] = NSPos[i];

    store_to_file(NS, fname);
}