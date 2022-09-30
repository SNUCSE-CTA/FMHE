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

#include "custom_int_vector.hpp"

custom_int_vector::custom_int_vector(char *file_name)
{
    read_file(file_name);
}

custom_int_vector::~custom_int_vector()
{
    if(ifs.is_open())
    {
        ifs.close();
        delete[] prev_buf;
        delete[] curr_buf;
    }
}

void custom_int_vector::read_file(char *file_name)
{
    read_file(std::string(file_name));
}

void custom_int_vector::read_file(std::string file_name)
{
    if(ifs.is_open())
    {
        ifs.close();
        delete[] prev_buf;
        delete[] curr_buf;
    }

    fname = file_name;
    ifs.open(fname, std::ifstream::binary);
    read_header();
    
    num_integers = (1 << 22) * (64 / m_width);
    buf_size = (1 << 16) * (64 / m_width) * m_width;
    
    prev_buf = new uint64_t[buf_size];
    curr_buf = new uint64_t[buf_size];
    
    prev_idx = curr_idx = num_buf_read = 0;
}

void custom_int_vector::read_header()
{
    ifs.read((char*)&m_size, sizeof(m_size));
    ifs.read((char*)&m_width, sizeof(m_width));
}

void custom_int_vector::read_buffer()
{
    std::swap(prev_buf, curr_buf);
    prev_idx = curr_idx;
    if(num_buf_read + buf_size < (capacity()>>6))
    {
        ifs.read((char*)curr_buf, buf_size * sizeof(uint64_t));
        num_buf_read += buf_size;
        curr_idx += num_integers;
    }
    else
    {
        ifs.read((char*)curr_buf, ((capacity()>>6)-num_buf_read) * sizeof(uint64_t));
        num_buf_read += (capacity()>>6) - num_buf_read;
        curr_idx = m_size / m_width;
    }
}

void custom_int_vector::reset()
{
    if(!ifs.is_open()) return;
    
    ifs.clear();
    ifs.seekg(9, std::ios::beg); // after buffer
    prev_idx = curr_idx = num_buf_read = 0;
}