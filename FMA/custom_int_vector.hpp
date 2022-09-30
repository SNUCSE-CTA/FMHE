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
/*! \file custom_int_vector.hpp
    \brief custom_int_vector.hpp contains a modified int_vector class of sdsl which uses a limited memory.
    \author Seunghwan Min
    \date 2017.1
*/

#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <sdsl/int_vector.hpp>

using namespace sdsl;

class custom_int_vector
{
    private:
        uint64_t        m_size;     // data total bit size
        uint8_t         m_width;    // width of the integer

        uint64_t        *prev_buf, *curr_buf;
        uint64_t        prev_idx, curr_idx;
        uint64_t        num_buf_read;   
        uint32_t        num_integers;   // number of integers in buffer
        uint32_t        buf_size;       
        
        std::string     fname;
        std::ifstream   ifs;

    
    public:
        custom_int_vector() {}
        custom_int_vector(char *file_name);
        ~custom_int_vector();
        
        void read_file(char *file_name);
        void read_file(std::string file_name);
        void read_header();
        void read_buffer();
        void reset();

        uint64_t capacity() const
        {
            return ((m_size+63)>>6)<<6;
        }

        uint64_t size() const
        {
            return m_size / m_width;
        }

        inline uint64_t operator[](const uint64_t& idx)
        {
            if(prev_idx !=0 && idx < prev_idx - num_integers) reset();
            while(idx >= curr_idx) read_buffer();
        
            uint64_t bits = (idx - prev_idx + num_integers) * m_width;
            if(idx >= prev_idx) bits -= num_integers * m_width;
        
            uint64_t *ptr;
            if(idx >= prev_idx) ptr = curr_buf;
            else ptr = prev_buf;
            ptr += bits >> 6;
        
            uint8_t offset = bits & 0x3F;
            uint8_t highbit = offset + m_width;
        
            uint64_t w1 = (*ptr) >> offset;
            if(highbit > 64)
            {
                return w1 | ((*(ptr+1) & bits::lo_set[highbit & 0x3F]) << (64 - offset));
            }
            else
            {
                return w1 & bits::lo_set[m_width];
            }
        } 
};