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
/*! \file saa_algorithm.hpp
    \brief saa_algorithm.hpp contains some helper classes for Suffix Array of Alignment (SAA), and pattern search(locate) and retrieval(extract) algorithms on SAAs.
    \author Hyunjoon Kim
    \reference 
    [1] J. Na, H. Kim, H. Park, T. Lecroq, M. Le'onard, L. Mouchar, and K. Park, FM-index of alignment: A compressed index for similar strings, TCS 638:159-170, 2016
    [2] J. Na, H. Kim, S. Min, H. Park, T. Lecroq, M. Le'onard, L. Mouchard, and K. Park, FM-index of alignment with gaps, TCS 710:148-157, 2018
    \date 2016.09
*/

#include <iterator>
#include <vector>
#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <utility>
#include <map>

#include <sdsl/suffix_array_helper.hpp>
#include <sdsl/iterators.hpp>
#include "construct_fma.hpp"
#include "gapinfo.hpp"

namespace sdsl
{
// psi[] trait
template<class t_csa,bool t_direction>
struct traverse_csaa_traits {
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::size_type size_type;
    static value_type access(const t_csa& csa,size_type i) {
        //char_type c = csa.F[i];
        //return csa.wavelet_tree.select(i - csa.C[csa.char2comp[c]] + 1 , c);
        return 0;
    }
    static value_type access(const t_csa& csa, size_type strn, size_type i) {
        return 0;
    }
};
// lf[] trait
template<class t_csa>
struct traverse_csaa_traits<t_csa,false> {
    typedef typename t_csa::value_type value_type;
    typedef typename t_csa::char_type char_type;
    typedef typename t_csa::size_type size_type;
    typedef typename t_csa::comp_char_type comp_char_type;
    static value_type access(const t_csa& csa, size_type strn, size_type i) {
        return csa.lf_value(strn, i);
    }
};

template<class t_csa,bool t_direction>
class traverse_csaa
{
    public:
        typedef typename t_csa::value_type                           value_type;
        typedef typename t_csa::size_type                             size_type;
        typedef typename t_csa::difference_type                 difference_type;
        typedef random_access_const_iterator<traverse_csaa>   const_iterator;

    private:
        const t_csa& m_csa;
    public:
        //! Constructor
        traverse_csaa(const t_csa& csaa) : m_csa(csaa) { }
        //! Copy constructor
        traverse_csaa(const traverse_csaa& tcsa) : m_csa(tcsa.m_csa) { }

        value_type value(size_type strn, size_type i) const 
        {
            assert(i < size());
            return traverse_csaa_traits<t_csa, t_direction>::access(m_csa, strn , i);
        }

        //! Returns the size of the \f$\LF\f$ function.
        size_type size() const {
            return m_csa.size();
        }

        //! Returns if the \f$\LF\f$ function is empty.
        size_type empty() const {
            return m_csa.empty();
        }

        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
            return const_iterator(this, size());
        }
};

template<class t_csa>
class isa_of_csaa
{
    public:
        typedef typename t_csa::value_type value_type;
        typedef typename t_csa::size_type size_type;
        typedef typename t_csa::difference_type difference_type;
        typedef random_access_const_iterator<isa_of_csaa> const_iterator;
    private:
        const t_csa& m_csa; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
        isa_of_csaa() {};    // disable default constructor
    public:
        //! Constructor
        isa_of_csaa(const t_csa& csaa) : m_csa(csaa) {}
        //! Calculate the Burrows Wheeler Transform (BWT) at position i.
        /*! \param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *  \par Time complexity
         *      \f$ \Order{\log |\Sigma|} \f$
         */
        value_type value(size_type strn, size_type i, GapInfo& gap_info)const {
            //GapInfo gap_info(m_csa.num_string, m_csa.gap_vector);
            size_type i0 = i;
            assert(i < m_csa.text_size);
            // get the leftmost sampled isa value to the right of i
            size_type ii = (i+m_csa.isa_sample_dens-1)/m_csa.isa_sample_dens;
            // std::cout << "isaa.value("<< strn<<","<<i<<"), leftmost align pos:"<<ii<<std::endl;
            value_type result = m_csa.index_table[ii];
            ii *= m_csa.isa_sample_dens;
            size_type ii_ind = gap_info.ConvertToInd(strn-1, ii);
            //bool leftToRight = false;
            if(result >= m_csa.size()){
                // std::cout << "result(non_shared): " << result << std::endl;
                result = m_csa.ns_isaa_sample[(result-m_csa.size())*m_csa.num_string+strn-1];
                 
                size_type ii = gap_info.ConvertToAlign(strn-1, ii_ind);
                
                //size_type ii_after = gap_info.ConvertToAlign(strn-1, ii_ind);
                //if(ii_after != ii) ii = ii_after;
                /*
                else{
                    ii = gap_info.ConvertToAlign(strn-1, ii_ind+1);
                    leftToRight = true;
                }
                */
                // endif
            }
            //value_type result = m_csa.isaa_sample[ ii*m_csa.num_string+2 ];
            // std::cout<<"isa_sample("<<strn<<","<<ii<<"): "<<result<<std::endl;
            if (ii < m_csa.text_size) {
                i = ii_ind - gap_info.ConvertToInd(strn-1, i);
                //if(leftToRight)   ++i;
            }  
            else {
                i = gap_info.ConvertToInd(strn-1, m_csa.text_size-1) - gap_info.ConvertToInd(strn-1, i);
            }
            size_type i1 = i;
            while (i--) {
                /*
                if(m_csa.has_multi_val[result] and strn != 0){
                    bit_vector::rank_1_type lf_ss_rank(&m_csa.has_multi_val);
                    result = m_csa.lf_ss[lf_ss_rank(result)*3+strn-1];
                }
                else{
                    result = m_csa.lf[result];
                }
                */
                result = m_csa.lf.value(strn, result); 
            }
            // std::cout << "After LF " << i1 << " times, isa.value(" << strn << ", " << i0 << "): " << result << std::endl;
            return result;
        }
        //! Returns the size of the BWT function.
        size_type size()const {
            return m_csa.size();
        }
        
        //! Returns if the BWT function is empty.
        size_type empty()const {
            return m_csa.empty();
        }
        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }
        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }
};

template<class t_csa, uint8_t int_width>
void set_isaa_samples(int_vector_buffer<int_width>& saa_buf, typename t_csa::isaa_sample_type& index_table, typename t_csa::isaa_sample_type& ns_isaa_sample, bit_vector& is_non_shared, typename t_csa::size_type& text_size, typename t_csa::size_type& num_string, typename t_csa::str_vector_type& strv, std::string gap_file)
{

    typedef typename t_csa::size_type size_type;
    auto saa_size = saa_buf.size()/2;
    bit_vector::rank_1_type non_shared_rank(&is_non_shared);
    size_type num_non_shared = non_shared_rank((size_type)is_non_shared.size());
    int_vector<> gap_vector;
    load_from_file(gap_vector, gap_file);
    GapInfo gap_info(num_string, gap_vector);
    
    index_table.width(bits::hi(saa_size+num_non_shared)+1);
    ns_isaa_sample.width(bits::hi(saa_size)+1);
    
    //std::cout<<"saa_size(set): "<<saa_size<<std::endl;
    if (saa_size >= 1) { // so n+t_csa::isa_sample_dens >= 2
        index_table.resize((text_size-1+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens + 1);
        //std::cout << "Suffix_array_helper: index_table size: " << index_table.size() << std::endl;
        ns_isaa_sample.resize((text_size/t_csa::isa_sample_dens)*num_string);
        // TODO: num_non_shared is too big and wrong. num_non_shared/isa_sample_dens+num_non_shared_region + 1
        //std::cout << "Suffix_array_helper: ns_isaa_sample size: " << text_size/t_csa::isa_sample_dens*num_string << std::endl;
        //std::cout << "isa_size: " << num_non_shared*num_string << std::endl;
    }
    util::set_to_value(index_table, 0);
    util::set_to_value(ns_isaa_sample, 0);
   
    // Mark the aligned positions converted from empty sampled positions
    size_type non_shared_cnt = 0;
    size_type isa_count = 0;
    std::map<size_type, size_type> aligned_to_saa_idx;
    size_type atos_cnt = 0;
    for(size_type i=0; i < index_table.size(); ++i){
        size_type pos = (i*t_csa::isa_sample_dens>text_size)?text_size-1:i*t_csa::isa_sample_dens;
        if(is_non_shared[pos]){
            index_table[i] = saa_size+non_shared_cnt;
            non_shared_cnt++;
            isa_count += num_string;
            for(size_type k = 0; k < num_string; ++k){
                size_type ind = gap_info.ConvertToInd(k, pos);
                size_type after = gap_info.ConvertToAlign(k, ind);
                if(pos != after){
                    if(!is_non_shared[after]){
                        ++atos_cnt;
                        aligned_to_saa_idx[after*(num_string+1)] = 1;
                    }
                    else{
                        ++atos_cnt;
                        aligned_to_saa_idx[after*(num_string+1)+k+1] = 1;
                    }
                }
            }
        }
    }
    //std::cout << "Suffix_array_helper: non_shared_cnt: " <<  non_shared_cnt << std::endl;
    //std::cout << "Suffix_array_helper: aligned_to_saa_idx size: " << aligned_to_saa_idx.size() << std::endl; 
    //std::cout << "Suffix_array_helper: aligned_to_saa_idx size: " << atos_cnt << std::endl; 

    // Search through a SAA for 1 and 2.
    // 1. Finding the corresponding SAA indexes to the marked aligned text positions
    // 2. Assigning values of ns_isaa_sample for non-shared positions and those of index_table for shared-positions
    atos_cnt = 0;
    for (size_type i=0; i < saa_buf.size()/2; ++i) {
        size_type strn  = saa_buf[2*i];
        size_type pos   = saa_buf[2*i+1];
        typename std::map<size_type, size_type>::iterator it;
        if(!is_non_shared[pos]){
            it = aligned_to_saa_idx.find(pos*(num_string+1));
            if(it!=aligned_to_saa_idx.end()){
                ++atos_cnt;
                aligned_to_saa_idx[pos*(num_string+1)] = i;
            }
        }
        else{
            for(size_type k = 0; k < num_string; ++k){
                if(strv[strn*num_string+k]){
                    it = aligned_to_saa_idx.find(pos*(num_string+1)+k+1);
                    if(it!=aligned_to_saa_idx.end()){
                        ++atos_cnt;
                        aligned_to_saa_idx[pos*(num_string+1)+k+1] = i;
                    }
                }
            }
        }
        size_type curr_idx;
        if ((pos % t_csa::isa_sample_dens) == 0)
            curr_idx = pos/t_csa::isa_sample_dens;
        else if (pos+1 == text_size)
            curr_idx = (pos+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens;
        else
            continue;
        
        if(is_non_shared[pos]){
            for(size_type k = 0; k < num_string; ++k){
                if(strv[strn*num_string+k]){
                    ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k] = i;
                }
            }
        }
        else{
            index_table[curr_idx] = i;
        }
    }
    //std::cout << "Suffix_array_helper: aligned_to_saa_idx size: " << aligned_to_saa_idx.size() << std::endl; 
    //std::cout << "Suffix_array_helper: aligned_to_saa_idx size: " << atos_cnt << std::endl; 

    // Assign values of ns_isaa_sample for not-yet-assigned non-shared positions (or empty aligned text positions).
    atos_cnt = 0;
    for(size_type i=0; i < saa_buf.size()/2; ++i){
        size_type strn  = saa_buf[2*i];
        size_type pos   = saa_buf[2*i+1];
        size_type curr_idx;
        if ((pos % t_csa::isa_sample_dens) == 0)
            curr_idx = pos/t_csa::isa_sample_dens;
        else if (pos+1 == text_size)
            curr_idx = (pos+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens;
        else
            continue;
        if(is_non_shared[pos]){
            for(size_type k = 0; k < num_string; ++k){
                if(!ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k]){
                    size_type ind = gap_info.ConvertToInd(k, pos);
                    size_type after = gap_info.ConvertToAlign(k, ind);
                    ++atos_cnt;
                    size_type answer = aligned_to_saa_idx[after*(num_string+1)+k+1]?aligned_to_saa_idx[after*(num_string+1)+k+1]:aligned_to_saa_idx[after*(num_string+1)];
                    ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k] = answer;
                }
            }
        }
    }
    ns_isaa_sample.resize(isa_count);
    //std::cout << "Suffix_array_helper: ns_isaa_sample size(answer): "<< ns_isaa_sample.size() << std::endl;
}


template<class t_csa>
void build_isaa_samples(typename t_csa::saa_sample_type& saa_sample, typename t_csa::isaa_sample_type& index_table, typename t_csa::isaa_sample_type& ns_isaa_sample, typename t_csa::size_type& text_size, typename t_csa::size_type& num_string, typename t_csa::gap_type& gap_vector)
{
    typedef typename t_csa::size_type size_type;
    GapInfo gap_info(num_string, gap_vector);

    auto saa_size = saa_sample.marked.size();
    size_type index_table_size = (text_size-1+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens + 1;
    size_type num_non_shared = 0; 
    // For aligned positions pos where pos % isa_sample_dens == 0 or pos == text_size - 1
    // is_non_shared[pos/isa_sample_dens] == 1 if (aligned) pos is a non-shared position
    // is_non_shared[pos/isa_sample_dens] == 0 otherwise
    bit_vector is_non_shared(index_table_size, 1);
    for(size_type i = 0; i < saa_sample.size(); ++i){
        size_type strn = saa_sample[i]/text_size;
        size_type pos  = saa_sample[i]%text_size; 
        if (pos % t_csa::isa_sample_dens == 0 || pos+1 == text_size){
            size_type curr_idx = (pos%t_csa::isa_sample_dens==0)?(pos/t_csa::isa_sample_dens):((pos+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens);
            if(strn == 0)
                is_non_shared[curr_idx] = 0;
            else
                ++num_non_shared;
        }
    }

    index_table.width(bits::hi(saa_size+num_non_shared)+1);
    ns_isaa_sample.width(bits::hi(saa_size)+1);
   
    if (saa_size >= 1) { // so n+t_csa::isa_sample_dens >= 2
        index_table.resize(index_table_size);
        //std::cout << "Suffix_array_helper: index_table size: " << index_table.size() << std::endl;
        ns_isaa_sample.resize((text_size/t_csa::isa_sample_dens)*num_string);
        //ns_isaa_sample.resize(num_non_shared*num_string);
        // TODO: num_non_shared is too big and wrong. num_non_shared/isa_sample_dens+num_non_shared_region + 1
        //std::cout << "Suffix_array_helper: ns_isaa_sample size: " << text_size/t_csa::isa_sample_dens*num_string << std::endl;
    }
    util::set_to_value(index_table, 0);
    util::set_to_value(ns_isaa_sample, 0);
    

    // For non-shared positions which are multiples of isa_sample_dens (pos), check if pos is a gap for each sequence
    // If pos is a gap, find the leftmost non-gap position on the right (called after).
    // For each sequence k, 
    // aligned_to_saa_idx[after*num_string+k] == 1 if pos is a gap
    // aligned_to_saa_idx[after*num_string+k] == 0 otherwise 
    
    // std::cout<<"For non-shared positions which are multiples of isa_sample_dens (pos), check if pos is a gap for each sequence"<<std::endl;
    size_type non_shared_cnt = 0;
    size_type isa_count = 0;
    std::map<size_type, size_type> aligned_to_saa_idx;
    size_type atos_cnt = 0;
    for(size_type i=0; i < index_table.size(); ++i){
        size_type pos = (i*t_csa::isa_sample_dens>text_size)?text_size-1:i*t_csa::isa_sample_dens;
        if(is_non_shared[i]){ // if pos is in non-shared region
            index_table[i] = saa_size+non_shared_cnt;
            ++non_shared_cnt;
            isa_count += num_string;
            for(size_type k = 0; k < num_string; ++k){
                size_type ind   = gap_info.ConvertToInd(k, pos);
                size_type after = gap_info.ConvertToAlign(k, ind);
                if(pos != after){
                    aligned_to_saa_idx[after*num_string+k] = 1;
                    ++atos_cnt;
                }
            }
        }
    }
    
    ns_isaa_sample.resize(isa_count);
    // std::cout << "Suffix_array_helper: ns_isaa_sample size(answer): "<< ns_isaa_sample.size() << std::endl;

    // Assign the corresponding SAA values to the (sequence, position) pairs in the aligned_to_saa_idx
    atos_cnt = 0; 
    // std::cout<<"Assign the corresponding SAA values to the (sequence, position) pairs in the aligned_to_saa_idx"<<std::endl;
    bit_vector::select_1_type select_saa_sample(&saa_sample.marked);
    for(size_type i = 0; i < saa_sample.size(); ++i){
        size_type strn = saa_sample[i]/text_size;
        size_type pos  = saa_sample[i]%text_size;
        if(strn != 0){ 
            typename std::map<size_type, size_type>::iterator it;
            for(size_type k = 0; k < num_string; ++k){
                if(saa_sample.strv[strn*num_string+k]){
                    it = aligned_to_saa_idx.find(pos*num_string+k);
                    if(it!=aligned_to_saa_idx.end()){
                        size_type ss = select_saa_sample(i+1);
                        aligned_to_saa_idx[pos*num_string+k] = ss;
                        ++atos_cnt;
                    }
                }
            }
        }
    }

    // Make index_table and ns_isaa_sample
    // std::cout<<"Make index_table and ns_isaa_sample"<<std::endl;
    for(size_type i = 0; i < saa_sample.size(); ++i){
        size_type strn = saa_sample[i]/text_size;
        size_type pos  = saa_sample[i]%text_size;
        if ((pos % t_csa::isa_sample_dens) == 0 or pos+1 == text_size){
            size_type curr_idx = (pos%t_csa::isa_sample_dens==0)?(pos/t_csa::isa_sample_dens):((pos+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens);
            size_type saa_idx = select_saa_sample(i+1);
            if(strn == 0){  // a shared position
                index_table[curr_idx] = saa_idx;
            }
            else{           // a non-shared position
                for(size_type k = 0; k < num_string; ++k){
                    if(saa_sample.strv[strn*num_string+k]){
                        ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k] = saa_idx;
                    }
                }
            }
        }
    }
    for(size_type i=0; i < saa_sample.size(); ++i){
        size_type strn = saa_sample[i]/text_size;
        size_type pos  = saa_sample[i]%text_size;
        size_type curr_idx;
        if ((pos % t_csa::isa_sample_dens) == 0 or pos+1 == text_size){
            size_type curr_idx = (pos%t_csa::isa_sample_dens==0)?(pos/t_csa::isa_sample_dens):((pos+t_csa::isa_sample_dens-1)/t_csa::isa_sample_dens);
            if(strn != 0){
                for(size_type k = 0; k < num_string; ++k){
                    if(!ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k]){
                        size_type ind = gap_info.ConvertToInd(k, pos);
                        size_type after = gap_info.ConvertToAlign(k, ind);
                        ns_isaa_sample[(index_table[curr_idx]-saa_size)*num_string+k] = aligned_to_saa_idx[after*num_string+k];
                    }
                }
            }
        }
    }
}

//! Backward search for a character c in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa CSA type.
 *
 * \param csa    The CSA object.
 * \param l      Left border of the interval \f$ [\ell..r]\f$.
 * \param r      Right border of the interval \f$ [\ell..r]\f$.
 * \param c      Character to be prepended to \f$\omega\f$.
 * \param l_res  New left border.
 * \param r_res  Right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *         \f$ \Order{ t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */

template<class t_csa>
bool intersection(bit_vector& Z, bit_vector& set, typename t_csa::size_type num_string)
{
    typename t_csa::size_type iteration = num_string/64;
    typename t_csa::size_type remainder = num_string%64;
    bool empty = true;
    uint64_t res;
    typename t_csa::size_type j;
    for(j = 0; j < iteration; ++j){
        res = set.get_int(j*64) & Z.get_int(j*64);
        Z.set_int(j*64, res);
        if(res != 0) empty = false;
    }
    if(remainder != 0){
        res = set.get_int(j*64, remainder) & Z.get_int(j*64, remainder);
        Z.set_int(j*64, res, remainder);
        if(res != 0) empty = false;
    }
    return empty;
}

//bool update_set(bit_vector& Zm, const t_csa& csa, typename t_csa::size_type j, typename t_csa::size_type size, typename t_csa::size_type& mapping, double& time)
// Return a flag to break or not
template<class t_csa>
bool update_set(bit_vector& Zm, const t_csa& csa, typename t_csa::size_type j, typename t_csa::size_type size)
{
    typename t_csa::pair_type p = csa[j];
    typename t_csa::size_type s = p.first;
        
    if(s > size){
        typename t_csa::size_type iteration = Zm.size()/64;
        typename t_csa::size_type remainder = Zm.size()%64;
        typename t_csa::size_type x;
        uint64_t res;
        typename t_csa::str_vector_type vec = csa.saa_sample.strv;
        for(x = 0; x < iteration; ++x){
            res = Zm.get_int(x*64) | vec.get_int(s*size+x*64);
            Zm.set_int(x*64, res);
        }
        if(remainder != 0){
            res = Zm.get_int(x*64, remainder) | vec.get_int(s*size+x*64, remainder);
            Zm.set_int(x*64, res, remainder);
        }
        return false;
    }
    else if(s != 0){
        Zm[s-1] = 1;
        return false;
    }
    else{
        for(typename t_csa::size_type x = 0; x < size; ++x)
            Zm[x] = 1;
        return true;
    }
}

template<class t_csa, class t_pat_iter>
typename t_csa::size_type backward_search(
    const t_csa& csa,
    t_pat_iter begin,
    t_pat_iter it,
    typename t_csa::size_type& first,
    typename t_csa::size_type& last,
    bit_vector& Z
)
{
    typename t_csa::size_type l, r, c_begin;
    typename t_csa::char_type c;
    bool isEmpty = false;
    while (begin < it) {
        --it;
        c = *it;
        l = first; r = last;
        // std::cout << "c: " << c << ", l: " << l << ", r: " << r << ", csa.size(): " << csa.size() << ", c_begin: " << csa.C[csa.char2comp[c]] << std::endl;
        assert(l <= r); assert(r < csa.size());
        if(l == 0  and r == 0) return 0;
        c_begin = csa.C[csa.char2comp[c]];
        first = c_begin + csa.bwt.rank(l, c); // count c in bwt[0..l-1] + 1
        last  = c_begin + csa.bwt.rank(r+1, c) - 1; // count c in bwt[0..r]
        if(first >= last)
        {      
            typename t_csa::size_type vec_size = Z.size();
            bit_vector Zm = bit_vector(vec_size, 0);
            if(l < r){
                if(csa.get_extra(l, c)){
                    update_set(Zm, csa, l, vec_size);
                    for(typename t_csa::size_type h = l+1; h <= r; ++h){
                        if(!csa.get_extra(h, c)) break;
                        if(update_set(Zm, csa, h, vec_size)) break;
                    }
                    isEmpty = intersection<t_csa>(Z, Zm, vec_size);
                    if(!isEmpty) // if Z is not empty
                        first = last;
                    else{ // if Z is empty
                        first = last + 1;
                        break;
                    }
                }
                else if(csa.get_extra(r, c)){
                    update_set(Zm, csa, r, vec_size);
                    for(typename t_csa::size_type h = r-1; h >= l; --h){
                        if(!csa.get_extra(h, c)) break;
                        if(update_set(Zm, csa, h, vec_size)) break;
                    }
                    isEmpty = intersection<t_csa>(Z, Zm, vec_size);
                    if(!isEmpty) // if Z is not empty
                        first = last;
                    else{ // if Z is empty
                        first = last + 1;
                        break;
                    }
                }
            }    
            else if (l == r){ 
                if(csa.get_extra(l, c)){
                    update_set(Zm, csa, l, vec_size);
                    isEmpty = intersection<t_csa>(Z, Zm, vec_size);
                    if(!isEmpty) // if Z is not empty
                        first = last;
                    else{ // if Z is empty
                        first = last + 1;
                        isEmpty = true;
                    }
                }
                break;
            }
        }
    }
    while(begin < it and !isEmpty){
        --it;
        c = *it;
        l = first; r = last;
        if(l == 0 and r == 0) return 0;
        c_begin = csa.C[csa.char2comp[c]];
        first = c_begin + csa.bwt.rank(l, c); // count c in bwt[0..l-1] + 1
        last  = c_begin + csa.bwt.rank(r+1, c) - 1; // count c in bwt[0..r]
        if(l <= r and csa.get_extra(l, c)){
            typename t_csa::size_type vec_size = Z.size();
            bit_vector Zm = bit_vector(vec_size, 0);
            //update_set(Zm, csa, l, vec_size, mapping, time);
            update_set(Zm, csa, l, vec_size);
            isEmpty = intersection<t_csa>(Z, Zm, vec_size);
            
            if(!isEmpty) // if Z is not empty
                first = last;
            else{ // if Z is empty
                first = last + 1;
                isEmpty = true;
                break;
            }
            
        }
    }

    //std::cout << "Backward_search: [" << c << "] First: " << first << ", Last: " << last << std::endl;
    assert(last-first+1 >= 0);
    //std::cout << "Backward_search return val = " << last+1-first << std::endl;
    return last-first+1;
}


//! Backward search for a pattern in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa      A CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa   The CSA object.
 * \param l     Left border of the lcp-interval \f$ [\ell..r]\f$.
 * \param r     Right border of the lcp-interval \f$ [\ell..r]\f$.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \param l_res New left border.
 * \param r_res New right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *       \f$ \Order{ len \cdot t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */

// Modified
template<class t_csa, class t_pat_iter>
typename t_csa::size_type
backward_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    bit_vector& Z,
    SDSL_UNUSED typename std::enable_if<std::is_same<csaa_tag, typename t_csa::index_category>::value, csaa_tag>::type x = csaa_tag()
)
{
    typename t_csa::size_type num_string = csa.num_string;
    t_pat_iter it = end;
    typename t_csa::comp_char_type cc = csa.char2comp[*(--it)];
    l = csa.C[cc];
    r = csa.C[cc+1] - 1;
    // std::cout << "Backward_search: start with " << *it << "[" << cc << "]: " << l << " ~ " << r << std::endl;
    backward_search(csa, begin, it, l, r, Z);
    if(l == 0 and r == 0) return 0;
    typename t_csa::size_type cnt = 0;
    
    typename t_csa::size_type iteration = num_string/64;
    typename t_csa::size_type remainder = num_string%64;
    
    for(typename t_csa::size_type i = l; i <= r; ++i){
        typename t_csa::size_type s = csa[i].first;
        //typename t_csa::size_type s = csa.access(i, mapping, time).first;
        if(l != r){
            uint64_t copy;
            typename t_csa::size_type j, idx = s*num_string;
            for(j = 0; j < iteration; ++j){
                copy = csa.saa_sample.strv.get_int(j*64+idx);
                cnt += __builtin_popcountl(copy);
            }
            if(remainder != 0){
                copy = csa.saa_sample.strv.get_int(j*64+idx, remainder);
                cnt += __builtin_popcountl(copy);
            }
        }
        else{
            uint64_t copy, zCopy;
            typename t_csa::size_type j, idx = s*num_string;
            for(j = 0; j < iteration; ++j){
                copy = csa.saa_sample.strv.get_int(j*64+idx);
                zCopy = Z.get_int(j*64);
                cnt += __builtin_popcountl(copy&zCopy);
            }
            if(remainder != 0){
                copy = csa.saa_sample.strv.get_int(j*64+idx, remainder);
                zCopy = Z.get_int(j*64, remainder);
                cnt += __builtin_popcountl(copy&zCopy);
            }
        }
    }

    l_res = l;
    r_res = r;
    //std::cout<<"Bachward_search: cnt = "<<cnt<<std::endl;
    return cnt;
}

template<class t_csa, class t_pat_iter, class t_rac=int_vector<64>>
t_rac pattern_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    bit_vector& Z 
)
{
    typename t_csa::size_type num_string = csa.num_string;
    t_pat_iter it = end;
    typename t_csa::comp_char_type cc = csa.char2comp[*(--it)];
    l = csa.C[cc];
    r = csa.C[cc+1] - 1;
    //std::cout << "Backward_search: start with " << *it << "[" << cc << "]: " << l << " ~ " << r << std::endl;
    backward_search(csa, begin, it, l, r, Z);
    
    t_rac occ((r-l+1)*2*num_string);
    typename t_csa::size_type cnt = 0;
    
    if(l < r){
        for(typename t_csa::size_type i = l; i <= r; ++i){
            //typename t_csa::pair_type p = csa.access(i, mapping, time);
            typename t_csa::pair_type p = csa[i];
            typename t_csa::size_type s = p.first;
            typename t_csa::size_type n = p.second;
            //std::cout <<"["<<i<<"]: ("<<s<<", "<<n<<")"<<std::endl;
            for(typename t_csa::size_type j = 0; j < num_string; ++j){
                if(csa.saa_sample.strv[s*num_string+j]){
                    occ[2*cnt]      = j+1;
                    occ[2*cnt+1]    = n;
                    ++cnt; 
                }
            }
        }
    }
    else if (l == r){
        //typename t_csa::pair_type p = csa.access(l, mapping, time);
        typename t_csa::pair_type p = csa[l];
        typename t_csa::size_type s = p.first;
        typename t_csa::size_type n = p.second;
        //std::cout <<"["<<l<<"]: ("<<s<<", "<<n<<")"<<std::endl;
        for(typename t_csa::size_type j = 0; j < num_string; ++j){
            if(csa.saa_sample.strv[s*num_string+j] and Z[j]){
                occ[2*cnt]      = j+1;
                occ[2*cnt+1]    = n;
                ++cnt; 
            }
        }
    }
    occ.resize(cnt*2);
    l_res = l;
    r_res = r;
    
    return occ;

}

//! Counts the number of occurrences of a pattern in a CSA.
/*!
 * \tparam t_csa      CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa   The CSA object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return The number of occurrences of the pattern in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} } \f$
 */

// Modified
template<class t_csa, class t_pat_iter>
typename t_csa::size_type count(
    const t_csa& csa,
    t_pat_iter begin,
    t_pat_iter end,
    csaa_tag
)
{
    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;
    if (end == begin) 
        return 0;
    typename t_csa::size_type occ_begin, occ_end; // dummy variable for the backward_search call
    bit_vector Z = bit_vector(csa.num_string, 1);
    typename t_csa::size_type result = backward_search(csa, 0, csa.size()-1, begin, end, occ_begin, occ_end, Z);
    return result;
}

/*
template<class t_csx, class t_pat_iter>
typename t_csx::size_type count(
    const t_csx& csx,
    t_pat_iter begin,
    t_pat_iter end
)
{
    typename t_csx::index_category tag;
    return count(csx, begin, end, tag);
}
*/

//! Counts the number of occurrences of a pattern in a CSA.
/*!
 * \tparam t_csa      CSA type.
 *
 * \param csa The CSA object.
 * \param pat The pattern.
 * \return The number of occurrences of the pattern in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} } \f$
 */

/*
template<class t_csx>
typename t_csx::size_type count(
    const t_csx& csx,
    const typename t_csx::string_type& pat
)
{
    typename t_csx::index_category tag;
    return count(csx, pat.begin(), pat.end(), tag);
}
*/

//! Calculates all occurrences of a pattern pat in a CSA.
/*!
 * \tparam t_csa      CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa   The CSA object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */

// Modified
template<class t_csa, class t_pat_iter, class t_rac=int_vector<64>>
t_rac locate(
    const t_csa&  csa,
    t_pat_iter begin,
    t_pat_iter end,
    SDSL_UNUSED typename std::enable_if<std::is_same<csaa_tag, typename t_csa::index_category>::value, csaa_tag>::type x = csaa_tag()
)
{
    typedef typename t_csa::size_type size_type;
    size_type num_string = csa.num_string;
    size_type occ_begin, occ_end, occs;
    typename t_csa::comp_char_type cc = csa.char2comp[*end];
    size_type first = csa.C[cc]-1;
    size_type last  = csa.C[cc+1]-1;
    //std::cout << *end << ", " << cc << ": " << first << " ~ " << last << std::endl;
    bit_vector Z = bit_vector(csa.num_string, 1);
    t_rac occ = pattern_search(csa, 0, csa.size()-1, begin, end, occ_begin, occ_end, Z);
    return occ;
}


//! Calculates all occurrences of a pattern pat in a CSA/CST.
/*!
 * \tparam t_csa      CSA/CST type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa  The CSA/CST object.
 * \param pat  The pattern.
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
/*
template<class t_csx, class t_rac=int_vector<64>>
t_rac locate(
    const t_csx&  csx,
    const typename t_csx::string_type& pat
)
{
    typename t_csx::index_category tag;
    return locate<t_csx, decltype(pat.begin()), t_rac>(csx, pat.begin(), pat.end(), tag);
}
*/


//! Writes the substring T[begin..end] of the original text T to text[0..end-begin+1].
/*!
 * \tparam t_csa       CSA type.
 * \tparam t_text_iter Random access iterator type.
 *
 * \param csa   The CSA object.
 * \param begin Position of the first character which should be extracted (inclusive).
 * \param end   Position of the last character which should be extracted (inclusive).
 * \param text  Random access iterator pointing to the start of an container, which can hold at least (end-begin+1) character.
 * \returns The length of the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *        \f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */
/*
template<class t_csa, class t_text_iter>
typename t_csa::size_type extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    t_text_iter text,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    typename t_csa::extract_category extract_tag;
    return extract(csa, begin, end, text, extract_tag);
}
*/

//! Specialization of extract for LF-function based CSAs

template<class t_csa, class t_text_iter>
typename t_csa::size_type extract(
    const t_csa& csa,
    typename t_csa::size_type strn,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    typename t_csa::size_type steps,
    t_text_iter text,
    GapInfo& gap_info
)
{
    assert(end < csa.text_size);
    assert(begin <= end);
    if(end >= csa.text_size)   
        end = csa.text_size-1;
    typename t_csa::size_type value = csa.isa.value(strn, end, gap_info);
    while (steps != 0) {
        text[steps-1] = first_row_symbol(value, csa);
        value = csa.lf_value(strn, value);
        --steps;
    }
    return end-begin+1;
}

//! Reconstructs the substring T[begin..end] of the original text T to text[0..end-begin+1].
/*!
 * \tparam t_rac Random access container which should hold the result.
 * \tparam t_csa CSA type.
 *
 * \param csa   The CSA object.
 * \param begin Position of the first character which should be extracted (inclusive).
 * \param end   Position of the last character which should be extracted (inclusive).
 * \return A t_rac object holding the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *        \f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */

template<class t_csa>
typename t_csa::string_type extract(
    const t_csa& csa,
    typename t_csa::size_type strn,
    typename t_csa::size_type begin,
    typename t_csa::size_type end, 
    typename t_csa::size_type steps,
    GapInfo& gap_info,
    SDSL_UNUSED typename std::enable_if<std::is_same<csaa_tag, typename t_csa::index_category>::value, csaa_tag>::type x = csaa_tag()
)
{
    //std::cout<<"extract("<<strn<<", "<<begin<<" ~ "<<end<<")"<<std::endl;
    assert(end <= csa.text_size);
    assert(begin <= end);
    typedef typename t_csa::string_type string_type;
    string_type result(end-begin+1, (typename string_type::value_type)0);
    extract(csa, strn, begin, end, steps, result.begin(), gap_info);
    return result;
}

} // end namespace
