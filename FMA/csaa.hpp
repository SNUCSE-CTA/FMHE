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
/*! \file csaa.hpp
    \brief csaa.hpp contains an implementation of FM-index of alignment with gaps.
    \author Hyunjoon Kim
    \reference 
    [1] J. Na, H. Kim, H. Park, T. Lecroq, M. Le'onard, L. Mouchar, and K. Park, FM-index of alignment: A compressed index for similar strings, TCS 638:159-170, 2016
    [2] J. Na, H. Kim, S. Min, H. Park, T. Lecroq, M. Le'onard, L. Mouchard, and K. Park, FM-index of alignment with gaps, TCS 710:148-157, 2018
    \date 2016.09
*/

#include <cstdlib>  // for exit(). should be removed after debugging.
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>
#include <utility>

#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/suffix_array_helper.hpp>
#include <sdsl/util.hpp>
#include <sdsl/csa_sampling_strategy.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>

#include "saa_algorithm.hpp"
#include "gapinfo.hpp"

// 160605 add temporarily
//#include <map>

// Add ctime
#include <ctime>

namespace sdsl
{

//! A class for the Compressed Suffix Array (CSA) proposed by Sadakane for practical implementation.
/*!
  *  \tparam t_enc_vec         Space-efficient vector for increasing integer sequences.
  *  \tparam t_dens            Sampling density of SA values
  *  \tparam t_int_dens        Sampling density of ISA values
  *  \tparam t_sa_sample_strat Policy of SA sampling. E.g. sample in SA-order or text-order.
  *  \tparam t_isa             Vector type for ISA sample values.
  *  \tparam t_alphabet_strat  Policy for alphabet representation.
  *
  *  \sa sdsl::csa_wt, sdsl::csa_bitcompressed
  * @ingroup csa
 */
template<class t_bitvector       = rrr_vector<127>,         // Vector type used to store the counters and extras
         uint32_t t_dens         = 32,                    // Sample density for suffix array (SA) values
         uint32_t t_inv_dens     = 64,                    // Sample density for inverse suffix array (ISA) values
         class t_sa_sample_strat = saa_sampling<>,// Policy class for the SA sampling. Alternative text_order_sa_sampling.
         class t_isa             = int_vector<>,          // Container for the ISA samples.
         class t_alphabet_strat  = saa_alphabet<>          // Policy class for the representation of the alphabet.
         >
class csaa
{
        friend class bwt_of_csa_psi<csaa>;
    public:
        enum { sa_sample_dens = t_dens,
               isa_sample_dens = t_inv_dens
             };

        typedef uint64_t                                                        value_type;
        typedef random_access_const_iterator<csaa>                              const_iterator;
        typedef const_iterator                                                  iterator;
        typedef const value_type                                                const_reference;
        typedef const_reference                                                 reference;
        typedef const_reference*                                                pointer;
        typedef const pointer                                                   const_pointer;
        typedef int_vector<>::size_type                                         size_type;
        typedef size_type                                                       csa_size_type;
        typedef ptrdiff_t                                                       difference_type;
        typedef enc_vector<>                                                    enc_vector_type;
        //typedef enc_vector_type                                                 lf_type;
        //typedef enc_vector_type                                                 lf_ss_type;
        typedef traverse_csaa<csaa,false>                                       lf_type;
        typedef bwt_of_csa_wt<csaa>                                             bwt_type;
        typedef isa_of_csaa<csaa>                                               isa_type;
        typedef text_of_csa<csaa>                                               text_type;
        typedef first_row_of_csa<csaa>                                          first_row_type;
        typedef typename t_sa_sample_strat::template type<csaa>::sample_type    saa_sample_type;
        typedef t_isa                                                           isaa_sample_type;
        typedef int_vector<>                                                    gap_type;
        typedef t_alphabet_strat                                                alphabet_type;
        typedef typename alphabet_type::alphabet_category                       alphabet_category;
        typedef typename alphabet_type::comp_char_type                          comp_char_type;
        typedef typename alphabet_type::char_type                               char_type; // Note: This is the char type of the CSA not the WT!
        typedef typename alphabet_type::string_type                             string_type;
        typedef csaa 	                                                        csa_type;
        typedef std::pair<size_type, size_type>                                 pair_type;
        typedef rrr_vector<127>                                                 str_vector_type;
        typedef t_bitvector                                                     bit_vector_type;
        typedef typename t_bitvector::rank_1_type                               rank_1_type;
        typedef typename t_bitvector::select_1_type                             select_1_type;
        typedef typename t_bitvector::select_0_type                             select_0_type;
        typedef csaa_tag                                                        index_category;
        typedef lf_csaa_tag                                                     extract_category;

        //friend class traverse_csaa<csaa,true>;
        friend class traverse_csaa<csaa,false>;

        static const uint32_t linear_decode_limit = 100000;
    private:
        gap_type              m_gap_vector;
        //lf_type               m_lf;         // lf function
        //lf_ss_type              m_lf_ss;      // lf_ss function
        saa_sample_type         m_saa_sample;  // suffix array samples
        isaa_sample_type        m_index_table; // indexes of inverse suffix array samples
        isaa_sample_type        m_ns_isaa_sample; // inverse suffix array samples of non-shared region
        
        // 160603 add
        //isaa_sample_type        m_index_table2; // indexes of inverse suffix array samples
        //isaa_sample_type        m_ns_isaa_sample2; // inverse suffix array samples of non-shared region
        
        alphabet_type           m_alphabet;   // alphabet component
		size_type 		        m_num_string;
        //size_type               m_lf_size;
        //bit_vector              m_shared_start, m_non_shared_start;
        //bit_vector              m_has_multi_val;
        bit_vector_type         m_bv[conf::CATEGORY_NUM][conf::CHAR_NUM];
        rank_1_type             m_rank[conf::CATEGORY_NUM][conf::CHAR_NUM];
        size_type               m_text_size; 
        size_type               m_saa_size;
        std::string             m_info_dir;
        //mutable std::vector<uint64_t> m_lf_buf; // buffer for decoded lf values

        void copy(const csaa& csa) {
            m_gap_vector        = csa.m_gap_vector;
            //m_lf                = csa.m_lf;
            //m_lf_ss             = csa.m_lf_ss;
            m_saa_sample        = csa.m_saa_sample;
            m_index_table       = csa.m_index_table;
            m_ns_isaa_sample    = csa.m_ns_isaa_sample;
            
            // 160603
            //m_index_table2       = csa.m_index_table2;
            //m_ns_isaa_sample2    = csa.m_ns_isaa_sample2;
            
            m_alphabet          = csa.m_alphabet;
			m_num_string        = csa.m_num_string;
            m_text_size         = csa.m_text_size;
            m_saa_size          = csa.m_saa_size;
        };
        /*
        void create_buffer() {
            if (enc_vector_type::sample_dens < linear_decode_limit) {
                m_lf_buf = std::vector<uint64_t>(enc_vector_type::sample_dens+1);
            }
        }
        */
    public:
        const typename alphabet_type::char2comp_type& char2comp         = m_alphabet.char2comp;
        const typename alphabet_type::comp2char_type& comp2char         = m_alphabet.comp2char;
        const typename alphabet_type::C_type&         C                 = m_alphabet.C;
        const typename alphabet_type::sigma_type&     sigma             = m_alphabet.sigma;
        //const lf_type&                                lf                = m_lf;
        //const lf_ss_type&                             lf_ss             = m_lf_ss;
        const lf_type                                 lf                = lf_type(*this);
        const bwt_type                                bwt               = bwt_type(*this);
        const isa_type                                isa               = isa_type(*this);
        const bwt_type                                L                 = bwt_type(*this);
        const first_row_type                          F                 = first_row_type(*this);
        const text_type                               text              = text_type(*this);
        const saa_sample_type&                        saa_sample        = m_saa_sample;
        const isaa_sample_type&                       index_table       = m_index_table;
        const isaa_sample_type&                       ns_isaa_sample    = m_ns_isaa_sample;
        const gap_type&                             gap_vector      = m_gap_vector;
        
        // 160603
        //const isaa_sample_type&                       index_table2       = m_index_table2;
        //const isaa_sample_type&                       ns_isaa_sample2    = m_ns_isaa_sample2;
		
        const size_type& 							  num_string        = m_num_string;
        //const size_type&                              lf_size           = m_lf_size;
        //const bit_vector&                             shared_start      = m_shared_start;
        //const bit_vector&                             non_shared_start  = m_non_shared_start;
        //const bit_vector&                             has_multi_val     = m_has_multi_val;
        const size_type&                              text_size         = m_text_size; 
        const size_type&                              saa_size          = m_saa_size; 
        const std::string&                            info_dir          = m_info_dir;
        //const (bit_vector_type&)                             bv[conf::CATEGORY_NUM][conf::CHAR_NUM];
        //const (rank_1_type&)                rank[conf::CATEGORY_NUM][conf::CHAR_NUM];
        /*
        for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
            for(size_type j = 0; j < conf::CHAR_NUM; ++j){
                bv[i][j] = m_bv[i][j];
                rank[i][j] = m_rank[i][j];
            }
        }
        */
        //! Default Constructor
        /*
        csaa() {
            //create_buffer();
        }
        */
        //csaa(int number, std::string str) {
        csaa(int number, std::string str, int length) {
            m_info_dir      = str;
            for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
                for(size_type j = 0; j < conf::CHAR_NUM; ++j){
                    m_rank[i][j] = rank_1_type();
                }
            }
            m_num_string = number;
            m_text_size = length;
        }
        //! Default Destructor
        ~csaa() { }
        
        //! Copy constructor
        csaa(const csaa& csa) {
            //create_buffer();
            copy(csa);
        }

        //! Move constructor
        csaa(csaa&& csa) {
            *this = std::move(csa);
        }

        csaa(
            cache_config& ref_config, 
            str_vector_type& strv,
            bit_vector& is_non_shared, 
            bit_vector& sampled, 
            bit_vector_type vector[][conf::CHAR_NUM],
            gap_type& gap_vector,
            const size_type& number
        );

        void set_env(size_type number, std::string str, size_type text_len){
            m_info_dir = str;
            m_num_string = number;
            m_text_size = text_len;
        }
        
        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const {
            //return m_has_multi_val.size();
            return m_saa_size; 
        }

        //! Returns the largest size that csaa can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return enc_vector_type::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            //return m_has_multi_val.empty();
            return (m_saa_size == 0);
        }

        //! Swap method for csaa
        /*! The swap method can be defined in terms of assignment.
            This requires three assignments, each of which, for a container type, is linear
            in the container's size. In a sense, then, a.swap(b) is redundant.
            This implementation guaranties a run-time complexity that is constant rather than linear.
            \param csa csaa to swap.

            Required for the Assignable Conecpt of the STL.
          */
        void swap(csaa& csa);

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
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

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Required for the STL Random Access Container Concept.
         * \par Time complexity
         *      \f$ \Order{s_{SA}\cdot t_{\Psi}} \f$, where every \f$s_{SA}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline pair_type operator[](size_type i)const;
        
        //inline pair_type access(size_type i, size_type& mapping, double& time)const;

        //! Assignment Copy Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        csaa& operator=(const csaa& csa) {
            if (this != &csa) {
                copy(csa);
            }
            return *this;
        }

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load from a stream.
        /*! \param in Input stream to load the data structure from.
         */
        void load(std::istream& in);

        uint32_t get_sample_dens() const {
            return t_dens;
        }
        
        size_type char_index(const char_type c) const {
            if(c == conf::CHAR[0])
                return 0;
            else if(c == conf::CHAR[1])
                return 1;
            else if(c == conf::CHAR[2])
                return 2;
            else if(c == conf::CHAR[3])
                return 3;
            else if(c == conf::CHAR[4])
                return 4;
            else{
                std::cout<<"csaa.hpp: char_index error ["<<(char)c<<"]"<<std::endl;
                exit(1);
            }
        }
        
        size_type get_extra(size_type i, const char_type c) const {
            return m_bv[1][char_index(c)][i];
        }
        size_type rank_extra(size_type i, const char_type c) const {
            return m_rank[1][char_index(c)](i);
        }

        select_1_type get_sel(const char_type c) const {
            return select_1_type(&m_bv[1][char_index(c)]);
        }
        
        value_type lf_value(size_type i) const {
           
            char_type ch; 
            size_type j, cnt = 0;
            // Case 1. Check whether the sum of the change of occ values is more than 1 or not
            if(i != 0){
                for(j = 1; j < sigma; ++j){
                    //select_0_type s = select_0_type(&m_bv[0][char_index(comp2char[j])]);
                    //size_type delta = s(i+1) - s(i) - 1;
                    // Case 1-2. If the sum is 1, 
                    // return the value calcalated from the char in which the count changes.
                    if(m_bv[0][char_index(comp2char[j])][i]){
                        ch = comp2char[j];
                        return C[j] + bwt.rank(i, ch);
                    }
                }
                
                // Case 1-1. The sum is 0. Check whether the extra value is equal to 1, 
                // If so, the current position is the starting position of a nonshared region.
                // Else, the char of the current position must be a unique final character, so return 0.
                for(j = 1; j < sigma; ++j){
                    ch = comp2char[j];
                    if(m_bv[1][char_index(ch)][i]){
                        return C[j] + bwt.rank(i+1, ch) - 1; 
                    }
                }
                return 0;
            }
            // Case 2. If the current position is 0,
            // find the char with nonzero occ value
            else{
                for(j = 1; j < sigma; ++j){
                    if(m_bv[0][j-1][0] != 0){
                        ch = comp2char[j];
                        //return C[j] + bwt.rank(i, ch);  
                        return C[j];  
                    }
                }
            }
            std::cout << "csaa.hpp: lf_value reaches end of non_void function" << std::endl;
            exit(1);
        }

        value_type lf_value(size_type strn, size_type i) const {
           
            //std::cout << "lf_value("<<strn<<", "<<i<<")"<<std::endl;
            char_type ch = -1; 
            size_type j, occCnt = 0, extCnt = 0;
            bool occur[sigma], extra[sigma];
            for(int k = 0; k < sigma; ++k){
                occur[k] = false;
                extra[k] = false;
            }
            // Check whether the sum of the change of occ is more than 1 or not
            if(i != 0){
                for(j = 1; j < sigma; ++j){
                    //select_0_type s = select_0_type(&m_bv[0][j-1]);
                    //size_type delta = s(i+1) - s(i) - 1;
                    size_type idx = char_index(comp2char[j]);
                    if(m_bv[0][idx][i]){
                        ch = comp2char[j];
                        ++occCnt;
                        occur[j] = true;
                    }
                    else if(m_bv[1][idx][i]){
                        ch = comp2char[j];
                        ++extCnt;
                        extra[j] = true;
                    }
                }
            }
            // If the current position is 0,
            // find the char with nonzero occ value
            else{
                for(j = 1; j < sigma; ++j){
                    size_type lf_val    = C[j] + bwt.rank(i+1, comp2char[j])-1;
                    size_type pair      = m_saa_sample.saa_value(lf_val);
                    size_type div       = pair/m_text_size;
                    // std::cout <<"Case 2. LF["<<i<<"]= "<<C[x]<<"+"<<bwt.rank(i+1, comp2char[x])<<"-1"<<std::endl;
                    // std::cout <<"Case 2. POS["<<lf_val<<"]: ("<<pair/m_text_size<<","<<pair%m_text_size<<"), pair: "<<pair<<", text_size: "<<m_text_size<<std::endl;
                    if(div == strn){
                        return lf_val;
                    }
                    else if(m_saa_sample.strv[div*num_string+strn-1]){ 
                        return lf_val;
                    }
                    else if(strn == 0){
                        std::cout << "csaa.hpp: lf_value error" << std::endl;
                        exit(1);
                    } 
                }
            }
            
            // Case 1-1. The occurrence sum is 0. Check whether the extra value is equal to 1, 
            // If so, the current position is the starting position of a nonshared region.
            // Else, the char of the current position must be a unique final character, so return 0.
            // Case 1-2. If the occurrence sum is 1, 
            // return the value calcalated from the char in which the occ changes.
            if(occCnt + extCnt == 1){
                // std::cout <<"Case 1. LF["<<i<<"]= "<<C[char2comp[ch]]<<"+"<<bwt.rank(i+1, ch)<<"-1"<<std::endl;
                return C[char2comp[ch]] + bwt.rank(i+1, ch) - 1; 
            }
            else if(occCnt == 0 and extCnt == 0){
                return 0;
            }
            // Case 2. If the occurrence sum is more than 1 and extra sum is 0,
            // the current position is the starting position of a shared region or in a non-shared region.
            else if(occCnt != 0 and extCnt == 0){
                for(size_type x = 1; x < sigma; ++x){
                    if(occur[x]){ 
                        size_type lf_val    = C[x] + bwt.rank(i+1, comp2char[x])-1;
                        size_type pair      = m_saa_sample.saa_value(lf_val);
                        size_type div       = pair/m_text_size;
                        // std::cout <<"Case 2. LF["<<i<<"]= "<<C[x]<<"+"<<bwt.rank(i+1, comp2char[x])<<"-1"<<std::endl;
                        // std::cout <<"Case 2. POS["<<lf_val<<"]: ("<<pair/m_text_size<<","<<pair%m_text_size<<"), pair: "<<pair<<", text_size: "<<m_text_size<<std::endl;
                        if(div == strn){
                            return lf_val;
                        }
                        else if(m_saa_sample.strv[div*num_string+strn-1]){ 
                            return lf_val;
                        }
                        else if(strn == 0){
                            std::cout << "csaa.hpp: lf_value error" << std::endl;
                            exit(1);
                        } 
                    }
                }
            }
            // Case 3. If the occurrence sum is more than 1 and extra sum is non-zero.
            // the current position is the starting position of a shared region or in a non-shared region.
            else{
                size_type shared_lf_val = 0, lf_val, pair, div;
                for(size_type x = 1; x < sigma; ++x){
                    if(occur[x] or extra[x]){ 
                        lf_val    = C[x] + bwt.rank(i+1, comp2char[x])-1;
                        pair      = m_saa_sample.saa_value(lf_val);
                        div       = pair/m_text_size;
                        // std::cout <<"Case 3. LF["<<i<<"]= "<<C[x]<<"+"<<bwt.rank(i+1, comp2char[x])<<"-1"<<std::endl;
                        // std::cout <<"Case 3. POS["<<lf_val<<"]: ("<<pair/m_text_size<<","<<pair%m_text_size<<"), pair: "<<pair<<", text_size: "<<m_text_size<<std::endl;
                        if(div == strn){
                            return lf_val;
                        }
                        else if(div > num_string and m_saa_sample.strv[div*num_string+strn-1]){ 
                            return lf_val;
                        }
                        else if(strn == 0){
                            std::cout << "csaa.hpp: lf_value error" << std::endl;
                            exit(1);
                        } 

                        if(extra[x] and div == 0) 
                            shared_lf_val = lf_val;
                    }
                }
                return shared_lf_val;
            }
            std::cout << "csaa.hpp: lf_value reaches end of non_void function" << std::endl;
            exit(1);
        }
        /*
        char_type bwt_value(size_type i) const {
            char_type result = 0;
            size_type j;
            for(j = 0; j < conf::CHAR_NUM; ++j){
                select_0_type s = select_0_type(&m_bv[0][j]);
                if((i == 0 and m_bv[0][j][0] != 0) or (i != 0 and s(i+1) - s(i) > 1))
                    break;
            }
            if(j == conf::CHAR_NUM){
                for(j = 0; j < conf::CHAR_NUM; ++j){
                    if(m_bv[1][j][i])
                        break;
                }
            }
            if(j == conf::CHAR_NUM)
                return result;
            else
                return conf::CHAR[j];
        }
        */
        void print() const {    
            std::string category[2];
            category[0] = "occ";
            category[1] = "extra";
            std::cout << "alphabet requires "  << size_in_mega_bytes(m_alphabet) << "MiB." << std::endl;
            double size = 0;
            for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
                for(size_type j = 0; j < conf::CHAR_NUM; ++j){
                    //std::cout << category[i] << conf::CHAR[j] << " requires " << size_in_mega_bytes(m_bv[i][j]) << "MiB." << std::endl;
                    //std::cout << "rank(" << category[i] << conf::CHAR[j] << ") requires " << size_in_mega_bytes(m_rank[i][j]) << "MiB." << std::endl; 
                    size += size_in_mega_bytes(m_bv[i][j]);
                    size += size_in_mega_bytes(m_rank[i][j]);
                }
            }
            std::cout << "occurrences and extras requires " << size << "MiB." << std::endl;
        }
        

        // Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *        \f$ \Order{\log n t_{\Psi}} \f$
         */
        size_type rank_bwt(size_type i, const char_type c)const {

            comp_char_type cc = char2comp[c];
            if (cc==0)  // character is not in the text or character is zero_symbol => return 0
                return 0;
            if (i == 0)
                return 0;
            // std::cout << i << " " << size() << std::endl;
            // assert(i < size());
            //select_0_type select_zero(&m_bv[0][cc-1]);
            //size_type temp = select_zero(i);
            size_type res = m_rank[0][char_index(c)](i);
            //size_type res = m_rank[0][cc-1]((i==0)?0:select_zero(i));
            //std::cout << "CSAA.rank_bwt(" << select_zero(i+1) << ", " << (char)c << ")=" << "rank[0][" << cc-1 << "]: " << res << std::endl;
            //std::cout << "CSAA.rank_bwt(" << i << ", " << (char)c << ") =" << "rank[0][" << cc-1 << "]: " << res << std::endl;
            
            return res; 

        }
};

// == template functions ==

template<class t_bitvector, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::csaa(cache_config& ref_config, str_vector_type& strv, bit_vector& is_non_shared, bit_vector& sampled, bit_vector_type vector[][conf::CHAR_NUM], gap_type& gap_vector, const size_type& number)
{
    //create_buffer();
    /*
    if (!cache_file_exists(key_trait<alphabet_type::int_width>::KEY_BWT, ref_config)) {
        return;
    }
    */
    m_text_size = is_non_shared.size();
    m_num_string = number;
    for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
        for(size_type j = 0; j < conf::CHAR_NUM; ++j){
            m_bv[i][j].swap(vector[i][j]);
            m_rank[i][j] = rank_1_type(&(m_bv[i][j]));
        }
    }
    m_gap_vector = gap_vector;

    //int_vector_buffer<alphabet_type::int_width> bwt_ref(cache_file_name(key_trait<alphabet_type::int_width>::KEY_BWT,ref_config));
    //int_vector_buffer<alphabet_type::int_width> bwt_ss(cache_file_name(conf::KEY_BSS, ref_config));
    {
        auto event = memory_monitor::event("construct csa-alpbabet");
        alphabet_type tmp_alphabet(m_bv[0], m_rank[0]);
        m_alphabet.swap(tmp_alphabet);
    }

    //int_vector<> cnt_chr(sigma, 0, bits::hi(n)+1);
    // std::cout << "C[i]" << std::endl;
    // for (typename alphabet_type::sigma_type i=0; i <= sigma; ++i) {
    //     cnt_chr[i] = C[i];
    //     std::cout << i << ":" << C[i] << std::endl;
    // }

    // calculate lf
    /*
    {
        auto event = memory_monitor::event("construct LF");
        //m_has_multi_val = bit_vector(n, 0);
        bit_vector::rank_1_type shared_rank(&shared_start);
        //int_vector<> lf(n, 0, bits::hi(n)+1);
        //int_vector<> lf_ss(shared_rank(shared_start.size())*m_num_string, 0, bits::hi(n)+1);
        //size_type lf_ss_cnt = 0;
        int_vector<64> non_shared_cnt(conf::CHAR_NUM+1, 0);
       
        //std::cout << "IDX\tBWT\tcc\tLF" << std::endl;
        for (size_type i=0; i < n; ++i) {
            size_type strn      = saa_buf[2*i];
            size_type pos       = saa_buf[2*i+1];
            comp_char_type cc   = char2comp[bwt_ref[i]];
            if(non_shared_start[pos]){
                //lf[i] = cnt_chr[cc];
                //std::cout <<  i << "\t" << (char)bwt_ref[i] << "\t" << (uint16_t)cc << "\t" << lf[i] << std::endl;
                non_shared_cnt[cc]++;
                if(non_shared_cnt[cc]%m_num_string == 0){
                    non_shared_cnt[cc] = 0;
                    cnt_chr[cc]++;
                }
            }
            else if(shared_start[pos]){
                //m_has_multi_val[i] = 1;
                //lf[i] = cnt_chr[ cc ];
                //std::cout << i << "\t" << (char)bwt_ref[i] << "\t" << (uint16_t)cc << "\t" << lf[i] << std::endl;
                for(size_type j = 0; j < m_num_string; ++j){
                    comp_char_type ccss = char2comp[bwt_ss[shared_rank(pos)*3+j]];
                    //lf_ss[lf_ss_cnt] = cnt_chr[char2comp[bwt_ss[shared_rank(pos)*3+j]]]++;
                    //std::cout << "\t\tBSS[" << j << "]: " << (char)bwt_ss[shared_rank(pos)*3+j] << "(cc: " << (uint16_t)ccss << ") | " << "LF_SS[" << j << "]:" << lf_ss[lf_ss_cnt] << std::endl;
                    //lf_ss_cnt++;
                }
            }
            else{
                //lf[i] = cnt_chr[ cc ];
                cnt_chr[ cc ]++;
                //std::cout << i << "\t" << (char)bwt_ref[i] << "\t" << (uint16_t)cc << "\t" << lf[i] << std::endl;
            }
        }
        */
        /*
        std::cout << "Has_multi_val" << std::endl;
        for(size_type i=0; i<n; ++i)
            std::cout << i << ": " << m_has_multi_val[i] << std::endl;
        */

        /*
        std::string lf_file = cache_file_name(conf::KEY_LF, ref_config);
        if (!store_to_cache(lf, conf::KEY_LF, ref_config)) {
            return;
        }
        */
        /*
        std::string lf_ss_file = cache_file_name(conf::KEY_LF_SS, ref_config);
        if (!store_to_cache(lf_ss, conf::KEY_LF_SS, ref_config)) {
            return;
        }
        
    }
    */
    /*
    {
        
        auto event = memory_monitor::event("encode LF");
        //int_vector_buffer<> lf_buf(cache_file_name(conf::KEY_LF, ref_config));
        //enc_vector_type tmp_lf(lf_buf);
        //m_lf.swap(tmp_lf);
        int_vector_buffer<> lf_ss_buf(cache_file_name(conf::KEY_LF_SS, ref_config));
        enc_vector_type tmp_lf_ss(lf_ss_buf);
        m_lf_ss.swap(tmp_lf_ss);
        
    }
    */
    {
        auto event = memory_monitor::event("sample SAA");
        saa_sample_type tmp_saa_sample(ref_config, strv, sampled, m_num_string, m_text_size);
        m_saa_sample.swap(tmp_saa_sample);
        
        /* 
        for(size_type i = 0; i < m_saa_sample.size(); ++i){
            size_type val = m_saa_sample[i];
            std::cout << i << ": (" << val/m_text_size << ", " << val%m_text_size << ")" << std::endl; 
        }
        */
        
    }
    {
        
        auto event = memory_monitor::event("sample ISA");
        int_vector_buffer<>  saa_buf(cache_file_name(conf::KEY_SA, ref_config));
        //bit_vector is_non_shared(non_shared_start.size(), 0);
        //bool non_shared = false;
        m_saa_size = saa_buf.size()/2;
        //bit_vector::rank_1_type r(&non_shared_start);
        //std::cout << "# of Non-shared region: " << r(non_shared_start.size()) << std::endl;
        /*
        for(size_type i = 0; i < is_non_shared.size(); ++i){
            if(shared_start[i])
                non_shared = false;
            else if(non_shared_start[i])
                non_shared = true;
            if(non_shared)
                is_non_shared[i] = 1;
        }
        */
        //bit_vector::rank_1_type rr(&is_non_shared);
        //std::cout << "# of Non-shared positions: " << rr(is_non_shared.size()) << std::endl;
        
        std::cout << "Text size: " << m_text_size << std::endl; 
        std::cout << "SAA size: " << m_saa_size << std::endl; 
        std::cout << "SAA sample size " << m_saa_sample.sample_size << std::endl;
        
        /*
        std::cout << "is_non_shared" << std::endl;
        for(size_type i = 0; i < is_non_shared.size(); ++i){
            std::cout << i << ": " << is_non_shared[i] << std::endl; 
        }
        */
        //std::map<size_type, size_type> map1, map2;
        //set_isaa_samples<csaa>(saa_buf, m_index_table, m_ns_isaa_sample, is_non_shared, m_text_size, m_num_string, strv, gap_file, map1);
        //bit_vector isns = build_isaa_samples<csaa>(m_saa_sample, m_index_table2, m_ns_isaa_sample2, m_text_size, m_num_string, strv, gap_file, map2);
        build_isaa_samples<csaa>(m_saa_sample, m_index_table, m_ns_isaa_sample, m_text_size, m_num_string, m_gap_vector);
        std::cout <<"csaa.hpp: Built isaa samples"<<std::endl;
        
        /*
        std::cout << "compare aligned_to_saa_idx: " << map1.size()<<", " <<map2.size()<<std::endl;
        std::map<size_type, size_type>::iterator it2 = map2.begin();
        for(std::map<size_type, size_type>::iterator it = map1.begin(); it != map1.end(); ++it, ++it2){
            std::cout<<it->first/(m_num_string+1)<<","<<it->first%(m_num_string+1)<<" - "<<it->second<<"\t"<<it2->first/m_num_string<<","<<(it2->first%m_num_string)+1<<" - "<<it2->second<<std::endl;
        }
        std::cout << "is_non_shared != isns" << std::endl; 
        for(size_type i = 0; i < isns.size(); ++i){
            size_type pos = (i*isa_sample_dens>text_size)?text_size-1:i*isa_sample_dens;
            if(is_non_shared[pos]!=isns[i])
                std::cout << "["<<pos<<", "<<i<<"]: "<<is_non_shared[pos]<<"\t"<<isns[i]<<std::endl;
        }
        */
        /*
        size_type m = saa_buf.size();
        std::cout  << "ISA sample" << std::endl;
        for(size_type i = 0; i < m_index_table.size(); ++i){
            size_type pos;
            if(i == m_index_table.size()-1 && i*isa_sample_dens > m_text_size)
                pos = m_text_size - 1;
            else
                pos = i*isa_sample_dens;
            size_type idx = m_index_table[i];
            if(idx < m){
                std::cout << "(0, " << pos << "): " << idx << std::endl;
            }
            else{
                for(size_type i = 0; i < m_num_string; ++i){
                    std::cout << "(" << i+1 << ", " << pos << "): " << m_ns_isaa_sample[m_num_string*(idx-m)+i] << std::endl;
                }
            }
        }
        */
    }
}

//inline auto csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::access(size_type i, size_type& mapping, double& time)const -> pair_type
template<class t_bitvector, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
inline auto csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::operator[](size_type i)const -> pair_type
{
    size_type off = 0;
    value_type v = 0;
    while (!m_saa_sample.is_sampled(i)) { // while i mod t_dens != 0 (SA[i] is not sampled)   SG: auf keinen Fall get_sample_dens nehmen, ist total langsam
        v = lf_value(i);       // go to the position where SAA[i]-1 is located
        ++off;              // add 1 to the offset
        //std::cout << "      LF[" << i << "]: " << v << std::endl;
        i = v;
    }
    size_type  pair     = m_saa_sample.saa_value(i);
    value_type strn     = pair/m_text_size;
    value_type result   = pair%m_text_size;
    if (result + off >= m_text_size) {
        return std::make_pair(strn, off + result - m_saa_size);
    } else{
        return std::make_pair(strn, result + off);
    }
}


    template<class t_bitvector, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
auto csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const -> size_type
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    //written_bytes += m_lf.serialize(out, child, "lf");
    //written_bytes += m_lf_ss.serialize(out, child, "lf_ss");
    written_bytes += m_gap_vector.serialize(out, child, "gap_vector");
    written_bytes += m_saa_sample.serialize(out, child, "saa_samples");
    //written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
    
    // 160601 deleted
    //written_bytes += m_ns_isaa_sample.serialize(out, child, "m_ns_isa_samples");
    //written_bytes += m_index_table.serialize(out, child, "m_index_table");
    
    // 160603 add for verification
    //written_bytes += m_ns_isaa_sample2.serialize(out, child, "m_ns_isa_samples2");
    //written_bytes += m_index_table2.serialize(out, child, "m_index_table2");
    
    written_bytes += m_alphabet.serialize(out, child, "alphabet");
    //written_bytes += m_shared_start.serialize(out, child, "shared_start");
    //written_bytes += m_non_shared_start.serialize(out, child, "non_shared_start");
    //written_bytes += m_has_multi_val.serialize(out, child, "has_multi_val");
    written_bytes += write_member(m_num_string, out, child, "num_string");
    written_bytes += write_member(m_text_size, out, child, "text_size");
    written_bytes += write_member(m_saa_size, out, child, "saa_size");
    const std::string category[2] = {"count", "rank"};
    const std::string rank_string = "rank_";
    std::string temp;
    for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
        for(size_type j = 0; j < conf::CHAR_NUM; ++j){
            temp = category[i] + std::to_string(j);
            written_bytes += m_bv[i][j].serialize(out, child, temp);
            temp = rank_string + temp;
            written_bytes += m_rank[i][j].serialize(out, child, temp); 
        }
    }
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_bitvector, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::load(std::istream& in)
{
    //const std::string gap_name = "gapinfo.vec"; 
	//std::string gap_file = info_dir + gap_name;
    m_gap_vector.load(in);
    //m_lf.load(in);
    //m_lf_ss.load(in);
    m_saa_sample.load(in);
    //m_isa_sample.load(in);
    
    // 160601 deleted
    //m_ns_isaa_sample.load(in);
    //m_index_table.load(in);
    
    // 160603 add for verification
    //m_ns_isaa_sample2.load(in);
    //m_index_table2.load(in);
    
    m_alphabet.load(in);
    //m_shared_start.load(in);
    //m_non_shared_start.load(in);
    //m_has_multi_val.load(in);
    read_member(m_num_string, in);
    read_member(m_text_size, in);
    read_member(m_saa_size, in);
    for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
        for(size_type j = 0; j < conf::CHAR_NUM; ++j){
            m_bv[i][j].load(in);
            m_rank[i][j].load(in, &(m_bv[i][j]));
        }
    }
    
    // 160601 build a ns_isaa_sample and index table
    //build_isaa_samples<csaa>(m_saa_sample, m_index_table, m_ns_isaa_sample, m_text_size, m_num_string, m_saa_sample.strv, gap_file);
    build_isaa_samples<csaa>(m_saa_sample, m_index_table, m_ns_isaa_sample, m_text_size, m_num_string, m_gap_vector);
    std::cout <<"csaa.hpp: Built isaa samples"<<std::endl;
}

template<class t_bitvector, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::swap(csaa<t_bitvector, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa)
{
    if (this != &csa) {
        //m_lf.swap(csa.m_lf);
        //m_lf_ss.swap(csa.m_lf_ss);
        m_gap_vector.swap(csa.m_gap_vector);
        m_saa_sample.swap(csa.m_saa_sample);
        //m_isa_sample.swap(csa.m_isa_sample);
        m_ns_isaa_sample.swap(csa.m_ns_isaa_sample);
        m_index_table.swap(csa.m_index_table);
        
        // 160603 add for verification
        //m_ns_isaa_sample2.swap(csa.m_ns_isaa_sample2);
        //m_index_table2.swap(csa.m_index_table2);
        
        m_alphabet.swap(csa.m_alphabet);
        //m_shared_start.swap(csa.m_shared_start);
        //m_non_shared_start.swap(csa.m_non_shared_start);
        //m_has_multi_val.swap(csa.m_has_multi_val);
        std::swap(m_num_string, csa.m_num_string);
        std::swap(m_text_size, csa.m_text_size);
        std::swap(m_saa_size, csa.m_saa_size);
        for(size_type i = 0; i < conf::CATEGORY_NUM; ++i){
            for(size_type j = 0; j < conf::CHAR_NUM; ++j){
                m_bv[i][j].swap(csa.m_bv[i][j]);
                util::swap_support(m_rank[i][j], csa.m_rank[i][j],
                                   &m_bv[i][j], &(csa.m_bv[i][j]));
                //m_rank[i][j].swap(csa.m_rank[i][j]);
            }
        }
    }
}

} // end namespace sdsl
