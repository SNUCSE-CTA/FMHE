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
/*! \file construct_fma.hpp
    \brief construct_fma.hpp contains a construction method for the Burrows Wheeler Transform (BWT) of Suffix Array of Alignment (SAA), FM-index of alignments with gaps, and includes the strategy class for an alphabet and SAA sampling in an SAA.
    \author Hyunjoon Kim
    \reference 
    [1] J. Na, H. Kim, H. Park, T. Lecroq, M. Le'onard, L. Mouchar, and K. Park, FM-index of alignment: A compressed index for similar strings, TCS 638:159-170, 2016
    [2] J. Na, H. Kim, S. Min, H. Park, T. Lecroq, M. Le'onard, L. Mouchard, and K. Park, FM-index of alignment with gaps, TCS 710:148-157, 2018
    \date 2016.09
*/

#include <string>
#include <iostream>
#include <stdexcept>
#include <list>
#include <algorithm>
#include <vector>
#include <cctype>
#include <set>

#include <sdsl/uintx_t.hpp> // for uint8_t
#include <sdsl/config.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include "divsufsort.h"
#include "divsufsort64.h"
#include <sdsl/qsufsort.hpp>
#include <sdsl/construct_sa_se.hpp>
#include <sdsl/construct_config.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_bwt.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/util.hpp>
#include <sdsl/sfstream.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/csa_alphabet_strategy.hpp> // for key_trait

namespace sdsl
{
struct csaa_tag {};// compressed suffix array of alignment (CSAAs) tag
struct lf_csaa_tag {}; // tag for CSAAs based on the LF function
struct saa_alphabet_tag { static const uint8_t WIDTH=8; };

namespace conf  // namespace for library constant
{
// size of the buffer for reading and writing data in elements (not in bytes)

const char KEY_BSS[] 		= "bwt_ss"; //BWT for shared starting positions
//const char KEY_ISA[] 		= "isa";
const char KEY_LF[] 		= "lf";
const char KEY_LF_SS[] 		= "lf_ss";
const uint8_t CATEGORY_NUM  = 2;
const uint8_t CHAR_NUM  = 5;
const char CHAR[5] = {'A', 'C', 'G', 'N', 'T'};
}

//! Helper classes to transform width=0 and width=8 to corresponding bss key
template<uint8_t width>
struct key_bwt_ss_trait {
    static const char* KEY_BSS;
};

//! Internal function used by csXprintf
template<class t_csa>
const t_csa& _idx_csa(const t_csa& t, csaa_tag)
{
    return t;
};

//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
 *         * conf::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * conf::KEY_BWT for t_width=8 or conf::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt_of_saa(cache_config& ref_config, cache_config* str_config, bit_vector& shared_start, int_vector<>::size_type num_string)
{
    static_assert(t_width == 0 or t_width == 8 , "construct_bwt: width must be `0` for integer alphabet and `8` for byte alphabet");

    typedef int_vector<>::size_type size_type;
    typedef int_vector<t_width> text_type;
    typedef int_vector_buffer<t_width> bwt_type;
    const char* KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    const char* KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;
    //const char* KEY_BSS = key_bwt_ss_trait<t_width>::KEY_BSS;
	
    //  (1) Load text from disk
    text_type text[num_string];
    for(size_type i = 0; i < num_string; ++i)
    	load_from_cache(text[i], KEY_TEXT, str_config[i]);
    size_type n = text[0].size();
    uint8_t bwt_width = text[0].width();
    
    //  (2) Prepare to stream SA from disc and BWT to disc
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
    int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, ref_config), std::ios::in, buffer_size);
    std::string bwt_file = cache_file_name(KEY_BWT, ref_config);
    bwt_type bwt_ref(bwt_file, std::ios::out, buffer_size, bwt_width);
	
    bwt_type bwt_ss;
    for(size_type i = 0; i < num_string; ++i){
	bwt_file = cache_file_name(conf::KEY_BSS, ref_config);
	bwt_ss = bwt_type(bwt_file, std::ios::out, buffer_size, bwt_width);
    }

    //  (3) Construct BWT sequentially by streaming SA and random access to text
    size_type to_add[2] = {(size_type)-1,n-1};
    size_type cnt = 0;
    std::cout << "Construct BWT..." << std::endl;
    for (size_type i=0; i < sa_buf.size()/2; ++i) {
	size_type strn 	= sa_buf[2*i];
	size_type pos 	= sa_buf[2*i+1];
	if(strn != 0){
	    bwt_ref[i] = text[strn-1][pos+to_add[pos==0]];
	}
	else if(strn == 0 and !shared_start[pos]){
	    bwt_ref[i] = text[strn][pos+to_add[pos==0]];
	}
	else{
            bwt_ref[i] = text[strn][pos+to_add[pos==0]];
	    for(size_type j = 0; j < num_string; ++j){
	        bwt_ss[cnt++] = text[j][pos+to_add[pos==0]];
	    }
	}
    }
    bwt_ss.close();
    bwt_ref.close();
    register_cache_file(conf::KEY_BSS, ref_config);
    register_cache_file(KEY_BWT, ref_config);	
}

template<uint8_t t_width>
void construct_saa(cache_config& config, const std::string& saa_file)
{
	static_assert(t_width == 0 or t_width == 8 , "construct_saa: width must be `0` for integer alphabet and `8` for byte alphabet");
	int_vector<> saa;
	load_from_file(saa, saa_file);
	store_to_cache(saa, conf::KEY_SA, config);
}


/*
template<class t_index>
void construct(t_index& idx, std::string file, uint8_t num_bytes=0)
{
    tMSS file_map;
    cache_config config;
    if (is_ram_file(file)) {
        config.dir = "@";
    }
    construct(idx, file, config, num_bytes);
}
*/


//! Constructs an index object of type t_index for a text stored on disk.
/*!
 * \param idx       	t_index object.  Any sdsl suffix array of suffix tree.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */
/*
template<class t_index>
void construct(t_index& idx, const std::string& file, cache_config& config, uint8_t num_bytes=0)
{
    // delegate to CSA or CST construction
    typename t_index::index_category 		index_tag;
    construct(idx, file, config, num_bytes, index_tag);
}
*/


std::string get_text_file_name(uint16_t n)
{
	std::string nn = std::to_string(n);
	int_vector_size_type len = nn.length();
	for(int_vector_size_type i = 0; i < 3 - len; ++i)
		nn = "0" + nn;
	return "S" + nn;
}

// Specialization for SAAs
template<class t_index>
void construct(t_index& idx, const std::string& info_dir, cache_config& config, uint8_t num_bytes, csaa_tag)
{
    const typename t_index::size_type num_string = idx.num_string;
    const std::string vec_suffix = ".vec";
    const std::string saa_name = "finSaa";
    const std::string alpha_name = "alpha";
    const std::string non_shared_name = "nonshared" ;
    const std::string count_name = "occ";
    const std::string extra_name = "ext";
    const std::string str_vec_name = "strv";
    const std::string sampled_name = "sampled";
    const std::string cnt_name = "cnt";
    const std::string gap_name = "gapinfo"; 
    const std::string num_string_name = std::to_string(num_string);
    const std::string dens_string = std::to_string(t_index::sa_sample_dens);
    const std::string char_name[5] = {"A", "C", "G", "N", "T"};
	
    std::string non_shared_file = info_dir + non_shared_name + num_string_name + vec_suffix;
    std::string alpha_file = info_dir + alpha_name + num_string_name + vec_suffix;
    std::string cnt_file = info_dir + cnt_name + num_string_name + vec_suffix;
    std::string sampled_file = info_dir + dens_string + sampled_name + num_string_name + vec_suffix;
    std::string gap_file = info_dir + gap_name + vec_suffix;
    int_vector<> non_shared, alpha, C;
    bit_vector picked;
    load_from_file(non_shared, non_shared_file);
    load_from_file(alpha, alpha_file);
    load_from_file(C, cnt_file);
    load_from_file(picked, sampled_file);
    

    typedef typename t_index::str_vector_type str_vector_type;
    typedef typename t_index::gap_type gap_type;
    auto event = memory_monitor::event("construct SAA");
    int_vector_size_type n;
    n = idx.text_size;
    std::cout << "SIZE(text)= " << n << std::endl;

    // (1) make is_non_shared
    int_vector_size_type dens = t_index::sa_sample_dens;
    typename t_index::size_type sample_cnt = 0;
    std::cout << "#(NonCommon): " << non_shared.size()/2 << std::endl;
    std::cout << "Sampling dens: " << dens << std::endl;
    std::cout << "#(String): " << num_string << std::endl;
    
    bit_vector is_non_shared(n, 0);
    for(int_vector_size_type i = 0; i < non_shared.size()/2;++i){
        int_vector_size_type start = non_shared[2*i];
        int_vector_size_type end = non_shared[2*i+1];
        for(int_vector_size_type j = start; j <= end; ++j){
            is_non_shared[j] = 1;
        }
    }

    {
        // (2) check, if the suffix array is cached
        auto event = memory_monitor::event("SAA");
		std::string saa_file = info_dir + saa_name + num_string_name + vec_suffix;
        if (!cache_file_exists(conf::KEY_SA, config)) {
            construct_saa<t_index::alphabet_category::WIDTH>(config, saa_file);
        }
        register_cache_file(conf::KEY_SA, config);
    }
    {
        //  (3) load bit vectors
        auto event = memory_monitor::event("construct CSAA");
        const int_vector_size_type row_size = conf::CATEGORY_NUM, col_size = conf::CHAR_NUM;
        typename t_index::bit_vector_type bv[row_size][col_size];
        for(int_vector_size_type i = 0; i < row_size; ++i){
            for(int_vector_size_type j = 0; j < col_size; ++j){
                std::string category_name = i?extra_name:count_name;
                std::string bv_path = info_dir + category_name + num_string_name + char_name[j] + vec_suffix;
                load_from_file(bv[i][j], bv_path);
            }
        }
          
	std::string str_vec_file = info_dir + str_vec_name + num_string_name + vec_suffix;
        str_vector_type strv;
        gap_type gap_vector;
        load_from_file(strv, str_vec_file);
        load_from_file(gap_vector, gap_file);
        t_index tmp(config, strv, is_non_shared, picked, bv, gap_vector, num_string);
        std::cout << "Construct: Finished making a csaa." << std::endl;
        idx.swap(tmp);
    }
    if (config.delete_files) {
        auto event = memory_monitor::event("delete temporary files");
        util::delete_all_files(config.file_map);
    }
}
};

// TODO: Strategy with 1-to-1 mapping and C_array type as template parameter
//       This can be also used for a integer based CSA.

/* A alphabet strategy provides the following features:
 *   * Member `sigma` which contains the size (=umber of unique symbols) of the alphabet.
 *   * Method `is_present(char_type c)` which indicates if character c occurs in the text.
 *   * Container `char2comp` which maps a symbol to a number [0..sigma-1]. The alphabetic
 *     order is preserved.
 *   * Container `comp2char` which is the inverse mapping of char2comp.
 *   * Container `C` contains the cumulative counts of occurrences. C[i] is the cumulative
 *     count of occurrences of symbols `comp2char[0]` to `comp2char[i-1]` in the text.
 *   * Typedefs for the four above members:
 *       * char2comp_type
 *       * comp2char_type
 *       * C_type
 *       * sigma_type
 *   * Constructor. Takes a int_vector_buffer<8> for byte-alphabets
 *     and int_vector_buffer<0> for integer-alphabets.
 *
 *    \par Note
 *   sigma_type has to be large enough to represent the alphabet size 2*sigma,
 *   since there is code which will perform a binary search on array `C`.
 */

namespace sdsl
{

template<class bit_vector_type     = rrr_vector<127>,
         class rank_rupport_type   = typename rrr_vector<127>::rank_1_type>
class saa_alphabet;

template<class bit_vector_type, class rank_support_type>
class saa_alphabet
{
    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<8>           char2comp_type;
        typedef int_vector<8>           comp2char_type;
        typedef int_vector<64>          C_type;
        typedef uint16_t                sigma_type;
        typedef uint8_t                 char_type;
        typedef uint8_t                 comp_char_type;
        typedef std::string             string_type;
        enum { int_width = 8 };

        typedef saa_alphabet_tag       alphabet_category;
    private:
        char2comp_type m_char2comp; // Mapping from a character into the compact alphabet.
        comp2char_type m_comp2char; // Inverse mapping of m_char2comp.
        C_type         m_C;         // Cumulative counts for the compact alphabet [0..sigma].
        sigma_type     m_sigma;     // Effective size of the alphabet.

        void copy(const saa_alphabet& saaa)
        {
            m_char2comp  = saaa.m_char2comp;
            m_comp2char  = saaa.m_comp2char;
            m_C          = saaa.m_C;
            m_sigma      = saaa.m_sigma;
        }
    public:

        const char2comp_type& char2comp;
        const comp2char_type& comp2char;
        const C_type&         C;
        const sigma_type&     sigma;

        //! Default constructor
        saa_alphabet(): char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
        {
            m_sigma = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        saa_alphabet(const bit_vector_type count[conf::CHAR_NUM], const rank_support_type count_rank[conf::CHAR_NUM]):
            char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
        {
            util::assign(m_C, int_vector<64>(257, 0));
            util::assign(m_char2comp, int_vector<8>(256, 0));
            util::assign(m_comp2char, int_vector<8>(256, 0));
            
            size_type len = 0; 
            for(size_type i = 0; i < conf::CHAR_NUM; ++i){
                m_C[conf::CHAR[i]] = count_rank[i](count[i].size());
                len += m_C[conf::CHAR[i]];
            }
            m_C[0] = 1;
            m_sigma = 0;
            for(size_type i = 0; i < 256; ++i){
                if(m_C[i]){
                    m_char2comp[i]      = m_sigma;
                    m_comp2char[sigma]  = i;
                    m_C[m_sigma]        = m_C[i];
                    if(i != 0)
                    std::cout << m_sigma << ":" << m_C[m_sigma] << " (<-m_C[" << i << "])" <<std::endl;
                    ++m_sigma;
                }
            }
            m_comp2char.resize(m_sigma);
            m_C.resize(m_sigma+1);
            for(size_type i = m_sigma; i > 0; --i){  
                m_C[i] = m_C[i-1];   
            }

            m_C[0] = 0;
            for(size_type i = 1; i <= m_sigma; ++i) m_C[i] += m_C[i-1];
            // assert(C[sigma]==len);

            std::cout << "m_C[i](final)" << std::endl;
            for(size_type i = 0; i <= m_sigma; i++)
                std::cout << i << ": " << m_C[i] << std::endl;

            std::cout << "Comp2Char..." << std::endl;
            for(size_type i = 0; i < m_sigma; i++)
                std::cout << i << ": " << (uint16_t) m_comp2char[i] << std::endl;
        }

        saa_alphabet(const saa_alphabet& saaa): char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
        {
            copy(saaa);
        }

        saa_alphabet& operator=(const saa_alphabet& saaa)
        {
            if(this != &saaa){
                copy(saaa);
            }
            return *this;
        }

        void swap(saa_alphabet& saaa)
        {
            m_char2comp.swap(saaa.m_char2comp);
            m_comp2char.swap(saaa.m_comp2char);
            m_C.swap(saaa.m_C);
            std::swap(m_sigma, saaa.m_sigma);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));

            size_type written_bytes = 0;
            written_bytes += m_char2comp.serialize(out, child, "m_char2comp");
            written_bytes += m_comp2char.serialize(out, child, "m_comp2char");
            written_bytes += m_C.serialize(out, child, "m_C");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            m_char2comp.load(in);
            m_comp2char.load(in);
            m_C.load(in);
            read_member(m_sigma, in);
        }
};

template<class t_csa,
         class bit_vector_type=bit_vector,
         //class bit_vector_type=rrr_vector<127>,
         class rank_support_type=typename bit_vector_type::rank_1_type,
         uint8_t t_width=0
         >
class _saa_sampling : public int_vector<t_width>
{
    public:
        typedef int_vector<t_width> base_type;
        typedef rrr_vector<127>     marked_type;
        typedef typename base_type::size_type  size_type;	// make typedefs of base_type visible
        typedef typename base_type::value_type value_type;	//
        typedef typename t_csa::str_vector_type  str_vector_type;
        enum { sample_dens = t_csa::sa_sample_dens };
    
    private:
        bit_vector_type		m_marked;
        rank_support_type	m_rank_marked;
        str_vector_type     m_strv;
        size_type           m_sample_size; 
    public:
        const bit_vector_type&      marked      = m_marked;
        const rank_support_type&    rank_marked = m_rank_marked;
        const str_vector_type&      strv        = m_strv;
        const size_type&            sample_size = m_sample_size;
        //! Default constructor
        _saa_sampling() {}

        //! Constructor
        /*
         * \param cconfig Cache configuration (SA is expected to be cached.).
         * \param csa    Pointer to the corresponding CSA. Not used in this class.
         * \par Time complexity
         *      Linear in the size of the suffix array.
         */
        //_saa_sampling(const cache_config& cconfig, bit_vector& sampled, size_type& sample_cnt, size_type& num, SDSL_UNUSED const t_csa* csa=nullptr) 
        _saa_sampling(const cache_config& cconfig, str_vector_type& str_vector, bit_vector& sampled, size_type& num, size_type& text_size, SDSL_UNUSED const t_csa* csa=nullptr) 
        {
            m_marked = sampled;
            util::init_support(m_rank_marked, &m_marked);
            int_vector_buffer<>  sa_buf(cache_file_name(conf::KEY_SA, cconfig));
            size_type sa_size = sa_buf.size()/2;
            this->width(bits::hi(text_size * (num + 1 + (str_vector.size()/num)))+1);
            this->resize(m_rank_marked(m_marked.size()));
            
            size_type sa_cnt = 0;
            for (size_type i=0; i < sa_size; ++i) {
                if (m_marked[i]) {
                    size_type strn  = sa_buf[2*i];
                    size_type pos   = sa_buf[2*i+1];
                    (*this)[sa_cnt++] = strn * text_size + pos;
                }
            }

            
            //this->resize(sa_cnt+1);
            m_sample_size = sa_cnt;
            m_strv = str_vector;

            //std::cout<< "Csa_sampling_strategy: width is " << (int)this->width() << std::endl;
            //std::cout<< "Csa_sampling_strategy: bit_size is " << this->bit_size() << std::endl;
            //std::cout<< "Csa_sampling_strategy: size is " << this->size() << std::endl;
            //std::cout<< "Csa_sampling_strategy: sample_cnt is " << m_sample_size << std::endl;
        }

        //! Copy constructor
        _saa_sampling(const _saa_sampling& st) : base_type(st) {
            m_marked = st.m_marked;
            m_rank_marked = st.m_rank_marked;
            m_rank_marked.set_vector(&m_marked);
        }

        //! Determine if index i is sampled or not
        inline bool is_sampled(size_type i) const {
            return m_marked[i];
        }
       
        //! Return the suffix array value for the sampled index i
        inline size_type saa_value(size_type i) const {
            size_type idx = m_rank_marked(i);
            return (*this)[idx];
        }
         
        inline size_type size() const {
            return m_sample_size;
        }
        
        //! Assignment operation
        _saa_sampling& operator=(const _saa_sampling& st) {
            if (this != &st) {
                base_type::operator=(st);
                m_marked = st.m_marked;
                m_rank_marked = st.m_rank_marked;
                m_rank_marked.set_vector(&m_marked);
            }
            return *this;
        }

        //! Swap operation
        void swap(_saa_sampling& st) {
            base_type::swap(st);
            m_marked.swap(st.m_marked);
            util::swap_support(m_rank_marked, st.m_rank_marked, &m_marked, &(st.m_marked));
            m_strv.swap(st.m_strv);
            std::swap(m_sample_size, st.m_sample_size);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0, temp_bytes = 0;
            
            written_bytes += base_type::serialize(out, child, "samples");
            written_bytes += m_marked.serialize(out, child, "marked");
            written_bytes += m_rank_marked.serialize(out, child, "rank_marked");
            written_bytes += m_strv.serialize(out, child, "strv");
            written_bytes += write_member(m_sample_size, out, child, "sample_size");
            
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            base_type::load(in);
            m_marked.load(in);
            m_rank_marked.load(in);
            m_rank_marked.set_vector(&m_marked);
            m_strv.load(in);
            read_member(m_sample_size, in);
        }
};

template<class t_bit_vec=bit_vector, class t_rank_sup=typename t_bit_vec::rank_1_type, uint8_t t_width=0>
//template<class t_bit_vec=rrr_vector<127>, class t_rank_sup=typename t_bit_vec::rank_1_type, uint8_t t_width=0>
class saa_sampling
{
    public:
        template<class t_csa> // template inner class which is used in CSAs to parametrize the
        class type            // sampling strategy class with the Sampling density of the CSA
        {
            public:
                typedef _saa_sampling<t_csa,
                        t_bit_vec,
                        t_rank_sup,
                        t_width
                        > sample_type;
        };
};

} // end namespace sdsl
