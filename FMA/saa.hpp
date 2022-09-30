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
/*! \file saa.hpp
    \brief saa.hpp contains an implementation suffix array of alignment.
    \author Seunghwan Min
    \reference
    [1] J. Na, H. Park, S. Lee, M. Hong, T. Lecroq, L. Mouchard, and K. Park, Suffix array of alignment: A practical index for similar data, In International Symposium on String Processing and Information Retrieval, 8214:243-254, Springer, Cham, 2013.
    \date 2017.11
*/

#include <iostream>
#include <fstream>
#include <mutex>
#include <thread>
#include <algorithm>
#include <tuple>
#include <set>

#include <sdsl/construct.hpp>
#include <sdsl/bit_vectors.hpp>
#include "gapinfo.hpp"
#include "io.hpp"
#include "custom_int_vector.hpp"

#pragma once

const int ALPHABETSIZE = 5;
const int TYPENUM = 2;
const int OCCID = 0;
const int EXTRAID = 1;
const char BASE[ALPHABETSIZE] = {'A', 'C', 'G', 'N', 'T'};

using namespace sdsl;

struct SAAConfig
{
    int numThread;
    int samplingDensity;

    SAAConfig(int nt = 8, int sd = 512) 
        : numThread(nt), samplingDensity(sd) {}
};

class SAA
{
private:
    // Sequence info
    int numSeq;
    std::string refSeq;
    std::vector<std::vector<std::string>> noncommons, ancgs;
    std::vector<int> ancgLength;
    std::vector<int> ancgLengthSum;

    // Configuration
    SAAConfig config;
    std::string seqDir;
    std::string outDir;

    std::mutex mtx;

    // Gap info
    GapInfo gapInfo;
    int alignLength;

    // Nonshared info
    int NSSize;
    std::vector<int> prevNSPos, NSPos;
    std::vector<int> prevrightNS, rightNS;
    bit_vector previsNS, isNS;
    bit_vector alphaPlus;
    std::vector<int> alphaLen, gammaLen;

    // SAA
    custom_int_vector gsa;
    custom_int_vector saa;
    std::vector<int> saaSeq;

    long long skip;

    std::map<bit_vector, int> vectorMap;
    int vectorID;
    bit_vector strv;

    int_vector<> C;
    bit_vector bv[TYPENUM][ALPHABETSIZE];
    bit_vector sampled;


    inline int ConvertToIndex(char ch)
    {
        switch(ch)
        {
            case '#':
                return 1;
            case 'A':
            case 'a':
                return 2;
            case 'C':
            case 'c':
                return 3;
            case 'G':
            case 'g':
                return 4;
            case 'N':
            case 'n':
                return 5;
            case 'T':
            case 't':
                return 6;
            default:
                return 0;
        }
    }

    void LoadNonCommon();
    void StoreData();
    void ClearData();

    void FindAlphaPlus();
    void FindAlphaPlusInner(int tid);
    void FindGammaPlus();
    void FindGammaPlusInner(int tid);
    void UpdateNonShared();
    void BuildANCG();
    void ConstructGSA();

    std::pair<int, int> GSAToSeqPos(long long idx);
    bool FromGamma(long long idx, int alignPos);
    int SeqComp(int s1, int p1, int s2, int p2, int compLen);
    std::string MakeSeq(int s, int p, int len);
    char AlignPosToChar(int seq, int alignPos);

    void ConstructSAA();
    int FindBindNumber(long long idx, int n);
    int GetSequenceID(long long idx, int numBinds, int seq_i);

    void CalculateCount();
    void MakeOcc();


public:
    SAA(int n, SAAConfig c, const std::string& seqDir_, const std::string& outDir_)
        : numSeq(n), config(c), seqDir(seqDir_), outDir(outDir_)
    {
        gapInfo = GapInfo(numSeq);
        gapInfo.LoadGapInfo(seqDir + "gapinfo");
        alignLength = gapInfo.GetAlignLength();

        rightNS = std::vector<int>(alignLength);
        isNS = bit_vector(alignLength, 0);
        alphaPlus = bit_vector(alignLength, 0);

        noncommons = std::vector<std::vector<std::string>>(numSeq);

        vectorID = numSeq * 2 + 1;
        skip = 0;
    }

    ~SAA() {}

    void BuildSAA();
};
