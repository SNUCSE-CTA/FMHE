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
/*! \file gapinfo.hpp
    \brief gapinfo.hpp contains an implementation of class which stores information of gaps.
    \author Seunghwan Min
    \reference
    [1] J. Na, H. Kim, S. Min, H. Park, T. Lecroq, M. Le'onard, L. Mouchard, and K. Park, FM-index of alignment with gaps, TCS 710:148-157, 2018
    \date 2017.11
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sdsl/int_vector.hpp>

#pragma once

using namespace sdsl;

class GapInfo
{
private:
    struct SeqGapInfo
    {
        std::vector<int> pos;
        std::vector<int> delta;

        int sumDelta;

        SeqGapInfo() : sumDelta(0) {}

        void Insert(int p, int d)
        {
            sumDelta += d;
            pos.push_back(p);
            delta.push_back(sumDelta);
        }
    };

    std::vector<SeqGapInfo> SeqGapInfos;
    int numSeq;
    int alignLength;

public:
    GapInfo() {}
    GapInfo(int n) : SeqGapInfos(std::vector<SeqGapInfo>(n)), numSeq(n)
    {
        for(int i = 0; i < n; ++i)
            SeqGapInfos[i].Insert(-1, 0);
    }
    GapInfo(int n, int_vector<> info) : SeqGapInfos(std::vector<SeqGapInfo>(n)), numSeq(n)
    {
        int cnt = 0;
        alignLength = info[cnt++];
        int negativebit = bits::hi(alignLength) + 1;

        for(int i = 0; i < numSeq; ++i)
        {
            int infosize = info[cnt++];
            int delta = 0;
            for(int j = 0; j < infosize; ++j)
            {
                int p = info[cnt++];
                int d;
                if(info[cnt] & (1 << negativebit)) d = -(info[cnt++] & (~(1 << negativebit)));
                else d = info[cnt++];

                SeqGapInfos[i].Insert(p, d - delta);
                delta = d;
            }
        }
    }
    ~GapInfo() {}

    inline int GetSize() const {return numSeq;}
    inline int GetSize(int n) const {return SeqGapInfos[n].pos.size();}
    inline int GetAlignLength() const {return alignLength;}
    inline int GetPosition(int n, int i) const {return SeqGapInfos[n].pos[i];}
    inline int GetDelta(int n, int i) const {return SeqGapInfos[n].delta[i];}

    void SetAlignLength(const int al) {alignLength = al;}
    void SetPosition(int n, int i, int p) {SeqGapInfos[n].pos[i] = p;}
    void SetDelta(int n, int i, int d) {SeqGapInfos[n].delta[i] = d;}

    void Insert(int n, int p, int d);
    int ConvertToAlign(int n, int x);
    int ConvertToAlign(int n, int x, int p, int d);
    int ConvertToInd(int n, int y);
    int ConvertToInd(int n, int y, int p, int d);

    void RemoveDuplication();

    void LoadGapInfo(std::string fname);
    void StoreGapInfo(std::string fname);
    void StoreIVGapInfo(std::string fname);

    void Print()
    {
        for(int i = 0; i < numSeq; ++i)
        {
            std::cout << i << " : ";
            for(int j = 0; j < GetSize(i); ++j) 
                std::cout << "(" << SeqGapInfos[i].pos[j] << ", " << SeqGapInfos[i].delta[j] << ") ";
            std::cout << "\n";
        }
    }
};