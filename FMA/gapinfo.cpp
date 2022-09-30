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

#include "gapinfo.hpp"

void GapInfo::Insert(int n, int p, int d)
{
    SeqGapInfos[n].Insert(p, d);
}

int GapInfo::ConvertToAlign(int n, int x)
{
    int lo = 0, hi = GetSize(n) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midx = SeqGapInfos[n].pos[mid];

        if(x < midx) hi = mid - 1;
        else lo = mid;
    }

    //int x1 = SeqGapInfos[n].pos[lo];
    int dx1 = SeqGapInfos[n].delta[lo];

    // Reference
    if(n == 0) return x + dx1;


    // Check insertion
    int x2, dx2;
    bool isInsertion = false;
    if(lo != GetSize(n) - 1)
    {
        x2 = SeqGapInfos[n].pos[lo + 1];
        dx2 = SeqGapInfos[n].delta[lo + 1];

        if(x2 + (dx2 - dx1) <= x) isInsertion = true;
    }

    int z = (isInsertion) ? x2 + dx2 : x + dx1;
    
    lo = 0, hi = GetSize(0) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midz = SeqGapInfos[0].pos[mid];

        if(z < midz) hi = mid - 1;
        else lo = mid;
    }
    
    //int z1 = SeqGapInfos[0].pos[lo];
    int dz1 = SeqGapInfos[0].delta[lo];

    int y = z + dz1;
    if(isInsertion) y -= (x2 - x);
    

    return y;
}

int GapInfo::ConvertToAlign(int n, int x, int p, int d)
{
    int lo = 0, hi = GetSize(n) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midx = SeqGapInfos[n].pos[mid];

        if(x < midx) hi = mid - 1;
        else lo = mid;
    }

    while(lo != GetSize(n) - 1 && 
        SeqGapInfos[n].pos[lo] == SeqGapInfos[n].pos[lo + 1])
        ++lo;

    //int x1 = SeqGapInfos[n].pos[lo];
    int dx1 = SeqGapInfos[n].delta[lo];
    if(n == 0 && lo > p) dx1 += d;

    // Reference
    if(n == 0) return x + dx1; 


    // Check insertion
    int x2, dx2;
    bool isInsertion = false;
    if(lo != GetSize(n) - 1)
    {
        x2 = SeqGapInfos[n].pos[lo + 1];
        dx2 = SeqGapInfos[n].delta[lo + 1];

        if(x2 + (dx2 - dx1) <= x) isInsertion = true;
    }

    int z = (isInsertion) ? x2 + dx2 : x + dx1;
    
    lo = 0, hi = GetSize(0) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midz = SeqGapInfos[0].pos[mid];

        if(z < midz) hi = mid - 1;
        else lo = mid;
    }
    
    while(lo != GetSize(0) - 1 && 
        SeqGapInfos[0].pos[lo] == SeqGapInfos[0].pos[lo + 1])
        ++lo;

    //int z1 = SeqGapInfos[0].pos[lo];
    int dz1 = SeqGapInfos[0].delta[lo];
    if(lo > p) dz1 += d;

    int y = z + dz1;
    if(isInsertion) y -= (x2 - x);
    

    return y;
}

int GapInfo::ConvertToInd(int n, int y)
{
    int lo = 0, hi = GetSize(0) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midz = SeqGapInfos[0].pos[mid];
        int middz = SeqGapInfos[0].delta[mid];

        if(y < midz + middz) hi = mid - 1;
        else lo = mid;
    }

    //int z1 = SeqGapInfos[0].pos[lo];
    int dz1 = SeqGapInfos[0].delta[lo];

    // Check insertion
    int z2, dz2;
    bool isInsertion = false;
    if(lo != GetSize(0) - 1)
    {
        z2 = SeqGapInfos[0].pos[lo + 1];
        dz2 = SeqGapInfos[0].delta[lo + 1];

        if(z2 + dz1 <= y) isInsertion = true;
    }

    int z = (isInsertion) ? z2 : y - dz1;
    if(n == 0) return z;

    lo = 0, hi = GetSize(n) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midx = SeqGapInfos[n].pos[mid];
        int middx = SeqGapInfos[n].delta[mid];

        if(z < midx + middx) hi = mid - 1;
        else lo = mid;
    }

    //int x1 = SeqGapInfos[n].pos[lo];
    int dx1 = SeqGapInfos[n].delta[lo];

    int x = z - dx1;
    if(isInsertion)
    {
        int dist = z2 + dz2 - y;
        // Check invaild position(other sequences' insertion)
        if(lo != 0)
        {
            int lo2 = 0, hi2 = GetSize(n) - 1;
            while(lo2 < hi2)
            {
                int mid = (lo2 + hi2 + 1) / 2;
                int midx = SeqGapInfos[n].pos[mid];
                int middx = SeqGapInfos[n].delta[mid];

                if(z - 1 < midx + middx) hi2 = mid - 1;
                else lo2 = mid;
            }

            int l = x - ((z -1) - SeqGapInfos[n].delta[lo2]) - 1;
            dist = std::min(dist, l);
        }
        x -= dist;
    }

    // Check invaild position(deletion)
    int x2, dx2;
    bool isDeletion = false;
    if(lo != GetSize(n) - 1)
    {
        x2 = SeqGapInfos[n].pos[lo + 1];
        dx2 = SeqGapInfos[n].delta[lo + 1];

        if(x2 + dx1 <= z) isDeletion = true;
    }
    if(isDeletion) x = x2;

    return x;
}

int GapInfo::ConvertToInd(int n, int y, int p, int d)
{
    int lo = 0, hi = GetSize(0) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midz = SeqGapInfos[0].pos[mid];
        int middz = SeqGapInfos[0].delta[mid];
        if(mid > p) middz += d;

        if(y < midz + middz) hi = mid - 1;
        else lo = mid;
    }

    while(lo != GetSize(0) - 1 && 
        SeqGapInfos[0].pos[lo] == SeqGapInfos[0].pos[lo + 1])
        ++lo;

    //int z1 = SeqGapInfos[0].pos[lo];
    int dz1 = SeqGapInfos[0].delta[lo];
    if(lo > p) dz1 += d;

    if(n == 0) return y - dz1;

    // Check insertion
    int z2, dz2;
    bool isInsertion = false;
    if(lo != GetSize(0) - 1)
    {
        z2 = SeqGapInfos[0].pos[lo + 1];
        dz2 = SeqGapInfos[0].delta[lo + 1];
        if(lo + 1 > p) dz2 += d;

        if(z2 + dz1 <= y) isInsertion = true;
    }

    int z = (isInsertion) ? z2 : y - dz1;

    lo = 0, hi = GetSize(n) - 1;
    while(lo < hi)
    {
        int mid = (lo + hi + 1) / 2;
        int midx = SeqGapInfos[n].pos[mid];
        int middx = SeqGapInfos[n].delta[mid];

        if(z < midx + middx) hi = mid - 1;
        else lo = mid;
    }

    while(lo != GetSize(n) - 1 && 
        SeqGapInfos[n].pos[lo] == SeqGapInfos[n].pos[lo + 1])
        ++lo;

    //int x1 = SeqGapInfos[n].pos[lo];
    int dx1 = SeqGapInfos[n].delta[lo];

    int x = z - dx1;
    if(isInsertion) x -= (z2 + dz2 - y);


    return x;
}

void GapInfo::RemoveDuplication()
{
    for(int i = 0; i < numSeq; ++i)
    {
        std::vector<int>& pos = SeqGapInfos[i].pos;
        pos.erase(std::unique(pos.begin(), pos.end()), pos.end());
        std::vector<int>& delta = SeqGapInfos[i].delta;
        delta.erase(std::unique(delta.begin(), delta.end()), delta.end());
    }
}

void GapInfo::LoadGapInfo(std::string fname)
{
    std::ifstream in(fname, std::ifstream::in);

    // File open error
    if(!in.is_open()) return;

    for(auto& SeqGapInfo : SeqGapInfos)
    {
        SeqGapInfo.sumDelta = 0;
        SeqGapInfo.pos.clear();
        SeqGapInfo.delta.clear();
    }
    
    in >> alignLength;
    int n, p, d;
    for(int i = 0; i < numSeq; ++i)
    {
        in >> n;
        int delta = 0;
        for(int j = 0; j < n; ++j)
        {
            in >> p >> d;
            SeqGapInfos[i].Insert(p, d - delta);
            delta = d;
        }
    }

    in.close();
}

void GapInfo::StoreGapInfo(std::string fname)
{
    std::ofstream out(fname, std::ofstream::out);

    // File open error
    if(!out.is_open()) return;

    out << alignLength << std::endl;
    for(int i = 0; i < numSeq; ++i)
    {
        out << GetSize(i) << std::endl;
        for(int j = 0; j < GetSize(i); ++j)
        {
            out << SeqGapInfos[i].pos[j] << " " << SeqGapInfos[i].delta[j] << std::endl;
        }
    }

    out.close();
}

void GapInfo::StoreIVGapInfo(std::string fname)
{
    int totalSize = 1;
    for(int i = 0; i < numSeq; ++i)
        totalSize += GetSize(i) * 2 + 1;

    int_vector<> info(totalSize, 0, bits::hi(alignLength) + 2);
    int negativebit = bits::hi(alignLength) + 1;
    int cnt = 0;
    info[cnt++] = alignLength;

    for(int i = 0; i < numSeq; ++i)
    {
        int n = GetSize(i);
        info[cnt++] = n;
        for(int j = 0; j < n; ++j)
        {
            info[cnt++] = SeqGapInfos[i].pos[j];
            if(SeqGapInfos[i].delta[j] < 0)
                info[cnt++] = (-SeqGapInfos[i].delta[j]) | (1 << negativebit);
            else
                info[cnt++] = SeqGapInfos[i].delta[j];
        }
    }

    store_to_file(info, fname);
}



