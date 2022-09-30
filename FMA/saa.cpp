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

#include "saa.hpp"

void SAA::BuildSAA()
{
    LoadNonCommon();
    std::cout << "LoadNonCommon Finished" << std::endl;

    FindAlphaPlus();
    std::cout << "FindAlphaPlus Finished" << std::endl;
    UpdateNonShared();
    std::cout << "UpdateNonShared Finished" << std::endl;
    FindGammaPlus();
    std::cout << "FindGammaPlus Finished" << std::endl;
    BuildANCG();
    std::cout << "BuildANCG Finished" << std::endl;

    ConstructGSA();
    std::cout << "ConstructGSA Finished" << std::endl;
    ConstructSAA();
    std::cout << "ConstructSAA Finished" << std::endl;

    CalculateCount();
    std::cout << "CalculateCount Finished" << std::endl;
    MakeOcc();
    std::cout << "MakeOcc Finished" << std::endl;

    StoreData();
    std::cout << "StoreData Finished" << std::endl;
    ClearData();
}

void SAA::LoadNonCommon()
{
    std::ifstream in(seqDir + "noncommon", std::ifstream::in);

    // File open error
    if(!in.is_open()) return;

    int start, end, last = 0, cnt = 0;
    while(in >> start >> end)
    {
        NSPos.push_back(start);
        NSPos.push_back(end);

        for(int i = start; i <= end; ++i) isNS[i] = 1;
        for(int i = last; i <= end; ++i) rightNS[i] = cnt;
        last = end + 1;
        ++cnt;
    }
    for(int i = last; i < alignLength; ++i) rightNS[i] = -1;

    NSSize = NSPos.size() / 2;

    in.close();
}

void SAA::StoreData()
{
    std::string fname = outDir + "strv" + std::to_string(numSeq) + ".vec";
    rrr_vector<127> rrr(strv);
    store_to_file(rrr, fname);

    fname = outDir + "nonshared" + std::to_string(numSeq) + ".vec";
    StoreNonShared(fname, NSPos, alignLength, numSeq);

    gapInfo.StoreGapInfo(outDir + "gapinfo");
    gapInfo.StoreIVGapInfo(outDir + "gapinfo.vec");

    fname = outDir + "cnt" + std::to_string(numSeq) + ".vec";
    store_to_file(C, fname);

    fname = outDir + std::to_string(config.samplingDensity) + "sampled" + std::to_string(numSeq) + ".vec";
    store_to_file(sampled, fname);

    for(int i = 0; i < TYPENUM; ++i)
    {
        for(int j = 0; j < ALPHABETSIZE; ++j)
        {
            rrr_vector<127> b(bv[i][j]);
            std::string temp = std::string(1, BASE[j]);
            std::string type = (i == 0) ? "occ" : "ext";

            fname = outDir + type + std::to_string(numSeq) + temp[0] + ".vec";
            store_to_file(b, fname);
        }
    }
}

void SAA::ClearData()
{
    util::clear(previsNS);
    util::clear(isNS);
    util::clear(strv);
    util::clear(C);
    for(int i = 0; i < TYPENUM; ++i)
        for(int j = 0; j < ALPHABETSIZE; ++j)
            util::clear(bv[i][j]);
    util::clear(sampled);
}

void SAA::FindAlphaPlus()
{ 
    alphaLen = std::vector<int>(NSSize, 0);

    std::vector<std::thread> threads;

    for(int i = 0; i < config.numThread - 1; ++i)
    {
        threads.push_back(std::thread([=] { FindAlphaPlusInner(i); }));
    }
    FindAlphaPlusInner(config.numThread - 1);

    for(int i = 0; i < config.numThread - 1; ++i)
    {
        threads[i].join();
    }

    if(NSPos[0] < alphaLen[0])
    {
        alphaLen[0] = NSPos[0];
    }


    for(int i = 0; i < NSSize; ++i)
    {
        int NSStart = NSPos[i * 2];
        for(int j = NSStart - 1; j > NSStart - 1 - alphaLen[i]; --j)
            alphaPlus[j] = 1;
    }
}

void SAA::FindAlphaPlusInner(int tid)
{
    int start = (numSeq / config.numThread) * tid + std::min(numSeq % config.numThread, tid);
    int end = (numSeq / config.numThread) * (tid+1) + std::min(numSeq % config.numThread, tid+1);

    for(int sn = start; sn < end; ++sn)
    {
        // Load sequence
        char fname[256];
        snprintf(fname, 256, "%sS%04d", seqDir.c_str(), sn);
        int length = GetLength(fname) + 1;        
        char *reverseSeq = new char[length]; 
        LoadSequence(fname, reverseSeq, length - 1);
        std::reverse(reverseSeq, reverseSeq + length - 1);
        reverseSeq[length - 1] = 0;
        for(int i = 0; i < length; ++i) reverseSeq[i] = ConvertToIndex(reverseSeq[i]);

        // Construct suffix array
        int_vector<> reverseSA;
        reverseSA.width(bits::hi(length) + 1);
        reverseSA.resize(length);
        algorithm::calculate_sa((const unsigned char*)reverseSeq, length, reverseSA);
        
        for(int i = 1; i < reverseSA.size(); ++i)
        {
            int pos = length - reverseSA[i] - 2;
            int alignPos = gapInfo.ConvertToAlign(sn, pos);

            if(!isNS[alignPos])
            {
                int target = rightNS[alignPos];
                if(target != -1 && alignPos == NSPos[target * 2] - 1)
                {
                    int alphaLength = 0;
                    if(i != 1)
                    {
                        int cnt = 0;
                        while(reverseSeq[reverseSA[i] + cnt] == reverseSeq[reverseSA[i-1] + cnt]) ++cnt;

                        alphaLength = std::max(alphaLength, cnt);
                    }

                    if(i + 1 < reverseSA.size())
                    {
                        int cnt = 0;
                        while(reverseSeq[reverseSA[i] + cnt] == reverseSeq[reverseSA[i+1] + cnt]) ++cnt;

                        alphaLength = std::max(alphaLength, cnt);
                    }

                    mtx.lock();
                    alphaLen[target] = std::max(alphaLen[target], alphaLength + 1);
                    mtx.unlock();
                }
            }
        }

        util::clear(reverseSA);
        delete[] reverseSeq;
    }
}

void SAA::FindGammaPlus()
{
    gammaLen = std::vector<int>(NSSize, 0);

    std::vector<std::thread> threads;

    for(int i = 0; i < config.numThread - 1; ++i)
    {
        threads.push_back(std::thread([=] { FindGammaPlusInner(i); }));
    }
    FindGammaPlusInner(config.numThread - 1);

    for(int i = 0; i < config.numThread - 1; ++i)
    {
        threads[i].join();
    }
}

void SAA::FindGammaPlusInner(int tid)
{
    int start = (numSeq / config.numThread) * tid + std::min(numSeq % config.numThread, tid);
    int end = (numSeq / config.numThread) * (tid+1) + std::min(numSeq % config.numThread, tid+1);

    for(int sn = start; sn < end; ++sn)
    {
        // Load sequence
        char fname[256];
        snprintf(fname, 256, "%sS%04d", seqDir.c_str(), sn);
        int length = GetLength(fname) + 1;        
        char *seq = new char[length]; 
        LoadSequence(fname, seq, length - 1);
        seq[length - 1] = 0;
        for(int i = 0; i < length; ++i) seq[i] = ConvertToIndex(seq[i]);

        // Construct suffix array
        int_vector<> SA;
        SA.width(bits::hi(length) + 1);
        SA.resize(length);
        algorithm::calculate_sa((const unsigned char*)seq, length, SA);

        for(int i = 1; i < SA.size(); ++i)
        {
            int pos = SA[i];
            int alignPos = gapInfo.ConvertToAlign(sn, pos);

            if(!isNS[alignPos])
            {
                int leftTarget = rightNS[alignPos] - 1;
                if(leftTarget == -2) leftTarget = NSSize - 1;
                if(leftTarget != -1 && alignPos == NSPos[leftTarget * 2 + 1] + 1)
                {
                    int gammaLength = 0;
                    if(i != 1)
                    {
                        int pos_prev = SA[i-1];

                        int cnt = 0;
                        while(seq[pos + cnt] == seq[pos_prev + cnt]) ++cnt;

                        gammaLength = std::max(gammaLength, cnt);
                    }

                    if(i + 1 < SA.size())
                    {
                        int pos_next = SA[i+1];

                        int cnt = 0;
                        while(seq[pos + cnt] == seq[pos_next + cnt]) ++cnt;

                        gammaLength = std::max(gammaLength, cnt);
                    }

                    mtx.lock();
                    gammaLen[leftTarget] = std::max(gammaLen[leftTarget], gammaLength + 1);
                    mtx.unlock();
                }
            }
        }

        util::clear(SA);
        delete[] seq;
    }
}

void SAA::UpdateNonShared()
{
    // Combine NC & alphaPlus
    int start = -1, end = -1;
    int last = 0, curr = 0;
    for(int nc = 0; nc < NSSize; ++nc)
    {
        if(start == -1) 
        {
            start = NSPos[nc * 2];
            while((start - 1) >= 0 && alphaPlus[start - 1]) --start;
        }
        end = NSPos[nc * 2 + 1];

        if(end + 1 == alignLength || alphaPlus[end + 1] == 0)
        {
            if(nc != NSSize - 1)
                alphaLen[curr + 1] = alphaLen[nc + 1];
            NSPos[curr * 2] = start;
            NSPos[curr * 2 + 1] = end;

            last = end + 1;
            ++curr;
            start = -1;
        }
    }
    alphaLen.resize(curr);
    NSPos.resize(curr * 2);
    NSSize = curr;

    // Align right
    std::vector<int> previdx(numSeq, 0), curridx(numSeq, 0), indelSize(numSeq, 0);
    int refUpdatePos = 0, refUpdateVal = 0;
    last = curr = 0;
    for(int ns = 0; ns < NSSize; ++ns)
    {
        for(int i = last; i < NSPos[ns * 2]; ++i) 
        {
            rightNS[curr] = ns;
            isNS[curr++] = alphaPlus[i];
        }
        last = NSPos[ns * 2 + 1] + 1;

        NSPos[ns * 2] += refUpdateVal; NSPos[ns * 2 + 1] += refUpdateVal;
        int start = NSPos[ns * 2], end = NSPos[ns * 2 + 1];
        
        for(int i = 0; i < numSeq; ++i)
        {
            previdx[i] = curridx[i];
            while(true)
            {
                if(curridx[i] == gapInfo.GetSize(i) - 1) break;

                int p = gapInfo.ConvertToAlign(i, gapInfo.GetPosition(i, curridx[i] + 1), refUpdatePos, refUpdateVal);
                if(start <= p && p <= end + 1) ++curridx[i];
                else break;
            }

            indelSize[i] = -(gapInfo.GetDelta(i, curridx[i]) - gapInfo.GetDelta(i, previdx[i]));
            if(i == 0 && curridx[i] > refUpdatePos) indelSize[i] -= refUpdateVal;
        }
        
        
        int maxInsertion = *(std::max_element(indelSize.begin(), indelSize.end()));
        maxInsertion = std::max(maxInsertion, 0);

        for(int i = 1;; ++i)
        {
            if(i == numSeq) i = 0;

            int p = gapInfo.ConvertToInd(i, start - 1, refUpdatePos, refUpdateVal) + 1;
            if(indelSize[i] > 0) p += indelSize[i];
            if(indelSize[i] == 0) p = gapInfo.GetPosition(i, previdx[i]);
            if(i == 0 && maxInsertion == 0) p = gapInfo.GetPosition(i, previdx[i]);

            int d = gapInfo.GetDelta(i, previdx[i]);
            if(i == 0) d += maxInsertion;
            else d -= indelSize[i];

            for(int j = previdx[i] + 1; j <= curridx[i]; ++j)
            {
                gapInfo.SetPosition(i, j, p);
                gapInfo.SetDelta(i, j, d);
            }

            if(i == 0)
            {
                refUpdateVal += (indelSize[0] + maxInsertion);
                refUpdatePos = curridx[0];
                break;
            }
        }

        NSPos[ns * 2 + 1] += (indelSize[0] + maxInsertion);
        for(int i = NSPos[ns * 2]; i <= NSPos[ns * 2 + 1]; ++i)
        {
            rightNS[curr] = ns;
            isNS[curr++] = 1;
        }
    }


    alignLength += refUpdateVal;
    for(int i = curr; i < alignLength; ++i)
    {
        rightNS[curr] = -1;
        isNS[curr++] = 0;
    }
    rightNS.resize(alignLength);
    isNS.resize(alignLength);

    util::clear(alphaPlus);

    // Remove duplication
    gapInfo.RemoveDuplication();
    gapInfo.SetAlignLength(alignLength);
}

void SAA::BuildANCG()
{
    ancgLength = std::vector<int>(NSSize, 0);
    ancgLengthSum = std::vector<int>(NSSize, 0);

    LoadSequence(seqDir + "S0000", refSeq);
    std::string indSeq;
    char fname[256];

    std::string noncommon;
    for(int i = 1; i < numSeq; ++i)
    {
        noncommons[i].reserve(NSSize);
        snprintf(fname, 256, "%sS%04d", seqDir.c_str(), i);
        LoadSequence(fname, indSeq);
        for(int j = 0; j < NSSize; ++j)
        {
            int NSStart = NSPos[2 * j];
            int NSEnd = NSPos[2 * j + 1];
            int NSStartInd = gapInfo.ConvertToInd(i, NSStart - 1) + 1;
            int NSEndInd = gapInfo.ConvertToInd(i, NSEnd + 1) - 1;

            noncommon = indSeq.substr(NSStartInd + alphaLen[j], NSEndInd - NSStartInd + 1 - alphaLen[j]);
            noncommons[i].push_back(noncommon);
            ancgLength[j] = std::max(ancgLength[j], (int)noncommon.length() + alphaLen[j] + gammaLen[j]);
        }
    }

    int sum = 0;
    for(int i = 0; i < NSSize; ++i)
    {
        sum += ancgLength[i] + 1;
        ancgLengthSum[i] = sum;
    }
}

void SAA::ConstructGSA()
{
    long long totalLength = refSeq.length() + 1;
    totalLength += (long long)ancgLengthSum[NSSize - 1] * (numSeq - 1);
    totalLength += 1;
    
    int_vector<8> text(totalLength, 0);
    long long cnt = 0;

    for(auto c : refSeq) 
    {
        text[cnt] = ConvertToIndex(c) + 1;
        if(text[cnt++] < 3) ++skip;
    }
    text[cnt++] = 1;
    ++skip;
    
    for(int idx = 1 ; idx < numSeq; ++idx)
    {
        for(int i = 0; i < NSSize; ++i)
        {
            int NSStart = NSPos[2 * i];
            int NSEnd = NSPos[2 * i + 1];
            int NSStartRef = gapInfo.ConvertToInd(0, NSStart - 1) + 1;
            int NSEndRef = gapInfo.ConvertToInd(0, NSEnd + 1) - 1;

            for(int j = NSStartRef; j < NSStartRef + alphaLen[i]; ++j)
            {
                text[cnt] = ConvertToIndex(refSeq[j]) + 1;
                if(text[cnt++] < 3) ++skip;
            }

            for(auto c : noncommons[idx][i])
            {
                text[cnt] = ConvertToIndex(c) + 1;
                if(text[cnt++] < 3) ++skip;
            }

            for(int j = NSEndRef + 1; j < NSEndRef + 1 + gammaLen[i]; ++j)
            {
                text[cnt] = ConvertToIndex(refSeq[j]) + 1;
                if(text[cnt++] < 3) ++skip;
            }

            for(int j = 0; j < ancgLength[i] - (noncommons[idx][i].length() + alphaLen[i] + gammaLen[i]) + 1; ++j) 
            {
                text[cnt++] = 1;
                ++skip;
            }
        }
    }
    text[cnt++] = 0;
    ++skip;

    cache_config config(false, outDir, std::to_string(numSeq));
    _construct_sa_se<int_vector<8>>(text, cache_file_name(conf::KEY_SA, config), 256, 0);
    util::clear(text);
    //store_to_cache(text, conf::KEY_TEXT, config);
    //util::clear(text);

    //construct_sa_se(config);
    //remove(cache_file_name(conf::KEY_TEXT, config));

    std::string fname = outDir + "sa_" + std::to_string(numSeq) + ".sdsl";
    gsa.read_file(fname);
}

std::pair<int, int> SAA::GSAToSeqPos(long long idx)
{
    int seq, pos;
    unsigned long long gsaVal = gsa[idx];
    
    if(gsaVal < refSeq.length() + 1)
    {
        seq = 0;
        pos = gsaVal;
    }
    else
    {
        gsaVal -= (refSeq.length() + 1);
        seq = (gsaVal / ancgLengthSum[NSSize - 1]) + 1;
        int ancgPos = gsaVal % ancgLengthSum[NSSize - 1];

        int lo = 0, hi = NSSize - 1;
        while(lo < hi)
        {
            int mid = (lo + hi) / 2;
            if(ancgLengthSum[mid] > ancgPos) hi = mid;
            else lo = mid + 1;
        }

        pos = gapInfo.ConvertToInd(seq, NSPos[2 * lo] - 1) + 1;
        pos += ancgPos;
        if(lo != 0) pos -= ancgLengthSum[lo - 1];
    }

    return std::make_pair(seq, pos);
}

bool SAA::FromGamma(long long idx, int alignPos)
{
    unsigned long long gsaVal = gsa[idx];
    gsaVal -= (refSeq.length() + 1);

    int ancgPos = gsaVal % ancgLengthSum[NSSize - 1];

    int lo = 0, hi = NSSize - 1;
    while(lo < hi)
    {
        int mid = (lo + hi) / 2;
        if(ancgLengthSum[mid] > ancgPos) hi = mid;
        else lo = mid + 1;
    }

    if(rightNS[alignPos] == lo) return false;
    else return true;
}

int SAA::SeqComp(int s1, int p1, int s2, int p2, int compLen)
{
    std::string str1 = MakeSeq(s1, p1, compLen);
    std::string str2 = MakeSeq(s2, p2, compLen);

    return strncmp(str1.c_str(), str2.c_str(), compLen);
}

std::string SAA::MakeSeq(int s, int p, int len)
{
    if(s == 0)
    {
        return refSeq.substr(p, len);
    }
    else
    {
        std::string ret;
        int remainLen = len;
        while(remainLen > 0)
        {
            int pos = gapInfo.ConvertToAlign(s, p);
            int NS = prevrightNS[pos];

            if(isNS[pos])
            {
                int NSStart = prevNSPos[NS * 2];
                int NSEnd = prevNSPos[NS * 2 + 1];
                int NSStartInd = gapInfo.ConvertToInd(s, NSStart - 1) + 1;
                int NSEndInd = gapInfo.ConvertToInd(s, NSEnd + 1) - 1;

                if(NSStartInd + alphaLen[NS] > p)
                {
                    int l = std::min(remainLen, NSStartInd + alphaLen[NS] - p);
                    int NSStartRef = gapInfo.ConvertToInd(0, NSStart - 1) + 1;
                    ret += refSeq.substr(NSStartRef + p - NSStartInd, l);

                    remainLen -= l;
                    p += l;
                }


                int l = std::min(remainLen, NSEndInd - p + 1);
                ret += noncommons[s][NS].substr(p - NSStartInd - alphaLen[NS], l);
                
                remainLen -= l;
                p += l;
            }
            else
            {
                int pRef = gapInfo.ConvertToInd(0, pos);
                int SEnd = prevNSPos[NS * 2] - 1;
                int SEndRef = gapInfo.ConvertToInd(0, SEnd);

                int l = std::min(remainLen, SEndRef - pRef + 1);
                ret += refSeq.substr(pRef, l);

                remainLen -= l;
                p += l;
            }
        }
        return ret;
    }
}

char SAA::AlignPosToChar(int seq, int alignPos)
{
    if(seq == 0)
    {
        int pos = gapInfo.ConvertToInd(seq, alignPos);
        return refSeq[pos];
    }
    else
    {
        if(previsNS[alignPos])
        {
            int pos = gapInfo.ConvertToInd(seq, alignPos);
            int NS = prevrightNS[alignPos];
            int NSStartInd = gapInfo.ConvertToInd(seq, prevNSPos[2 * NS] - 1) + 1;

            
            if(NSStartInd + alphaLen[NS] > pos)
            {
                int NSStartRef = gapInfo.ConvertToInd(0, prevNSPos[2 * NS] - 1) + 1;
                return refSeq[NSStartRef + pos - NSStartInd];
            }
            else
            {
                return noncommons[seq][NS][pos-NSStartInd-alphaLen[NS]];
            }
        }
        else
        {
            int pos = gapInfo.ConvertToInd(0, alignPos);
            return refSeq[pos];
        }
    }
}

void SAA::ConstructSAA()
{
    // Save non-shared
    prevNSPos = NSPos;
    previsNS.resize(isNS.size());
    for(int i = 0; i < isNS.size(); ++i) previsNS[i] = isNS[i];
    prevrightNS.resize(rightNS.size());
    for(int i = 0; i < rightNS.size(); ++i) prevrightNS[i] = rightNS[i];

    std::string fname = outDir + "finSaa" + std::to_string(numSeq) + ".vec";
    std::fstream fs(fname, std::ios::out | std::ios::binary);

    uint64_t m_size = 0;
    uint8_t m_width = 32;
    fs.write((char*)&m_size, sizeof(m_size));
    fs.write((char*)&m_width, sizeof(m_width));

    int cmp, curr = 0;
    int seqID, numBinds;

    int tmp = 0;
    fs.write((char*)&tmp, sizeof(tmp));
    tmp = alignLength - 1;
    fs.write((char*)&tmp, sizeof(tmp));
    m_size += m_width * 2;

    for(long long idx = skip; idx < gsa.size(); )
    {
        int seq, pos;
        std::tie(seq, pos) = GSAToSeqPos(idx);
        int alignPos = gapInfo.ConvertToAlign(seq, pos);

        // Remove index start from gamma region
        int leftTarget = rightNS[alignPos] - 1;
        if(leftTarget == -2) leftTarget = NSSize - 1;

        if(leftTarget != -1)
        {
            int leftNSEndInd = gapInfo.ConvertToInd(seq, NSPos[2 * leftTarget + 1] + 1) - 1;
            if(pos <= leftNSEndInd + gammaLen[leftTarget])
            {
                if(isNS[alignPos])
                {
                    if(seq != 0 && FromGamma(idx, alignPos))
                    {
                        idx += 1;
                        continue;
                    }
                }
                else
                {
                    if(seq == 0)
                    {
                        seqID = 0;
                        fs.write((char*)&seqID, sizeof(seqID));
                        fs.write((char*)&alignPos, sizeof(alignPos));
                        m_size += m_width * 2;
                    }
                    idx += 1;
                    continue;
                }
            }
        }

        // Reamin < numSeq
        if(idx + numSeq > gsa.size())
        {
            int remain = gsa.size() - idx;
            if(!isNS[alignPos])
            {
                numBinds = 1;
                seqID = 0;
            }
            else
            {
                if(remain == 1)
                {
                    numBinds = 1;
                    seqID = seq + 1;
                }
                else
                {
                    numBinds = FindBindNumber(idx, remain);
                    seqID = GetSequenceID(idx, numBinds, seq);
                }
            }

            fs.write((char*)&seqID, sizeof(seqID));
            fs.write((char*)&alignPos, sizeof(alignPos));
            m_size += m_width * 2;
            idx += numBinds;
            continue;
        }

        if(!isNS[alignPos])
        {
            numBinds = 1;
            seqID = 0;
        }
        else
        {
            numBinds = FindBindNumber(idx, numSeq);
            seqID = GetSequenceID(idx, numBinds, seq);
        }

        fs.write((char*)&seqID, sizeof(seqID));
        fs.write((char*)&alignPos, sizeof(alignPos));
        m_size += m_width * 2;
        idx += numBinds;
    }

    fs.seekp(0, std::ios::beg);
    fs.write((char*)&m_size, sizeof(m_size));
    fs.close();

    saa.read_file(fname);

    // Make strv
    long long num = vectorMap.size() + numSeq * 2 + 1;
    strv = bit_vector(num * numSeq, 0);

    for(int i = 0; i < numSeq; ++i) 
        strv[i] = 1;
    for(int i = 0; i < numSeq; ++i)
    {
        strv[(i+1) * numSeq + i] = 1;
        for(int j = 0; j < numSeq; ++j)
            if(i != j)
                strv[(i+numSeq+1) * numSeq +j] = 1;
    }

    for(auto it = vectorMap.begin(); it != vectorMap.end(); ++it)
    {
        bit_vector vec = it->first;
        int id = it->second;
        for(int i = 0; i < numSeq; ++i)
        {
            long long strvidx = (long long)id * numSeq + i;
            if(vec[i]) strv[strvidx] = 1;
        }
        util::clear(vec);
    }

    // Consider special char of sequence
    // alignLength -= 1;
    gapInfo.SetAlignLength(alignLength);

    // Update non-shared
    int cnt = 0;
    NSPos.clear();
    for(int i = 1; i < alignLength; ++i)
    {
        rightNS[i] = cnt;
        // Start of non-shared
        if(isNS[i] && !isNS[i-1])
        {
            NSPos.push_back(i);
        }
        // End of non-shared
        else if(!isNS[i] && isNS[i-1])
        {
            NSPos.push_back(i - 1);
            ++cnt;
        }
    }
    if(NSPos.size() % 2 == 1) NSPos.push_back(alignLength - 1);
    NSSize = NSPos.size() / 2;
}

int SAA::FindBindNumber(long long idx, int n)
{
    int cmp, seq_i, pos_i;
    std::tie(seq_i, pos_i) = GSAToSeqPos(idx);
    int alignPos_i = gapInfo.ConvertToAlign(seq_i, pos_i);

    int targetNS = rightNS[alignPos_i];
    int compareLength = NSPos[targetNS * 2 + 1] - alignPos_i + 1;

    // Compare idx and idx+1
    int seq_j, pos_j;
    std::tie(seq_j, pos_j) = GSAToSeqPos(idx + 1);
    int alignPos_j = gapInfo.ConvertToAlign(seq_j, pos_j);

    if(alignPos_i != alignPos_j) return 1;
    cmp = SeqComp(seq_i, pos_i, seq_j, pos_j, compareLength);
    if(cmp != 0) return 1;

    // Compare idx and idx+numSeq-1
    std::tie(seq_j, pos_j) = GSAToSeqPos(idx + n - 1);
    alignPos_j = gapInfo.ConvertToAlign(seq_j, pos_j);


    if(alignPos_i == alignPos_j)
    {
        cmp = SeqComp(seq_i, pos_i, seq_j, pos_j, compareLength);
        if(cmp == 0) return n;
    }
    
    // Compare idx and idx+numSeq-2
    std::tie(seq_j, pos_j) = GSAToSeqPos(idx + n - 2);
    alignPos_j = gapInfo.ConvertToAlign(seq_j, pos_j);

    if(alignPos_i == alignPos_j)
    {
        cmp = SeqComp(seq_i, pos_i, seq_j, pos_j, compareLength);
        if(cmp == 0) return n - 1;
    }

    int lo = 0, hi= n - 1;
    while(lo < hi)
    {
        int mid = (lo + hi) / 2;
        std::tie(seq_j, pos_j) = GSAToSeqPos(idx + mid);
        alignPos_j = gapInfo.ConvertToAlign(seq_j, pos_j);

        if(alignPos_i != alignPos_j)
        {
            hi = mid;
            continue;
        }

        cmp = SeqComp(seq_i, pos_i, seq_j, pos_j, compareLength);
        if(cmp != 0) hi = mid;
        else lo = mid + 1;
    }

    return lo;
}

int SAA::GetSequenceID(long long idx, int numBinds, int seq_i)
{
    // Update nonshared
    if(numBinds == numSeq)
    {
        int pos_i;
        std::tie(seq_i, pos_i) = GSAToSeqPos(idx);
        int alignPos_i = gapInfo.ConvertToAlign(seq_i, pos_i);
        isNS[alignPos_i] = 0;

        return 0;
    }

    if(numBinds == 1)
        return seq_i + 1;

    if(numBinds == numSeq - 1)
    {
        int totalSum = numSeq * (numSeq - 1) / 2;
        int cnt = seq_i;
        for(int i = 1; i < numBinds; ++i)
        {
            int seq, pos;
            std::tie(seq, pos) = GSAToSeqPos(idx + i);
            cnt += seq;
        }
        return totalSum - cnt + 1 + numSeq;
    }

    int ret;
    bit_vector vec(numSeq, 0);
    vec[seq_i] = 1;

    for(int i = 1; i < numBinds; ++i)
    {
        int seq, pos;
        std::tie(seq, pos) = GSAToSeqPos(idx + i);
        vec[seq] = 1;
    }

    auto it = vectorMap.find(vec);

    if(it == vectorMap.end())
    {
        ret = vectorID;
        vectorMap[vec] = vectorID++;
    }
    else
    {
        ret = it->second;
    }

    util::clear(vec);
    return ret;
}

void SAA::CalculateCount()
{
    C.width(bits::hi(saa.size())+1);
    C.resize(ALPHABETSIZE+2);
    util::set_to_value(C, 0);

    for(long long x = 0; x < saa.size() / 2; ++x)
    {
        // Special chars
        if(x == 0) continue;

        int seq, seqID = saa[2 * x], alignPos = saa[2 * x + 1];
        for(int i = 0; i < numSeq; ++i)
        {
            long long strvidx = (long long)seqID * numSeq + i;
            if(strv[strvidx])
            {
                seq = i;
                break;
            }
        }

        char ch = AlignPosToChar(seq, alignPos);
        ++C[ConvertToIndex(ch)];
    }

    C[0] = 0;
    C[1] = 1;
    for(int c = 2; c <= ALPHABETSIZE + 1; ++c) C[c] += C[c-1];

    
    for(unsigned int c = 0; c < C.size(); ++c) 
        std::cout << "C[" << c << "]: " << C[c] << std::endl;
    
}

void SAA::MakeOcc()
{
    for(int i = 0; i < TYPENUM; ++i)
    {
        for(int j = 0; j < ALPHABETSIZE; ++j)
        {
            bv[i][j].resize(saa.size() / 2);
            util::set_to_value(bv[i][j], 0);
        }
    }
    sampled.resize(saa.size() / 2);
    util::set_to_value(sampled, 0);

    int CurrOcc[ALPHABETSIZE] = {0};
    int Cluster[ALPHABETSIZE] = {0};

    for(long long i = 0; i < saa.size() / 2; ++i)
    {
        int seqID = saa[2 * i], alignPos = saa[2 * i + 1];

        if(alignPos % config.samplingDensity == 0) sampled[i] = 1;
        if(alignPos == alignLength - 1) sampled[i] = 1;
        if(alignPos == 0) continue;

        bool isCurrNS = (alignPos != alignLength) ? isNS[alignPos] : false;

        // curr char is in shared region
        if(!isCurrNS)
        {
            // Case 1
            // Both the curr char and the prev char are in shared region
            if(!isNS[alignPos - 1])
            {
                char ch = AlignPosToChar(0, alignPos - 1);
                int chIndex = ConvertToIndex(ch) - 2;
                bv[OCCID][chIndex][i] = 1;
                ++CurrOcc[chIndex];
            }
            // Case 2
            // The curr char is in shared region and the prev char is in non-shared region
            else
            {
                std::set<char> charSet;
                for(int j = 0; j < numSeq; ++j)
                {
                    if(gapInfo.ConvertToInd(j, alignPos) == 0) continue;
                    char ch = AlignPosToChar(j, alignPos - 1);
                    charSet.insert(ch);
                }

                sampled[i] = 1;

                for(auto it = charSet.begin(); it != charSet.end(); ++it)
                {
                    int chIndex = ConvertToIndex(*it) - 2;
                    int countID = ConvertToIndex(*it) - 1;

                    bv[OCCID][chIndex][i] = 1;
                    ++CurrOcc[chIndex];

                    int lf = C[countID] + CurrOcc[chIndex] - 1;
                    sampled[lf] = 1;
                }
            }
        }
        // curr char is in non-shared region
        else
        {
            std::set<char> charSet;

            bool case3 = false;
            for(int j = 0; j < numSeq; ++j)
            {
                long long strvidx = (long long)seqID * numSeq + j;
                if(strv[strvidx])
                {
                    if(gapInfo.ConvertToInd(j, alignPos) == 0)
                    {
                        sampled[i] = 0;
                        continue;
                    }

                    int prevIndPos_j = gapInfo.ConvertToInd(j, alignPos) - 1;
                    int prevAlignPos_j = gapInfo.ConvertToAlign(j, prevIndPos_j);

                    // Case 3
                    // The curr char is in non-shared region and the prev char is in shared region
                    if(!isNS[prevAlignPos_j])
                    {
                        char ch = AlignPosToChar(j, prevAlignPos_j);
                        int chIndex = ConvertToIndex(ch) - 2;

                        bv[EXTRAID][chIndex][i] = 1;
                        sampled[i] = 1;

                        if(Cluster[chIndex] % numSeq == 0)
                        {
                            Cluster[chIndex] = 0;
                            bv[OCCID][chIndex][i] = 1;
                            ++CurrOcc[chIndex];
                        }
                        ++Cluster[chIndex];

                        case3 = true;
                    }
                    // Case 4
                    // Both the curr char and the prev char are in non-shared region
                    else
                    {
                        char ch = AlignPosToChar(j, prevAlignPos_j);
                        charSet.insert(ch);
                    }
                }
            }

            bool valid = (charSet.size() + case3 > 1) ? true : false;
            if(valid) sampled[i] = 1;

            for(auto it = charSet.begin(); it != charSet.end(); ++it)
            {
                int chIndex = ConvertToIndex(*it) - 2;
                int countID = ConvertToIndex(*it) - 1;
                bv[OCCID][chIndex][i] = 1;
                ++CurrOcc[chIndex];

                if(valid)
                {
                    int lf = C[countID] + CurrOcc[chIndex] - 1;
                    sampled[lf] = 1;
                }
            }
        }
    }
}
