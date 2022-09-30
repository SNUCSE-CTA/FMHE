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

#include "vcf.hpp"

void split(const std::string &s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
}

void VCF::OpenFile(const std::string fname)
{
    if(vcf.is_open()) vcf.close();
    vcf.open(fname, std::ifstream::in);
}

void VCF::ReadHeader()
{
    std::getline(vcf, line);
    while(line[0] == '#')
    {
        if(line.compare(1, 5, "CHROM") == 0)
        {
            std::vector<std::string> tokens;
            split(line, '\t', tokens);
            numChrom = tokens.size() - 9;
            tokens.clear();

            ReadNext(true);

            numSeq = curInfo.AllelInfos.size();

            break;
        }

        else
            std::getline(vcf, line);
    }
}

void VCF::ReadNext(bool first = false)
{
    if(vcf >> curInfo.Chrom)
    {
        curInfo.AllelInfos.clear();
        curInfo.AllelInfos.push_back(""); // for reference sequence
        curInfo.Alts.clear();

        vcf >> curInfo.Pos
            >> curInfo.ID
            >> curInfo.Ref
            >> curInfo.Alt
            >> curInfo.Qual
            >> curInfo.Filter
            >> curInfo.Info
            >> curInfo.Format;

        
        for(int i = 0; i < numChrom; ++i)
        {
            std::string s;
            vcf >> s;
            s = s.substr(0, s.find(':'));

            if(first && i == 0)
            {
                if(s.find('/') != std::string::npos) delim = '/';
                else delim = '|';
            }
            split(s, delim, curInfo.AllelInfos);
        }
    }

    // for 0-index
    curInfo.Pos -= 1;
    split(curInfo.Alt, ',', curInfo.Alts);
}

void VCF::Parsing(const std::string& vcfPath, const std::string& refPath, const std::string& outDir)
{
    // Read VCF file header
    OpenFile(vcfPath);
    ReadHeader();

    std::string ref;
    LoadSequence(refPath, ref);
    ref += '#'; // Add end symbol
    GapInfo gapInfo(numSeq);
    std::vector<int> ncPos;


    // Sequences output
    std::vector<std::string> seqs(numSeq);
    std::vector<std::ofstream> seqFiles(numSeq);
    for(int i = 0; i < numSeq; ++i)
    {
        char fname[256];
        snprintf(fname, 256, "%sS%04d", outDir.c_str(), i);
        seqFiles[i].open(fname, std::ofstream::out);
    }
    seqFiles[0] << ref;

    std::vector<int> indelSize(numSeq, 0);
    std::vector<int> endPos(numSeq, 0);
    std::vector<int> delta(numSeq, 0);

    int totalInsertion = 0;
    int startNC, endNC = -1;

    bool inNC = false;
    while(!vcf.eof())
    {
        int prevendNC = endNC;
        std::vector<int> indelLen(numSeq, 0);
        bool isNC = false, isExtend = false, isSkip = true;
        for(int i = 1; i < numSeq; ++i)
        {
            if(endPos[i] > curInfo.Pos) continue;

            seqs[i] += ref.substr(endPos[i], curInfo.Pos - endPos[i]);
            endPos[i] = curInfo.Pos;

            std::string GT = curInfo.AllelInfos[i];
            int nr = (GT.compare(".") == 0) ? -1 : atoi(GT.c_str()) - 1;
            if(nr != -1) isSkip = false;

            // Case 1: No change
            if(nr == -1)
            {
                seqs[i] += ref[curInfo.Pos];
                endPos[i] += 1;
            }
            // Case 2 : Indel
            else if(curInfo.Info.compare(0, 5, "INDEL") == 0)
            {
                int refLen = curInfo.Ref.length();
                int altLen = curInfo.Alts[nr].length();
                indelLen[i] = altLen - refLen;


                isNC = true;
                startNC = curInfo.Pos + 1;
                if(startNC <= prevendNC + 1) isExtend = true;

                seqs[i] += curInfo.Alts[nr];
                endPos[i] += refLen;
                // Insertion
                if(altLen - refLen > 0)
                {
                    endNC = std::max(endNC, curInfo.Pos);
                }
                // Deletion
                else
                {
                    endNC = std::max(endNC, curInfo.Pos + abs(indelLen[i]));
                }
            }
            // Case 3 : Subtitution
            else
            {
                seqs[i] += curInfo.Alts[nr][0];
                endPos[i] += 1;

                isNC = true;
                startNC = curInfo.Pos;
                if(startNC <= prevendNC + 1) isExtend = true;
                endNC = std::max(endNC, curInfo.Pos);
            }

            // File wrtie
            if(seqs[i].length() >= BUFFERSIZE)
            {
                seqFiles[i] << seqs[i];
                seqs[i].clear();
            }
        }

        // Exit noncommon
        if(inNC && !isExtend && !isSkip)
        {
            int enterNCRefPos = ncPos[ncPos.size() - 1] - totalInsertion;
            int maxInsertion = *(std::max_element(indelSize.begin(), indelSize.end()));

            // Reference sequence
            if(maxInsertion)
            {
                totalInsertion += maxInsertion;
                gapInfo.Insert(0, enterNCRefPos, maxInsertion);
            }

            // Individual sequences
            for(int i = 1; i < numSeq; ++i)
            {
                if(indelSize[i] != 0)
                {
                    int p = enterNCRefPos + delta[i];
                    if(indelSize[i] < 0) p -= indelSize[i];

                    gapInfo.Insert(i, p, -indelSize[i]);
                    indelSize[i] = 0;
                }
            }

            inNC = false;
            ncPos.push_back(prevendNC + totalInsertion);
        }

        for(int i = 0; i < numSeq; ++i)
        {
            indelSize[i] += indelLen[i];
            delta[i] += indelLen[i];
        }
        
        
        // Enter noncommon
        if(!inNC && isNC)
        {
            inNC = true;
            ncPos.push_back(startNC + totalInsertion);
        }

        int nowPos = curInfo.Pos;
        while(!vcf.eof() && curInfo.Pos == nowPos) ReadNext();
    }


    if(inNC)
    {
        int maxInsertion = *(std::max_element(indelSize.begin(), indelSize.end()));
        int enterNCRefPos = ncPos[ncPos.size() - 1] - totalInsertion;

        // Reference sequence
        if(maxInsertion)
        {
            totalInsertion += maxInsertion;
            gapInfo.Insert(0, enterNCRefPos, maxInsertion);
        }
        // Individual sequences
        for(int i = 1; i < numSeq; ++i)
        {
            if(indelSize[i] != 0)
            {
                int p = enterNCRefPos + delta[i];
                if(indelSize[i] < 0) p -= indelSize[i];

                gapInfo.Insert(i, p, -indelSize[i]);
                indelSize[i] = 0;
            }
        }
        inNC = false;
        ncPos.push_back(endNC + totalInsertion);
    }

    for(int i = 1; i < numSeq; ++i)
    {
        seqs[i] += ref.substr(endPos[i]);
        seqFiles[i] << seqs[i];
        seqs[i].clear();
        seqFiles[i].close();
    }

    gapInfo.SetAlignLength(ref.size() + totalInsertion);
    gapInfo.StoreGapInfo(outDir + "gapinfo");
    StoreNonCommon(outDir + "noncommon", ncPos);
}





