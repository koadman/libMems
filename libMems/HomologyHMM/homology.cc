/* Code generated by HMMoC version 0.5, Copyright (C) 2006 Gerton Lunter */
/* Generated from file homology.xml (author: Aaron Darling) on Fri Jul 06 14:27:49 EST 2007 */

/*
This file is a work based on HMMoC 0.5, a hidden Markov model compiler.
Copyright (C) 2006 by Gerton Lunter, Oxford University.

HMMoC and works based on it are free software; you can redistribute 
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

HMMOC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HMMoC; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include "homology.h"


const extern string _HomologystateId[];
const extern string _HomologyemissionId[];
const extern string _HomologytransitionId[];
const extern string _HomologytransF[];
const extern string _HomologytransT[];
const extern string _HomologytransP[];
const extern string _HomologytransE[];
const extern string _HomologyoutputId[];
const extern string _Homologyempty;
const extern int _HomologystateNum;
const extern int _HomologyemitNum;
const extern int _HomologytransNum;
const extern int _HomologyoutputNum;

HomologyDPTable::HomologyDPTable(int iLen) : isInCharge(true), stateId(_HomologystateId), emissionId(_HomologyemissionId), transitionId(_HomologytransitionId), transitionFrom(_HomologytransF), transitionTo(_HomologytransT), transitionProb(_HomologytransP), transitionEmit(_HomologytransE), outputId(_HomologyoutputId) {
    // init code:
    this->iLen = iLen;
    StateMemoryblock2.allocate(0+iLen);
    StateMemoryblock1.allocate();
    StateMemoryblock3.allocate();
}


HomologyDPTable::~HomologyDPTable() {
    if (!isInCharge) {
        // make sure data does not get deleted:
        StateMemoryblock2.absolve();
        StateMemoryblock1.absolve();
        StateMemoryblock3.absolve();
    } // if(!isInCharge)
} // destructor

const string& HomologyDPTable::getTransitionId(int id) { return id>=0 && id<_HomologytransNum ? _HomologytransitionId[id] : _Homologyempty; }
const string& HomologyDPTable::getEmissionId(int id) { return id>=0 && id<_HomologyemitNum ? _HomologyemissionId[id] : _Homologyempty; }
const string& HomologyDPTable::getStateId(int id) { return id>=0 && id<_HomologystateNum ? _HomologystateId[id] : _Homologyempty; }
const string& HomologyDPTable::getOutputId(int id) { return id>=0 && id<_HomologyoutputNum ? _HomologyoutputId[id] : _Homologyempty; }
int HomologyDPTable::getId(const string& sId)
{
    static bool bInit = false;
    static map<string,int>* pmId;
    if (!bInit) {
        pmId = new map<string,int>();
        for (int i=0;i<_HomologystateNum;i++) {
            (*pmId)[_HomologystateId[i]] = i;         // add state identifiers
        }
        for (int i=0; i<_HomologyemitNum; i++) {
            (*pmId)[_HomologyemissionId[i]] = i;      // add emission identifiers
        }
        for (int i=0; i<_HomologytransNum; i++) {  
            (*pmId)[_HomologytransitionId[i]] = i;    // add transition identifiers
        }
        for (int i=0; i<_HomologyoutputNum; i++) {
            (*pmId)[_HomologyoutputId[i]] = i;        // finally, add output identifiers
        }
        bInit = true;
    }
    map<string,int>::iterator iter = pmId->find(sId);
    if (iter == pmId->end()) {
        cout << "HomologyDPTable::getId: WARNING: identifier '" << sId << "' not found." << endl;
        return -1;
    }
    return iter->second;
}


bfloat HomologyDPTable::getProb(const string sState ,int iPos0) const
{
    return getProb(getId(sState) ,iPos0);
}


bfloat HomologyDPTable::getProb(int iState ,int iPos0) const
{
    const bfloat *CurStateMemoryblock1Secondary;
    const bfloat *CurStateMemoryblock2Secondary;
    const bfloat *CurStateMemoryblock3Secondary;
    static const int blockTable[] = {0, 1, 1, 2};
    static const int stateTable[] = {0, 0, 1, 0};
    switch (blockTable[iState]) {
        default:
        return 0.0;
        break;
        case 0:
        if ((iPos0+0>=0)&&(iPos0+0<=0)) {
            CurStateMemoryblock1Secondary = this->StateMemoryblock1.read();
            return CurStateMemoryblock1Secondary[stateTable[iState]];
        } 
        else { 
            return 0.0;
            
        }
        break;
        case 1:
        if ((iPos0+0>=1)&&(iPos0+0<=iLen+0)) {
            CurStateMemoryblock2Secondary = this->StateMemoryblock2.read((iPos0-(0))-(1));
            return CurStateMemoryblock2Secondary[stateTable[iState]];
        } 
        else { 
            return 0.0;
            
        }
        break;
        case 2:
        if ((iPos0+0>=iLen+0)&&(iPos0+0<=iLen+0)) {
            CurStateMemoryblock3Secondary = this->StateMemoryblock3.read();
            return CurStateMemoryblock3Secondary[stateTable[iState]];
        } 
        else { 
            return 0.0;
            
        }
    } // switch
} // DPTable...::getProb(int,...)

int HomologyBaumWelch::transitionIndex(string strId) const {
    map<const string,int>::const_iterator iter = mId.find(strId);
    if (iter == mId.end()) {
        cout << "HomologyBaumWelch::transitionIndex: WARNING: identifier '" << strId << "' not found." << endl;
        return -1;
    }
    return iter->second;
}


int HomologyBaumWelch::emissionIndex(string strId) const {
    map<const string,int>::const_iterator iter = mId.find(strId);
    if (iter == mId.end()) {
        cout << "HomologyBaumWelch::emissionIndex: WARNING: identifier '" << strId << "' not found." << endl;
        return -1;
    }
    return iter->second;
}


void HomologyBaumWelch::resetCounts() {
    static bool bInited = false;
    if (!bInited) {
        static const int aTemp[] = {0, 1, 2, 3, 4, 5, 6, 7};
        for (int i=0; i<8; i++) {
            transitionIdentifier0[i] = aTemp[i];
            atransitionIdx[aTemp[i]] = i;
            mId[_HomologytransitionId[aTemp[i]]] = i;
        }
    }
    for (int i=0; i<8; i++) {
        
        transitionBaumWelchCount0[i] = 0.0;
    }
    if (!bInited) {
        static const int aTemp[] = {1};
        for (int i=0; i<1; i++) {
            emissionIdentifier0[i] = aTemp[i];
            aemissionIdx[aTemp[i]] = i;
            mId[_HomologyemissionId[aTemp[i]]] = i;
        }
    }
    for (int i=0; i<1; i++) {
        
        emissionBaumWelchCount0[i] = 0.0;
    }
    if (!bInited) {
        static const int aTemp[] = {0, 2};
        for (int i=0; i<2; i++) {
            emissionIdentifier1[i] = aTemp[i];
            aemissionIdx[aTemp[i]] = i;
            mId[_HomologyemissionId[aTemp[i]]] = i;
        }
    }
    for (int i=0; i<2; i++) {
        for(int v00=0;v00<8;v00++)
        emissionBaumWelchCount1[v00][i] = 0.0;
    }
    bInited = true;
};


int HomologyBaumWelch::transitionIdentifier0[];
int HomologyBaumWelch::emissionIdentifier0[];
int HomologyBaumWelch::emissionIdentifier1[];

void HomologyBaumWelch::scaleCounts(bfloat scale) {
    for (int i=0; i<8; i++) {
        
        transitionBaumWelchCount0[i] *= scale;
    }
    for (int i=0; i<1; i++) {
        
        emissionBaumWelchCount0[i] *= scale;
    }
    for (int i=0; i<2; i++) {
        for(int v00=0;v00<8;v00++)
        emissionBaumWelchCount1[v00][i] *= scale;
    }
}


map<const string,int> HomologyBaumWelch::mId;
int HomologyBaumWelch::atransitionIdx[];
int HomologyBaumWelch::aemissionIdx[];

bfloat hmmocMax(bfloat i, bfloat j) { return i>j ? i : j; }
ostream& operator<<(ostream& os, const Path& p)
{
    for (unsigned int i=0; i<p.size(); i++) {
        os << p.fromState(i) << "--{";
            for (unsigned int j=0; j<p.emission(i).size(); j++) {
                if (j) os<<",";
                os<<p.emission(i)[j];
            }
        os<<"};"<<p.prob(i)<<"-->"<<p.toState(i)<<endl;
    }
    return os;
}

void SimplePath::addEdge(int tr, double p, vector<int>& e, int f, int t) {
    transitions.push_back(tr);
    probs.push_back(p);
    emissions.push_back(e);
    froms.push_back(f);
    tos.push_back(t);
}

void SimplePath::reverse() 
{
    std::reverse(transitions.begin(),transitions.end());
    std::reverse(probs.begin(),probs.end());
    std::reverse(emissions.begin(),emissions.end());
    std::reverse(froms.begin(),froms.end());
    std::reverse(tos.begin(),tos.end());
}

double SimplePath::prob(int i) const {
    return probs[i];
}

int SimplePath::nextFrom(int i) const {
    if (i+1 < (int)transitions.size())
    return i+1;
    else
    return -1;
}

int SimplePath::nextTo(int i) const {
    return -1;
}

const vector<int>& SimplePath::emission(int i) const {
    return emissions[i];
}

int SimplePath::fromState(int i) const {
    return froms[i];
}

int SimplePath::toState(int i) const {
    return tos[i];
}

const string _HomologystateId[] = {"start","homologous","unrelated","end"};
const string _HomologyemissionId[] = {"emitHomologous","empty","emitUnrelated"};
const string _HomologytransitionId[] = {"id$13","id$14","id$15","id$16","id$17","id$18","id$19","id$20"};
const string _HomologytransF[] = {"start","start","homologous","homologous","unrelated","unrelated","homologous","unrelated"};
const string _HomologytransT[] = {"homologous","unrelated","homologous","unrelated","unrelated","homologous","end","end"};
const string _HomologytransP[] = {"startHomologous","startUnrelated","stayHomologous","goUnrelated","stayUnrelated","goHomologous","goStopFromHomologous","goStopFromUnrelated"};
const string _HomologytransE[] = {"emitHomologous","emitUnrelated","emitHomologous","emitUnrelated","emitUnrelated","emitHomologous","empty","empty"};
const string _HomologyoutputId[] = {"sequence"};
const string _Homologyempty = "";
const int _HomologystateNum = 4;
const int _HomologyemitNum = 3;
const int _HomologytransNum = 8;
const int _HomologyoutputNum = 1;




bfloat Forward(HomologyDPTable** ppOutTable,Params iPar,char *aSeq,int iLen) {
    double iTransition[8];
    bfloat *CurStateMemoryblock2To;
    const bfloat *CurStateMemoryblock1From;
    const bfloat *CurStateMemoryblock2From;
    bfloat *CurStateMemoryblock3To;
    const bfloat *CurStateMemoryblock3From;
    int iPrevSlowCoord;
    int iSymbol[1];
    if (false && iSymbol[0] == iSymbol[0]) {}   // avoid 'unused variable' warnings
    double iEmission[2];
    /* temporary storage for ordinary (double) reals */
    register double iTempResult[1];
    /* temporary storage for extended-exponent reals */
    register bfloat iTempProb[1];
    HomologyDPTable dp(iLen);
    iTransition[0] = iPar.iStartHomologous;
    
    iTransition[1] = 1.0 - iPar.iStartHomologous;
    
    iTransition[2] = 1.0 - iPar.iGoUnrelated - iPar.iGoStopFromHomologous;
    
    iTransition[3] = iPar.iGoUnrelated;
    
    iTransition[4] = 1.0 - iPar.iGoHomologous - iPar.iGoStopFromUnrelated;
    
    iTransition[5] = iPar.iGoHomologous;
    
    iTransition[6] = iPar.iGoStopFromHomologous;
    
    iTransition[7] = iPar.iGoStopFromUnrelated;
    dp.StateMemoryblock1.write()[0] = 1.0;
    dp.StateMemoryblock1.written();
    iPrevSlowCoord = -1;
    for (int iPos0=0; iPos0<iLen+1; ++iPos0) {
        if ((iPos0+0<=0)) {
        }
        if ((iPos0+0>=1)) {
            if ((iPos0+-1>=0)) {
                iSymbol[0] = aSeq[iPos0+-1];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock2To = dp.StateMemoryblock2.write((iPos0-(0))-(1));
            iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
            iEmission[0] = iTempResult[0];
            iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
            iEmission[1] = iTempResult[0];
            if ((iPos0+-1<=0)) {
                CurStateMemoryblock1From = dp.StateMemoryblock1.read();
                CurStateMemoryblock2To[0] += (iTransition[0]*iEmission[0])*CurStateMemoryblock1From[0];
                CurStateMemoryblock2To[1] += (iTransition[1]*iEmission[1])*CurStateMemoryblock1From[0];
            }
            if ((iPos0+-1>=1)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(1))-(1));
                CurStateMemoryblock2To[0] += (iTransition[2]*iEmission[0])*CurStateMemoryblock2From[0];
                CurStateMemoryblock2To[0] += (iTransition[5]*iEmission[0])*CurStateMemoryblock2From[1];
                CurStateMemoryblock2To[1] += (iTransition[3]*iEmission[1])*CurStateMemoryblock2From[0];
                CurStateMemoryblock2To[1] += (iTransition[4]*iEmission[1])*CurStateMemoryblock2From[1];
            }
            dp.StateMemoryblock2.written();
        }
        if ((iPos0+0>=iLen+0)) {
            CurStateMemoryblock3To = dp.StateMemoryblock3.write();
            iEmission[0] = 1.0;
            if ((iPos0+0>=1)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(0))-(1));
                CurStateMemoryblock3To[0] += (iTransition[6]*iEmission[0])*CurStateMemoryblock2From[0];
                CurStateMemoryblock3To[0] += (iTransition[7]*iEmission[0])*CurStateMemoryblock2From[1];
            }
            dp.StateMemoryblock3.written();
        }
        iPrevSlowCoord = iPos0;
    }
    iPrevSlowCoord = -1;
    {
        int iPos0=iLen+0;
        if (iPos0==iPos0) {} // avoid 'unused variable' warnings
        CurStateMemoryblock3From = dp.StateMemoryblock3.read();
        iTempProb[0] = CurStateMemoryblock3From[0];
    }
    *ppOutTable = new HomologyDPTable(dp);
    // make sure tables don't get deleted
    dp.isInCharge = false;
    return iTempProb[0];
};





bfloat Backward(HomologyBaumWelch& bw,HomologyDPTable* pInTable,HomologyDPTable** ppOutTable,Params iPar,char *aSeq,int iLen) {
    const bfloat *CurStateMemoryblock3Secondary;
    double iTransition[8];
    bfloat *CurStateMemoryblock2To;
    const bfloat *CurStateMemoryblock2Secondary;
    const bfloat *CurStateMemoryblock2From;
    unsigned char alphaSymbolsitepatterns[8] = {'1', '2', '3', '4', '5', '6', '7', '8'};
    unsigned char alphaIndexsitepatterns[256];
    const bfloat *CurStateMemoryblock3From;
    bfloat *CurStateMemoryblock1To;
    const bfloat *CurStateMemoryblock1Secondary;
    const bfloat *CurStateMemoryblock1From;
    int iPrevSlowCoord;
    int iSymbol[1];
    if (false && iSymbol[0] == iSymbol[0]) {}   // avoid 'unused variable' warnings
    double iEmission[2];
    /* temporary storage for ordinary (double) reals */
    register double iTempResult[1];
    /* temporary storage for extended-exponent reals */
    register bfloat iTempProb[3];
    HomologyDPTable dp(iLen);
    HomologyDPTable dp2(*pInTable);
    // make sure tables don't get deleted
    dp2.isInCharge = false;
    iTransition[0] = iPar.iStartHomologous;
    
    iTransition[1] = 1.0 - iPar.iStartHomologous;
    
    iTransition[2] = 1.0 - iPar.iGoUnrelated - iPar.iGoStopFromHomologous;
    
    iTransition[3] = iPar.iGoUnrelated;
    
    iTransition[4] = 1.0 - iPar.iGoHomologous - iPar.iGoStopFromUnrelated;
    
    iTransition[5] = iPar.iGoHomologous;
    
    iTransition[6] = iPar.iGoStopFromHomologous;
    
    iTransition[7] = iPar.iGoStopFromUnrelated;
    for (int i=0; i<256; i++) {
        alphaIndexsitepatterns[i]=0;
    }
    for (int i=0; i<8; i++) {
        alphaIndexsitepatterns[alphaSymbolsitepatterns[i]]=i;
    }
    dp.StateMemoryblock3.write()[0] = 1.0;
    dp.StateMemoryblock3.written();
    iPrevSlowCoord = -1;
    {
        int iPos0=iLen+0;
        if (iPos0==iPos0) {} // avoid 'unused variable' warnings
        CurStateMemoryblock3Secondary = dp2.StateMemoryblock3.read();
        iTempProb[2] = CurStateMemoryblock3Secondary[0];
        bw.scaleCounts(iTempProb[2]);
    }
    iPrevSlowCoord = -1;
    for (int iPos0=(iLen+1)-1; iPos0>=0; --iPos0) {
        if ((iPos0+0>=iLen+0)) {
        }
        if ((iPos0+0>=1)) {
            if ((iPos0+0<=iLen+-1)) {
                iSymbol[0] = aSeq[iPos0+0];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock2To = dp.StateMemoryblock2.write((iPos0-(0))-(1));
            CurStateMemoryblock2Secondary = dp2.StateMemoryblock2.read((iPos0-(0))-(1));
            iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
            iEmission[0] = iTempResult[0];
            iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
            iEmission[1] = iTempResult[0];
            if ((iPos0+1<=iLen+0)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(-1))-(1));
                CurStateMemoryblock2To[1] += iTempProb[1] = (iTransition[4]*iEmission[1])*CurStateMemoryblock2From[1];
                iTempProb[1] *= CurStateMemoryblock2Secondary[1];
                bw.transitionBaumWelchCount0[4] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][1] += iTempProb[1];
                CurStateMemoryblock2To[1] += iTempProb[1] = (iTransition[5]*iEmission[0])*CurStateMemoryblock2From[0];
                iTempProb[1] *= CurStateMemoryblock2Secondary[1];
                bw.transitionBaumWelchCount0[5] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][0] += iTempProb[1];
                CurStateMemoryblock2To[0] += iTempProb[1] = (iTransition[3]*iEmission[1])*CurStateMemoryblock2From[1];
                iTempProb[1] *= CurStateMemoryblock2Secondary[0];
                bw.transitionBaumWelchCount0[3] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][1] += iTempProb[1];
                CurStateMemoryblock2To[0] += iTempProb[1] = (iTransition[2]*iEmission[0])*CurStateMemoryblock2From[0];
                iTempProb[1] *= CurStateMemoryblock2Secondary[0];
                bw.transitionBaumWelchCount0[2] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][0] += iTempProb[1];
            }
            iEmission[0] = 1.0;
            if ((iPos0+0>=iLen+0)) {
                CurStateMemoryblock3From = dp.StateMemoryblock3.read();
                CurStateMemoryblock2To[1] += iTempProb[1] = (iTransition[7]*iEmission[0])*CurStateMemoryblock3From[0];
                iTempProb[1] *= CurStateMemoryblock2Secondary[1];
                bw.transitionBaumWelchCount0[7] += iTempProb[1];
                bw.emissionBaumWelchCount0[0] += iTempProb[1];
                CurStateMemoryblock2To[0] += iTempProb[1] = (iTransition[6]*iEmission[0])*CurStateMemoryblock3From[0];
                iTempProb[1] *= CurStateMemoryblock2Secondary[0];
                bw.transitionBaumWelchCount0[6] += iTempProb[1];
                bw.emissionBaumWelchCount0[0] += iTempProb[1];
            }
            dp.StateMemoryblock2.written();
        }
        if ((iPos0+0<=0)) {
            if ((iPos0+0<=iLen+-1)) {
                iSymbol[0] = aSeq[iPos0+0];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock1To = dp.StateMemoryblock1.write();
            CurStateMemoryblock1Secondary = dp2.StateMemoryblock1.read();
            iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
            iEmission[0] = iTempResult[0];
            iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
            iEmission[1] = iTempResult[0];
            if ((iPos0+1<=iLen+0)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(-1))-(1));
                CurStateMemoryblock1To[0] += iTempProb[1] = (iTransition[1]*iEmission[1])*CurStateMemoryblock2From[1];
                iTempProb[1] *= CurStateMemoryblock1Secondary[0];
                bw.transitionBaumWelchCount0[1] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][1] += iTempProb[1];
                CurStateMemoryblock1To[0] += iTempProb[1] = (iTransition[0]*iEmission[0])*CurStateMemoryblock2From[0];
                iTempProb[1] *= CurStateMemoryblock1Secondary[0];
                bw.transitionBaumWelchCount0[0] += iTempProb[1];
                bw.emissionBaumWelchCount1[alphaIndexsitepatterns[iSymbol[0]]][0] += iTempProb[1];
            }
            dp.StateMemoryblock1.written();
        }
        iPrevSlowCoord = iPos0;
    }
    bw.scaleCounts(1.0 / iTempProb[2]);
    iPrevSlowCoord = -1;
    {
        int iPos0=0;
        if (iPos0==iPos0) {} // avoid 'unused variable' warnings
        CurStateMemoryblock1From = dp.StateMemoryblock1.read();
        iTempProb[0] = CurStateMemoryblock1From[0];
    }
    *ppOutTable = new HomologyDPTable(dp);
    // make sure tables don't get deleted
    dp.isInCharge = false;
    return iTempProb[0];
};





bfloat Viterbi_recurse(HomologyDPTable** ppOutTable,Params iPar,char *aSeq,int iLen) {
    double iTransition[8];
    bfloat *CurStateMemoryblock2To;
    const bfloat *CurStateMemoryblock2From;
    const bfloat *CurStateMemoryblock3From;
    bfloat *CurStateMemoryblock1To;
    const bfloat *CurStateMemoryblock1From;
    int iPrevSlowCoord;
    int iSymbol[1];
    if (false && iSymbol[0] == iSymbol[0]) {}   // avoid 'unused variable' warnings
    double iEmission[2];
    /* temporary storage for ordinary (double) reals */
    register double iTempResult[1];
    /* temporary storage for extended-exponent reals */
    register bfloat iTempProb[1];
    HomologyDPTable dp(iLen);
    iTransition[0] = iPar.iStartHomologous;
    
    iTransition[1] = 1.0 - iPar.iStartHomologous;
    
    iTransition[2] = 1.0 - iPar.iGoUnrelated - iPar.iGoStopFromHomologous;
    
    iTransition[3] = iPar.iGoUnrelated;
    
    iTransition[4] = 1.0 - iPar.iGoHomologous - iPar.iGoStopFromUnrelated;
    
    iTransition[5] = iPar.iGoHomologous;
    
    iTransition[6] = iPar.iGoStopFromHomologous;
    
    iTransition[7] = iPar.iGoStopFromUnrelated;
    dp.StateMemoryblock3.write()[0] = 1.0;
    dp.StateMemoryblock3.written();
    iPrevSlowCoord = -1;
    for (int iPos0=(iLen+1)-1; iPos0>=0; --iPos0) {
        if ((iPos0+0>=iLen+0)) {
        }
        if ((iPos0+0>=1)) {
            if ((iPos0+0<=iLen+-1)) {
                iSymbol[0] = aSeq[iPos0+0];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock2To = dp.StateMemoryblock2.write((iPos0-(0))-(1));
            iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
            iEmission[0] = iTempResult[0];
            iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
            iEmission[1] = iTempResult[0];
            if ((iPos0+1<=iLen+0)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(-1))-(1));
                CurStateMemoryblock2To[1] = hmmocMax( CurStateMemoryblock2To[1], (iTransition[4]*iEmission[1])*CurStateMemoryblock2From[1] );
                CurStateMemoryblock2To[1] = hmmocMax( CurStateMemoryblock2To[1], (iTransition[5]*iEmission[0])*CurStateMemoryblock2From[0] );
                CurStateMemoryblock2To[0] = hmmocMax( CurStateMemoryblock2To[0], (iTransition[3]*iEmission[1])*CurStateMemoryblock2From[1] );
                CurStateMemoryblock2To[0] = hmmocMax( CurStateMemoryblock2To[0], (iTransition[2]*iEmission[0])*CurStateMemoryblock2From[0] );
            }
            iEmission[0] = 1.0;
            if ((iPos0+0>=iLen+0)) {
                CurStateMemoryblock3From = dp.StateMemoryblock3.read();
                CurStateMemoryblock2To[1] = hmmocMax( CurStateMemoryblock2To[1], (iTransition[7]*iEmission[0])*CurStateMemoryblock3From[0] );
                CurStateMemoryblock2To[0] = hmmocMax( CurStateMemoryblock2To[0], (iTransition[6]*iEmission[0])*CurStateMemoryblock3From[0] );
            }
            dp.StateMemoryblock2.written();
        }
        if ((iPos0+0<=0)) {
            if ((iPos0+0<=iLen+-1)) {
                iSymbol[0] = aSeq[iPos0+0];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock1To = dp.StateMemoryblock1.write();
            iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
            iEmission[0] = iTempResult[0];
            iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
            iEmission[1] = iTempResult[0];
            if ((iPos0+1<=iLen+0)) {
                CurStateMemoryblock2From = dp.StateMemoryblock2.read((iPos0-(-1))-(1));
                CurStateMemoryblock1To[0] = hmmocMax( CurStateMemoryblock1To[0], (iTransition[1]*iEmission[1])*CurStateMemoryblock2From[1] );
                CurStateMemoryblock1To[0] = hmmocMax( CurStateMemoryblock1To[0], (iTransition[0]*iEmission[0])*CurStateMemoryblock2From[0] );
            }
            dp.StateMemoryblock1.written();
        }
        iPrevSlowCoord = iPos0;
    }
    iPrevSlowCoord = -1;
    {
        int iPos0=0;
        if (iPos0==iPos0) {} // avoid 'unused variable' warnings
        CurStateMemoryblock1From = dp.StateMemoryblock1.read();
        iTempProb[0] = CurStateMemoryblock1From[0];
    }
    *ppOutTable = new HomologyDPTable(dp);
    // make sure tables don't get deleted
    dp.isInCharge = false;
    return iTempProb[0];
};





Path& Viterbi_trace(HomologyDPTable* pInTable,Params iPar,char *aSeq,int iLen) {
    double iTransition[8];
    const bfloat *CurStateMemoryblock1To;
    const bfloat *CurStateMemoryblock2To;
    const bfloat *CurStateMemoryblock3To;
    int iPrevSlowCoord;
    SimplePath* pPath = new SimplePath();
    vector<int> emit;
    int iSymbol[1];
    if (false && iSymbol[0] == iSymbol[0]) {}   // avoid 'unused variable' warnings
    double iEmission[2];
    /* temporary vector storage */
    bfloat iTempVector[9];
    /* temporary int vector storage */
    int iTempIntVec[6];
    /* temporary storage for ordinary (double) reals */
    register double iTempResult[1];
    iTransition[0] = iPar.iStartHomologous;
    
    iTransition[1] = 1.0 - iPar.iStartHomologous;
    
    iTransition[2] = 1.0 - iPar.iGoUnrelated - iPar.iGoStopFromHomologous;
    
    iTransition[3] = iPar.iGoUnrelated;
    
    iTransition[4] = 1.0 - iPar.iGoHomologous - iPar.iGoStopFromUnrelated;
    
    iTransition[5] = iPar.iGoHomologous;
    
    iTransition[6] = iPar.iGoStopFromHomologous;
    
    iTransition[7] = iPar.iGoStopFromUnrelated;
    static const int stateTable[] = {1, 2, 1, 2, 2, 1, 3, 3};
    static const int stateFromTable[] = {0, 0, 1, 1, 2, 2, 1, 2};
    static const int iPos0Table[] = {1, 1, 1, 1, 1, 1, 0, 0};
    HomologyDPTable dp(*pInTable);
    // make sure tables don't get deleted
    dp.isInCharge = false;
    dp.StateMemoryblock1.write()[0] = 1.0;
    dp.StateMemoryblock1.written();
    iPrevSlowCoord = -1;
    {
        int iPos0=0;
        if (iPos0==iPos0) {} // avoid 'unused variable' warnings
        iTempIntVec[0] = 0;
        while (iTempIntVec[0] != 3) {
            iTempIntVec[1] = 2;
            if ((iPos0+0<=iLen+-1)) {
                iSymbol[0] = aSeq[iPos0+0];
            } 
            else { 
                iSymbol[0] = '1' /* gets replaced by a dummy value from the alphabet */;
                
            }
            CurStateMemoryblock1To = dp.StateMemoryblock1.read();
            CurStateMemoryblock2To = dp.StateMemoryblock2.read((iPos0-(0))-(1));
            if ((iPos0+1<=iLen+0)) {
                iTempResult[0] = iPar.aEmitHomologous[ iSymbol[0] - '1' ];
                iEmission[0] = iTempResult[0];
                iTempResult[0] = iPar.aEmitUnrelated[ iSymbol[0] - '1' ];
                iEmission[1] = iTempResult[0];
                CurStateMemoryblock2To = dp.StateMemoryblock2.read((iPos0-(-1))-(1));
                switch (iTempIntVec[0]) {
                    default:
                    break;
                    case 0:
                    iTempVector[iTempIntVec[1]] = iTransition[0]*iEmission[0]*CurStateMemoryblock2To[0];
                    iTempVector[iTempIntVec[1]+3] = iTransition[0]*iEmission[0];
                    iTempIntVec[iTempIntVec[1]++] = 0;
                    iTempVector[iTempIntVec[1]] = iTransition[1]*iEmission[1]*CurStateMemoryblock2To[1];
                    iTempVector[iTempIntVec[1]+3] = iTransition[1]*iEmission[1];
                    iTempIntVec[iTempIntVec[1]++] = 1;
                    break;
                    case 1:
                    iTempVector[iTempIntVec[1]] = iTransition[2]*iEmission[0]*CurStateMemoryblock2To[0];
                    iTempVector[iTempIntVec[1]+3] = iTransition[2]*iEmission[0];
                    iTempIntVec[iTempIntVec[1]++] = 2;
                    iTempVector[iTempIntVec[1]] = iTransition[3]*iEmission[1]*CurStateMemoryblock2To[1];
                    iTempVector[iTempIntVec[1]+3] = iTransition[3]*iEmission[1];
                    iTempIntVec[iTempIntVec[1]++] = 3;
                    break;
                    case 2:
                    iTempVector[iTempIntVec[1]] = iTransition[5]*iEmission[0]*CurStateMemoryblock2To[0];
                    iTempVector[iTempIntVec[1]+3] = iTransition[5]*iEmission[0];
                    iTempIntVec[iTempIntVec[1]++] = 5;
                    iTempVector[iTempIntVec[1]] = iTransition[4]*iEmission[1]*CurStateMemoryblock2To[1];
                    iTempVector[iTempIntVec[1]+3] = iTransition[4]*iEmission[1];
                    iTempIntVec[iTempIntVec[1]++] = 4;
                    break;
                }
            }
            CurStateMemoryblock3To = dp.StateMemoryblock3.read();
            if ((iPos0+0>=iLen+0)) {
                iEmission[0] = 1.0;
                CurStateMemoryblock3To = dp.StateMemoryblock3.read();
                switch (iTempIntVec[0]) {
                    default:
                    break;
                    case 1:
                    iTempVector[iTempIntVec[1]] = iTransition[6]*iEmission[0]*CurStateMemoryblock3To[0];
                    iTempVector[iTempIntVec[1]+3] = iTransition[6]*iEmission[0];
                    iTempIntVec[iTempIntVec[1]++] = 6;
                    break;
                    case 2:
                    iTempVector[iTempIntVec[1]] = iTransition[7]*iEmission[0]*CurStateMemoryblock3To[0];
                    iTempVector[iTempIntVec[1]+3] = iTransition[7]*iEmission[0];
                    iTempIntVec[iTempIntVec[1]++] = 7;
                    break;
                }
            }
            iTempVector[0] = 0.0;
            for (int i=2; i<iTempIntVec[1]; i++) {
                if (iTempVector[i]>iTempVector[0]) {
                    iTempVector[0]=iTempVector[i];
                    iTempIntVec[0] = i;
                }
            }
            emit.resize(1);
            emit[0] = iPos0Table[iTempIntVec[iTempIntVec[0]]];
            pPath->addEdge(iTempIntVec[iTempIntVec[0]],iTempVector[iTempIntVec[0]+3],emit,stateFromTable[iTempIntVec[iTempIntVec[0]]],stateTable[iTempIntVec[iTempIntVec[0]]]);
            iPos0 += iPos0Table[iTempIntVec[iTempIntVec[0]]];
            iTempIntVec[0] = stateTable[iTempIntVec[iTempIntVec[0]]];
        }
    }
    return *pPath;
};



/* --- end of HMMoC-generated file --- */
