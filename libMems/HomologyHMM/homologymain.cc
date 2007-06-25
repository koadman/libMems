/*
 *    This file is part of HMMoC 0.5, a hidden Markov model compiler.
 *    Copyright (C) 2006 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
\*/
#include <stdlib.h>
#include "homology.h"


void run(std::string& sequence, std::string& prediction) {

  // The parameters of the model
  Params iPar;

  // Fill it with some values
  iPar.iGoUnrelated = 0.04;     // probability of going from Homologous to the Unrelated state
  iPar.iGoHomologous = 0.04;        // probability of going from Unrelated to the Homologous state
  iPar.iGoStop = 0.00001;       // probability of going from either to the End state

  // these values derived from the HOXD matrix of Chiaramonte et al 2002
  iPar.aEmitHomologous[0] = 0.4522941233821017820048370573372;		//a:a, t:t
  iPar.aEmitHomologous[1] = 0.44848066098577262840393576702477;		//c:c, g:g
  iPar.aEmitHomologous[2] = 0.0051952975642225538102111260336614;	//a:c, c:a, g:t, t:g
  iPar.aEmitHomologous[3] = 0.030436787995307373189004842137023;	//a:g, g:a, c:t, t:c
  iPar.aEmitHomologous[4] = 0.0042890021493835663524878678185533;	//a:t, t:a
  iPar.aEmitHomologous[5] = 0.0041101279232120962395233396480786;	//g:c, c:g
  iPar.aEmitHomologous[6] = 0.004461;	// gap open (from an e. coli y pestis alignment)
  // gap extend // 0.050733
  iPar.aEmitHomologous[7] = 1 - (iPar.aEmitHomologous[0] + iPar.aEmitHomologous[1] + iPar.aEmitHomologous[2] +
			iPar.aEmitHomologous[3] + iPar.aEmitHomologous[4] + iPar.aEmitHomologous[5] + iPar.aEmitHomologous[6]);

  iPar.aEmitUnrelated[0] = 0.12818742714404662781015820149872;	// a:a, t:t
  iPar.aEmitUnrelated[1] = 0.10493347210657785179017485428807;	// c:c, g:g
  iPar.aEmitUnrelated[2] = 0.11597910074937552039966694421313;	// a:c, c:a
  iPar.aEmitUnrelated[3] = iPar.aEmitUnrelated[2];
  iPar.aEmitUnrelated[4] = iPar.aEmitUnrelated[0];
  iPar.aEmitUnrelated[5] = iPar.aEmitUnrelated[1]; 
  iPar.aEmitUnrelated[6] = 0.0483;	// gap open
  // gap extend // 0.2535
  iPar.aEmitUnrelated[7] = 1 - (iPar.aEmitUnrelated[0] + iPar.aEmitUnrelated[1] + iPar.aEmitUnrelated[2] +
			iPar.aEmitUnrelated[3] + iPar.aEmitUnrelated[4] + iPar.aEmitUnrelated[5] + iPar.aEmitUnrelated[6]);

  //
  // Next, build an input emission sequence by sampling the emitted symbols according to true path
  //

  int iPathLength = sequence.length();
  char* aSequence = new char[ iPathLength ];
  memcpy(aSequence, sequence.data(), iPathLength);

  // Decode the emission sequence using Viterbi, and compute posteriors and Baum Welch counts using Forward and Backward
  HomologyDPTable *pViterbiDP, *pFWDP, *pBWDP;
  HomologyBaumWelch bw;

//  cout << "Calculating Viterbi probability..." << endl;
//  bfloat iVitProb = Viterbi_recurse(&pViterbiDP, iPar, aSequence, iPathLength );
//  cout << "Viterbi: "<<iVitProb<<endl;

//  cout << "Calculating Forward probability..." << endl;
  bfloat iFWProb = Forward(&pFWDP, iPar, aSequence, iPathLength );
//  cout << "Forward: "<<iFWProb<<endl;

//  cout << "Calculating Backward probability..." << endl;
  bfloat iBWProb = Backward(bw, pFWDP, &pBWDP, iPar, aSequence, iPathLength );
//  cout << "Backward:"<<iBWProb<<endl;
//  cout << "Calculating Viterbi path..." << endl;
  
//  Path& iViterbiPath = Viterbi_trace(pViterbiDP, iPar, aSequence, iPathLength );
/*  
  cout << "Baum Welch emission counts:"<<endl;

  cout << "\tHomologous\tUnrelated" << endl;
  for (int i=0; i<8; i++) {
    cout << i << ":\t"  
	 << bw.emissionBaumWelchCount1[i][ bw.emissionIndex("emitHomologous") ] << "\t" 
	 << bw.emissionBaumWelchCount1[i][ bw.emissionIndex("emitUnrelated") ] << endl;
  }
*/
  // Compare the true and Viterbi paths, and print the posterior probability of being in the homologous state
//  int iVHomologous = pViterbiDP->getId("homologous");

  prediction.resize(iPathLength);
  for (int i=0; i<iPathLength; i++) {

//    cout << " Viterbi:";
    double iPosterior = pFWDP->getProb("homologous",i+1)*pBWDP->getProb("homologous",i+1)/iFWProb;
//    if (iViterbiPath.toState(i) == iVHomologous) {
    if (iPosterior >= 0.5) {
      prediction[i] = 'H';
    } else {
      prediction[i] = 'N';
    }
//    cout << " " << iPosterior << endl;

  }
}
