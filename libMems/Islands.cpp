/*******************************************************************************
 * $Id: Islands.cpp,v 1.12 2004/04/19 23:11:19 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/Islands.h"
#include "libMems/Aligner.h"
#include "libMems/GappedAlignment.h"

using namespace std;
using namespace genome;
namespace mems {

/**
 * Identifies gaps in the alignment between pairs of sequences that are longer than
 * some number of base pairs in length.  Prints islands to an output stream
 */
void simpleFindIslands( IntervalList& iv_list, uint island_size, ostream& island_out ){
	vector< Island > island_list;
	simpleFindIslands( iv_list, island_size, island_list );
	for( size_t isleI = 0; isleI < island_list.size(); isleI++ )
	{
		Island& i = island_list[isleI];
		island_out << i.seqI << '\t' << i.leftI << '\t' << i.rightI << '\t' 
				<< i.seqJ << '\t' << i.leftJ << '\t' << i.rightJ << endl;
	}
}


void simpleFindIslands( IntervalList& iv_list, uint island_size, vector< Island >& island_list ){
	if( iv_list.size() == 0 )
		return;
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		Interval& iv = iv_list[ iv_listI ];
		gnAlignedSequences gnas;
		iv.GetAlignedSequences( gnas, iv_list.seq_table );
		uint seq_count = iv_list.seq_table.size();
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				uint columnI = 0;
				gnSeqI curI = 0;
				gnSeqI curJ = 0;
				gnSeqI lastI = 0;
				gnSeqI lastJ = 0;
				for( columnI = 0; columnI < gnas.alignedSeqsSize(); columnI++ ){
					if( gnas.sequences[ seqI ][ columnI ] != '-' )
						curI++;
					if( gnas.sequences[ seqJ ][ columnI ] != '-' )
						curJ++;
					if( toupper( gnas.sequences[ seqI ][ columnI ] ) == 
						toupper( gnas.sequences[ seqJ ][ columnI ] ) &&
						gnas.sequences[ seqJ ][ columnI ] != '-' ){
						// check for an island that was big enough
						if( curI - lastI > island_size ||
							curJ - lastJ > island_size ){
							int64 leftI = iv.Start( seqI );
							int64 rightI = leftI < 0 ? leftI - curI : leftI + curI;
							leftI = leftI < 0 ? leftI - lastI : leftI + lastI;
							int64 leftJ = iv.Start( seqJ );
							int64 rightJ = leftJ < 0 ? leftJ - curJ : leftJ + curJ;
							leftJ = leftJ < 0 ? leftJ - lastJ : leftJ + lastJ;
							Island isle;
							isle.leftI = leftI;
							isle.leftJ = leftJ;
							isle.rightI = rightI;
							isle.rightJ = rightJ;
							island_list.push_back(isle);
						}
						
						lastI = curI;
						lastJ = curJ;
					}
				}
			}
		}
	}
}

void computeSPScore( const vector<string>& alignment, const PairwiseScoringScheme& pss, vector<score_t>& scores, score_t& score )
{
	vector< score_t > cur_m_scores( alignment[0].size(), INVALID_SCORE );
	vector< score_t > cur_g_scores( alignment[0].size(), INVALID_SCORE );
	scores.resize(alignment[0].size());
	std::fill(scores.begin(), scores.end(), 0);
	score = 0;
	double w = 1;	// weight, to be determined later...
	for( size_t i = 0; i < alignment.size(); ++i )
	{
		for( size_t j = i+1; j < alignment.size(); ++j )
		{
			std::fill( cur_m_scores.begin(), cur_m_scores.end(), INVALID_SCORE );
			std::fill( cur_g_scores.begin(), cur_g_scores.end(), INVALID_SCORE );
			computeMatchScores( alignment.at(i), alignment.at(j), pss, cur_m_scores );
			computeGapScores( alignment.at(i), alignment.at(j), pss, cur_g_scores );
			for( size_t k = 0; k < cur_m_scores.size(); ++k )
			{
				score_t s = 0;
				if( cur_m_scores[k] != INVALID_SCORE )
					s += cur_m_scores[k];
				if( cur_g_scores[k] != INVALID_SCORE )
					s += cur_g_scores[k];
				scores[k] += (score_t)(w * (double)s);
			}
		}
	}
	for( size_t k = 0; k < scores.size(); ++k )
		score += scores[k];
}

//tjt: function to compute the consensus column score, consensus sequence, and cumulative consensus score from an alignment 
void computeConsensusScore( const vector<string>& alignment, const PairwiseScoringScheme& pss, vector<score_t>& scores, string& consensus, score_t& score )
{

	consensus.clear();
	std::vector< std::vector< score_t > > allscores;

	scores.resize( alignment.at(0).size() );
	std::fill(scores.begin(), scores.end(), INVALID_SCORE);

	score =	INVALID_SCORE;

	std::vector< string > nucleotides;
	nucleotides.push_back(std::string(alignment.at(0).size(),'A'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'G'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'C'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'T'));
	
	for( size_t i = 0; i < nucleotides.size(); i++)
	{
		//tjt: score alignment!
		//for each row in the alignment, compare to string of A,G,C,T and build consensus
		std::vector< score_t > consensus_scores(alignment.at(0).size(), 0);
		
		for( gnSeqI j = 0; j < alignment.size(); j++)
		{
			std::vector< score_t > tscores(alignment.at(0).size(), 0);
		
			computeMatchScores( alignment.at(j), nucleotides.at(i), pss, tscores );
			
			for( gnSeqI k = 0; k < alignment.at(j).size(); k++)
				if( tscores.at(k) != INVALID_SCORE )
					consensus_scores.at(k) += tscores.at(k);

			computeGapScores( alignment.at(j), nucleotides.at(i), pss, tscores );

			for( gnSeqI k = 0; k < alignment.at(j).size(); k++)
				if( tscores.at(k) != INVALID_SCORE )
					consensus_scores.at(k) += tscores.at(k);
			
		}
		allscores.push_back(consensus_scores);
	}
	
	//tjt: find maxvalue for each column
	// 0 = A, 1 = G, 2 = C, 3 = T
	
	std::vector< int > columnbp( alignment.at(0).size(), (std::numeric_limits<int>::min)());
	
	//for A,G,C,T
	for( size_t i = 0; i < nucleotides.size(); i++)
	{
		//for each column
		for( size_t j = 0; j < alignment.at(0).size(); j++)
		{
			if( allscores.at(i).at(j) == INVALID_SCORE )
				continue;
			if( i == 0  )
			{				
				scores.at(j) = allscores.at(i).at(j);
				columnbp.at(j) = 0;
			}
			else if (allscores.at(i).at(j) > scores.at(j))
			{
				scores.at(j) = allscores.at(i).at(j);
				columnbp.at(j) = i;
			}
		}
	}
	//update score with maxvalue from each column
	for( size_t j = 0; j < alignment.at(0).size(); j++)
	{
		if( scores.at(j) != INVALID_SCORE )
			score += scores.at(j);
		if (columnbp.at(j) == 0)
			consensus.append("A");
		else if (columnbp.at(j) == 1)
			consensus.append("G");
		else if (columnbp.at(j) == 2)
			consensus.append("C");
		else if (columnbp.at(j) == 3)
			consensus.append("T");
	
	}
}

void computeMatchScores( const string& seq1, const string& seq2, const PairwiseScoringScheme& scoring, vector<score_t>& scores )
{
	scores.resize( seq1.size() );
	std::fill(scores.begin(), scores.end(), INVALID_SCORE);
	const uint8* table = SortedMerList::BasicDNATable();

	for (unsigned uColIndex = 0; uColIndex < seq1.size(); ++uColIndex)
	{
		char c1 = seq1[uColIndex];
		char c2 = seq2[uColIndex];
		if( c1 == '-' || c2 == '-' )
			continue;
		unsigned uLetter1 = table[c1];
		unsigned uLetter2 = table[c2];

		score_t scoreMatch = scoring.matrix[uLetter1][uLetter2];
		scores[uColIndex] = scoreMatch;
	}
}

void computeGapScores( const string& seq1, const string& seq2, const PairwiseScoringScheme& scoring, vector<score_t>& scores )
{
	scores.resize(seq1.size());

	bool bGapping1 = false;
	bool bGapping2 = false;
	score_t gap_open_score = scoring.gap_open;
	score_t gap_extend_score = scoring.gap_extend;
	score_t term_gap_score = gap_open_score;

	unsigned uColCount = seq1.size();
	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < seq1.size(); ++uColIndex)
	{
		bool bGap1 = seq1[uColIndex] == '-';
		bool bGap2 = seq2[uColIndex] == '-';
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = seq1[iColIndex] == '-';
		bool bGap2 = seq2[iColIndex] == '-';
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

	unsigned gap_left_col = 0;
	score_t cur_gap_score = 0;
	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		bool bGap1 = seq1[uColIndex] == '-';
		bool bGap2 = seq2[uColIndex] == '-';

		if (bGap1 && bGap2)
			continue;

		if (bGap1)
			{
			if (!bGapping1)
				{
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					cur_gap_score += term_gap_score;
				}else{
					cur_gap_score += gap_open_score;
					}
				bGapping1 = true;
				}
			else
				{
				cur_gap_score += gap_extend_score;
				}
			continue;
			}

		else if (bGap2)
			{
			if (!bGapping2)
				{
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					cur_gap_score += term_gap_score;
				}else{
					cur_gap_score += gap_open_score;
					}
				bGapping2 = true;
				}
			else
				{
				cur_gap_score += gap_extend_score;
				}
			continue;
			}

		if( (bGapping1 || bGapping2) )
		{
			score_t valid_cols = 0;
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColIndex; ++uGapIndex )
				if( seq1[uGapIndex] != '-' || seq2[uGapIndex] != '-' )
					valid_cols++;
			// spread the total gap penalty evenly across all columns
			score_t per_site_penalty = cur_gap_score / valid_cols;
			score_t extra = cur_gap_score - (per_site_penalty * valid_cols);
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColIndex; ++uGapIndex )
			{
				if( seq1[uGapIndex] == '-' && seq2[uGapIndex] == '-' )
					continue;
				if( scores[uGapIndex] != INVALID_SCORE )
				{
					genome::breakHere();
					cerr << "asdgohasdoghasodgh\n";
				}
				scores[uGapIndex] = per_site_penalty;
			}
			if( scores[gap_left_col] == INVALID_SCORE )
			{
				cerr << "crap!\n";
				genome::breakHere();
			}
			scores[gap_left_col] += extra;
			gap_left_col = (std::numeric_limits<unsigned>::max)();
			cur_gap_score = 0;
		}
		bGapping1 = false;
		bGapping2 = false;
		}

	if (bGapping1 || bGapping2)
		{
		cur_gap_score -= gap_open_score;
		cur_gap_score += term_gap_score;

		score_t valid_cols = 0;
		for( unsigned uGapIndex = gap_left_col; uGapIndex < uColCount; ++uGapIndex )
			if( seq1[uGapIndex] != '-' || seq2[uGapIndex] != '-' )
				valid_cols++;
		// spread the total gap penalty evenly across all columns
		score_t per_site_penalty = cur_gap_score / valid_cols;
		score_t extra = cur_gap_score - (per_site_penalty * valid_cols);
		for( unsigned uGapIndex = gap_left_col; uGapIndex < uColCount; ++uGapIndex )
		{
			if( seq1[uGapIndex] == '-' && seq2[uGapIndex] == '-' )
				continue;
			scores[uGapIndex] = per_site_penalty;
		}
		if( valid_cols > 0 )
		{
			if( scores[gap_left_col] == INVALID_SCORE )
			{
				cerr << "crap!\n";
				genome::breakHere();
			}
			scores[gap_left_col] += extra;
		}
	}
}


/**
 * Identifies stretches of alignment existing in all sequences that doesn't
 * contain a gap larger than a particular size.  Such regions are considered
 * the backbone of the alignment.
 */
void simpleFindBackbone( IntervalList& iv_list, uint backbone_size, uint max_gap_size, vector< GappedAlignment >& backbone_regions ){
	if( iv_list.size() == 0 )
		return;
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		Interval& iv = iv_list[ iv_listI ];
		gnAlignedSequences gnas;
		uint seqI;
		uint seq_count = iv_list.seq_table.size();
		vector< int64 > positions( seq_count );
		vector< int64 > starts( seq_count );
		vector< int64 > ends( seq_count );
		vector< uint > gap_size( seq_count, 0 );
		uint seqJ;
		gnSeqI bb_start_col = 0;
		gnSeqI bb_end_col = 0;
		GappedAlignment cur_backbone( seq_count, 0 );
		
		// initialize positions and starts
		for( seqI = 0; seqI < seq_count; seqI++ ){
			positions[ seqI ] = iv_list[ iv_listI ].Start( seqI );
			if( positions[ seqI ] < 0 )
				positions[ seqI ] -= iv_list[ iv_listI ].Length( seqI ) + 1;
		}
		starts = positions;
		ends = positions;

		iv.GetAlignedSequences( gnas, iv_list.seq_table );
		bool backbone = true;	// assume we are starting out with a complete alignment column
		uint columnI = 0;
		vector< int64 > prev_positions;
		for( ; columnI < gnas.alignedSeqsSize(); columnI++ ){
			bool no_gaps = true;
			prev_positions = positions;
			for( seqI = 0; seqI < seq_count; seqI++ ){
				char cur_char = gnas.sequences[ seqI ][ columnI ];
				if( cur_char != '-' && toupper(cur_char) != 'N' ){
					if( gap_size[ seqI ] > max_gap_size && backbone ){
						// end a stretch of backbone here only
						// if the backbone meets size requirements in each
						// sequence.
						for( seqJ = 0; seqJ < seq_count; seqJ++ ){
							if( ends[ seqJ ] - starts[ seqJ ] < backbone_size ){
								break;
							}
						}
						if( seqJ == seq_count ) {
							// it's a legitimate stretch of backbone
							backbone_regions.push_back( cur_backbone );
							uint bbI = backbone_regions.size() - 1;
							vector< string > aln_mat( seq_count );
							for( seqJ = 0; seqJ < seq_count; seqJ++ ){
								if( starts[ seqJ ] < 0 )
									backbone_regions[ bbI ].SetStart( seqJ, ends[ seqJ ] + 1);
								else
									backbone_regions[ bbI ].SetStart( seqJ, starts[ seqJ ] );
								backbone_regions[ bbI ].SetLength( ends[ seqJ ] - starts[ seqJ ], seqJ );
								aln_mat[ seqJ ] = gnas.sequences[ seqJ ].substr( bb_start_col, bb_end_col - bb_start_col + 1);
							}
							backbone_regions[ bbI ].SetAlignment(aln_mat);
							
						}
						// we either just finished backbone or a short area that didn't
						// qualify as backbone
						// look for a new backbone region
						backbone = false;
					}
					positions[ seqI ]++;
					gap_size[ seqI ] = 0;
				}else{
					gap_size[ seqI ]++;
					no_gaps = false;
				}
			}
			if( no_gaps ){
				bb_end_col = columnI;
				ends = positions;
				if( !backbone ){
					starts = prev_positions;
					bb_start_col = columnI;
					backbone = true;
				}
			}
		}

		// check for backbone one last time
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( ends[ seqJ ] - starts[ seqJ ] < backbone_size ){
				break;
			}
		}
		if( seqJ == seq_count ) {
			// it's a legitimate stretch of backbone
			backbone_regions.push_back( cur_backbone );
			uint bbI = backbone_regions.size() - 1;
			vector< string > aln_mat( seq_count );
			for( seqJ = 0; seqJ < seq_count; seqJ++ ){
				if( starts[ seqJ ] < 0 )
					backbone_regions[ bbI ].SetStart( seqJ, ends[ seqJ ] + 1);
				else
					backbone_regions[ bbI ].SetStart( seqJ, starts[ seqJ ] );
				backbone_regions[ bbI ].SetLength( ends[ seqJ ] - starts[ seqJ ], seqJ );
				aln_mat[ seqJ ] = gnas.sequences[ seqJ ].substr( bb_start_col, bb_end_col - bb_start_col + 1);
			}
			backbone_regions[ bbI ].SetAlignment( aln_mat );
		}
	}
}


void outputBackbone( const vector< GappedAlignment >& backbone_regions, ostream& backbone_out ){
	for( uint bbI = 0; bbI < backbone_regions.size(); bbI++ ){
		for( uint seqJ = 0; seqJ < backbone_regions[ bbI ].SeqCount(); seqJ++ ){
			if( seqJ > 0 )
				backbone_out << '\t';
			int64 bb_rend = backbone_regions[ bbI ].Start( seqJ );
			if( backbone_regions[ bbI ].Start( seqJ ) < 0 )
				bb_rend -= (int64)backbone_regions[ bbI ].Length( seqJ );
			else
				bb_rend += (int64)backbone_regions[ bbI ].Length( seqJ );
			backbone_out << backbone_regions[ bbI ].Start( seqJ ) << '\t' <<  bb_rend;
		}
		backbone_out << endl;
	}
}


// always return the left end of the one to the left and the right of the one to the right

void getGapBounds( vector<gnSeqI>& seq_lengths, vector< LCB >& adjacencies, uint seqJ, int leftI, int rightI, int64& left_start, int64& right_start ){
	if( rightI != -1 )
		right_start = absolut( adjacencies[ rightI ].left_end[ seqJ ] );
	else
		right_start = seq_lengths[seqJ] + 1;
	
	if( leftI != -1 )
		left_start = absolut( adjacencies[ leftI ].right_end[ seqJ ] );
	else
		left_start = 1;
}


void addUnalignedIntervals( IntervalList& iv_list, set< uint > seq_set, vector<gnSeqI> seq_lengths ){
	vector< LCB > adjacencies;
	vector< int64 > weights;
	uint lcbI;
	uint seqI;
	if( seq_lengths.size() == 0 )
		for( seqI = 0; seqI < iv_list.seq_table.size(); seqI++ )
			seq_lengths.push_back(iv_list.seq_table[seqI]->length());

	uint seq_count = seq_lengths.size();


	if( seq_set.size() == 0 )
	{
		// if an empty seq set was passed then assume all seqs
		// should be processed
		for( seqI = 0; seqI < seq_count; seqI++ )
			seq_set.insert( seqI );
	}
	
	weights = vector< int64 >( iv_list.size(), 0 );
	computeLCBAdjacencies_v2( iv_list, weights, adjacencies );

	vector< int > rightmost;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		rightmost.push_back( -1 );
	}
	for( lcbI = 0; lcbI <= adjacencies.size(); lcbI++ ){
		set< uint >::iterator seq_set_iterator = seq_set.begin();
		for( ; seq_set_iterator != seq_set.end(); seq_set_iterator++ ){
			seqI = *seq_set_iterator;
			// scan left
			int leftI;
			if( lcbI < adjacencies.size() ){
// left is always to the left!!
				leftI = adjacencies[ lcbI ].left_adjacency[ seqI ];
			}else
				leftI = rightmost[ seqI ];

			int rightI = lcbI < adjacencies.size() ? lcbI : -1;
// right is always to the right!!
			if( lcbI < adjacencies.size() )
				if( adjacencies[ lcbI ].right_adjacency[ seqI ] == -1 )
					rightmost[ seqI ] = lcbI;
			
			int64 left_start, right_start;
			getGapBounds( seq_lengths, adjacencies, seqI, leftI, rightI, left_start, right_start );
			int64 gap_len =  absolut( right_start ) - absolut( left_start );
			if( gap_len > 0 ){
				Match mm( seq_count );
				Match* m = mm.Copy();
				for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
					m->SetStart( seqJ, 0 );
				}
				m->SetStart( seqI, left_start );
				m->SetLength( gap_len );
				vector<AbstractMatch*> tmp(1, m);
				iv_list.push_back( Interval(tmp.begin(), tmp.end()) );
				m->Free();
			}
		}
	}
}


void findIslandsBetweenLCBs( IntervalList& iv_list, uint island_size, ostream& island_out ){
	IntervalList iv_list_tmp = iv_list;
	addUnalignedIntervals( iv_list_tmp );
	uint seq_count = iv_list.seq_table.size();
	
	for( int ivI = iv_list.size(); ivI < iv_list_tmp.size(); ivI++ ){
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			if( iv_list_tmp[ ivI ].Length( seqI ) < island_size )
				continue;

			// this is an island, write the LCB island out
			gnSeqI left_end = absolut( iv_list_tmp[ ivI ].Start( seqI ) );
			gnSeqI right_end = left_end + iv_list_tmp[ ivI ].Length( seqI ) - 1;
			island_out << "LCB island:\t" << seqI << '\t' << left_end << '\t' << right_end << endl;
		}
	}
}


void allocateDetailList( MatchList& mlist, vector< pair< uint64, uint64 > >& detail_list ){
	IntervalList iv_list;
	iv_list.seq_table = mlist.seq_table;
	iv_list.seq_filename = mlist.seq_filename;
	allocateDetailList( iv_list, detail_list );
}

void getLCBDetailList( MatchList& mlist, vector< pair< uint64, uint64 > >& detail_list ){
	if( mlist.size() > 0 )
	{
		Interval iv( mlist.begin(), mlist.end() );
		getLCBDetailList( iv, detail_list );
	}
}


void allocateDetailList( IntervalList& iv_list, vector< pair< uint64, uint64 > >& detail_list ){
	uint seq_count = iv_list.seq_table.size();
	uint64 max_types = 1 << seq_count; // (uint64) pow( (double)2, (double)seq_count );
	detail_list.clear();
	uint64 seqI;
	for( seqI = 0; seqI < max_types; seqI++ )
		detail_list.push_back( pair< uint64, uint64 >( seqI, 0 ) );

	uint64 match_bits = 1;
	for( seqI = seq_count; seqI > 0; seqI-- ){
		detail_list[ match_bits ].second = iv_list.seq_table[ seqI - 1 ]->length();
		match_bits <<= 1;
	}
}

void getLCBDetailList( Interval& iv, vector< pair< uint64, uint64 > >& detail_list ) {
	cerr << "re-implement getLCBDetailList()\n";
	throw "Re-implement me!\n";
/*	if( iv.matches.size() == 0 )
		return;
	uint seq_count = iv.matches[ 0 ]->SeqCount();
	uint seqI;
	for( uint matchI = 0; matchI < iv.matches.size(); matchI++ ){
		//get the sequence matching entry into match bits
		const GappedAlignment* cr = dynamic_cast< const GappedAlignment* >( iv.matches[ matchI ] );
		const Match* match = dynamic_cast< const Match* >( iv.matches[ matchI ] );
		if( match != NULL ){
			uint64 match_bits = 0;
			for( seqI = 0; seqI < seq_count; seqI++ ){
				match_bits <<= 1;
				if((*match)[ seqI ] != NO_MATCH){
					match_bits |= 0x1;
					uint64 cur_matchnum = 1;
					cur_matchnum <<= (seq_count - seqI - 1);
					detail_list[ cur_matchnum ].second -= match->Length();
				}
			}
			detail_list[ match_bits ].second += match->Length();
		}else if( cr != NULL ){
			vector< gnSequence* > seq_table;
			const vector< string >& align_matrix = GetAlignment( *cr, seq_table );
			for( gnSeqI charI = 0; charI < cr->AlignmentLength(); charI++ ){
				for( seqI = 0; seqI < seq_count; seqI++ ){
					uint64 cur_matchnum = 0;
					uint seqJ = 0;
					uint64 match_bits = 1;
					match_bits <<= seqI;
					boolean matched = false;
					for( seqJ = 0; seqJ < seq_count; seqJ++ ){
						if( seqJ > seqI )
							match_bits <<= 1;
						if( toupper( align_matrix[ seqI ][ charI ] ) == 
							toupper( align_matrix[ seqJ ][ charI ] ) )
							if( seqI < seqJ && align_matrix[ seqI ][ charI ] != '-' ){
								cur_matchnum = 1;
								cur_matchnum <<= seq_count - seqJ - 1;
								detail_list[ cur_matchnum ].second -= 1;
								match_bits |= 0x1;
								matched = true;
							}else if( seqI > seqJ )
								break;
						
					}
					if( seqJ == seq_count && matched ){
						cur_matchnum = 1;
						cur_matchnum <<= seq_count - seqI - 1;
						detail_list[ cur_matchnum ].second -= 1;
						detail_list[ match_bits ].second += 1;
					}
						
				}
			}
		}
	}
	*/
}


void PrintDetailList( uint seq_count, const vector< pair< uint64, uint64 > >& detail_list, ostream& os, boolean skip_zeros ){
	uint64 detailMax = (uint64) pow((double)2, (double)seq_count);
	for(uint64 detailI = 0; detailI < detailMax; detailI++){
		string seq_string;
		if( skip_zeros && detail_list[detailI].second == 0 )
			continue;
		uint64 curDetail = detail_list[ detailI ].first;
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			if(curDetail & 0x1)
				seq_string = "1 " + seq_string;
			else
				seq_string = "0 " + seq_string;
			curDetail >>= 1;
		}
		os << seq_string << '\t' << detail_list[detailI].second << '\n';
	}
}

void DistanceMatrix( const MatchList& mlist, NumericMatrix<double>& distance ){
	IdentityMatrix(mlist, mlist.seq_table, distance );
	TransformDistanceIdentity( distance );
}

void TransformDistanceIdentity( NumericMatrix<double>& identity ){
	for( int i = 0; i < identity.cols(); i++ ){
		for( int j = 0; j < identity.rows(); j++ ){
			identity( i, j ) = 1 - identity( i, j );
		}
	}
}

void DistanceMatrix( uint seq_count, const vector< pair< uint64, uint64 > >& detail_list, NumericMatrix<double>& distance ){
	distance = NumericMatrix<double>( seq_count, seq_count );
	distance.init( 0 );
	uint seqI;
	uint seqJ;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		uint64 seqI_mask = 1;
		seqI_mask <<= seq_count - seqI - 1;
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			uint64 seqJ_mask = 1;
			seqJ_mask <<= seq_count - seqJ - 1;
			for( uint pairI = 0; pairI < detail_list.size(); pairI++ ){
				if( (detail_list[ pairI ].first & seqI_mask) != 0 &&
					(detail_list[ pairI ].first & seqJ_mask) != 0 ){
					distance( seqI, seqJ ) += detail_list[ pairI ].second;
				}
			}
		}
	}
	
	for( seqI = 0; seqI < seq_count; seqI++ ){
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( seqI == seqJ )
				continue;
			double avg_length = ( distance( seqI, seqI ) + distance( seqJ, seqJ ) ) / 2;
			distance( seqI, seqJ ) = 1.0 - ( distance( seqI, seqJ ) / avg_length );
			if( !(distance( seqI, seqJ ) == distance( seqI, seqJ )) ){
				distance( seqI, seqJ ) = 1.0;
			}
		}
	}

	// set the diagonal identical to itself
	for( seqI = 0; seqI < seq_count; seqI++ )
		distance( seqI, seqI ) = 0;
}


/*
void Aligner::GetLCBIdentityMatrix( MatchList& lcb, NumericMatrix< double >& identity, NumericMatrix< double >& range ){
	NumericMatrix< double > density( seq_count, seq_count );
	NumericMatrix< double > m_range( seq_count, seq_count );
	NumericMatrix< double > last_coord( seq_count, seq_count );
	uint seqI, seqJ;
	density.init( 0 );
	m_range.init( 0 );
	last_coord.init( 0 );

	vector< Match* >::iterator match_iter = lcb.begin();
	for( ; match_iter != lcb.end(); match_iter++ ){
		for( seqI = 0; seqI < seq_count; seqI++ ){
			if( (*match_iter)->Start( seqI ) == NO_MATCH )
				continue;
			for( seqJ = 0; seqJ < seq_count; seqJ++ ){
				if( (*match_iter)->Start( seqJ ) == NO_MATCH || seqI == seqJ )
					continue;
				last_coord( seqI, seqJ ) = (*match_iter)->End( seqI );
				density( seqI, seqJ ) += (*match_iter)->Length();
				if( m_range( seqI, seqJ ) == 0 ){
					m_range( seqI, seqJ ) = (*match_iter)->Start( seqI );
					if( m_range( seqI, seqJ ) < 0 )
						m_range( seqI, seqJ ) -= (*match_iter)->Length() - 1;
				}
			}
		}
	}
	for( seqI = 0; seqI < seq_count; seqI++ ){
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( seqI == seqJ )
				density( seqI, seqJ ) = 1;
			else if( last_coord( seqI, seqJ ) != 0 ){
				m_range( seqI, seqJ ) = absolut( last_coord( seqI, seqJ ) - m_range( seqI, seqJ ) ) + 1;
				density( seqI, seqJ ) = density( seqI, seqJ ) / (double)m_range( seqI, seqJ );
			}
		}
	}
	range = m_range;
	identity = density;
}
*/

}
