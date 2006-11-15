#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatHashCat.h"

//#include "libMems/MemSubsets.h"
using namespace std;
using namespace genome;
namespace mems {

// Constructor


RepeatHashCat::RepeatHashCat( )
{
}
RepeatHashCat::~RepeatHashCat( )
{
}
RepeatHashCat::RepeatHashCat(const RepeatHashCat& mh) : RepeatHash(mh){

}
SortedMerList* RepeatHashCat::GetSar(uint32 sarI) const{
	return sar_table[0];
}
//RepeatHashCat* RepeatHashCat::Clone() const{
//	return new RepeatHashCat(*this);
//}

boolean RepeatHashCat::CreateMatches(){
	if(seq_count == 1){
		MatchFinder::FindMatchSeeds();
		return true;
	}

	return false;
}

boolean RepeatHashCat::EnumerateMatches( IdmerList& match_list ){
	return HashMatch(match_list);
}


vector<uint32> RepeatHashCat::concatContigStart( void ) {
	STACK_TRACE_START
		vector<uint32> ccs = this->concat_contig_start;
		return ccs;
	STACK_TRACE_END
}
void RepeatHashCat::FindMatchesFromPosition( RepeatMatchList& ml, const vector<gnSeqI>& start_points ){
	
	
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); ++seqI ){
		if( !AddSequence( ml.sml_table[ seqI ], ml.seq_table[ seqI ] ) ){
			ErrorMsg( "Error adding " + ml.seq_filename[seqI] + "\n");
			return;
		}
	}
	MatchFinder::FindMatchSeeds( start_points );
	
	//punt: need to check here if a match spans a concatenated seq border
	//since already iterating through all matches
	uint32 ccs = this->concat_contig_start.at(0);
	uint32 seq_count = this->concat_contig_start.size();
	//vector<set<Match*, MheCompare>> mems_to_erase;
	bool mem_split = false;
	for(uint32 i=0; i < table_size; ++i)
	{

		//need to store sequence info with RepeatMatch

		//set<Match>::iterator iter = mem_table[i].begin();
		set<RepeatMatch*, MheCompare>::const_iterator mem_iter = mem_table[i].begin();
		RepeatMatch* split;
		//for each mem in mem_table
		for(;mem_iter != mem_table[i].end(); mem_iter++)
		{
			//for each occurence of mem
			for ( uint32 j = 0; j < (*mem_iter)->SeqCount(); j++)
			{
				//for each sequence boundary

				for( uint32 k = 0 ; k < seq_count; k++)
				{

					for ( uint32 l = k+1; k < seq_count; l++)
					{

						if ( absolut((*mem_iter)->Start(j)) <= this->concat_contig_start.at(l) && absolut((*mem_iter)->Start(j)) > this->concat_contig_start.at(l-1))
						{
							(*mem_iter)->FromSeq(j,l-1);
							break;

						}
					}

					//if start of mem + length > a sequence border, do not add to match list
					if( (absolut((*mem_iter)->Start(j)) <= this->concat_contig_start.at(k) && (*mem_iter)->Start(j)+(*mem_iter)->Length() >
						this->concat_contig_start.at(k)))
					{	

						//mem_table[i].erase(*mem_iter);
						//mems_to_erase.push_back(*mem_iter);
					    cout << " BORDER JUMPING MEM SPLIT" << endl;
						
						split = *mem_iter;
						(*mem_iter)->CropEnd( this->concat_contig_start.at(k) );
						split->CropStart( this->concat_contig_start.at(k) );
						mem_split = true;

						break;
					}
					
				

				}
				ml.insert(ml.end(), *mem_iter );
				if ( mem_split )
				{
					ml.insert(ml.end(), split );
					mem_split = false;
				}
					
			    break;
				
			}
		}
		//ml.insert(ml.end(), mem_table[i].begin(), mem_table[i].end() );

	}
	
}

}	// namespace mems
