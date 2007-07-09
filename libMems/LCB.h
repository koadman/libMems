#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __LCB_h__
#define __LCB_h__

#include <vector>
#include <libGenome/gnDefs.h>

namespace mems {

/** 
 * This class is used to track relationships between LCBs during the LCB determination process.
 */
class LCB{
public:
	std::vector< int64 > left_end;	/**< The left end position of the LCB in each sequence */
	std::vector< int64 > right_end;  /**< The right end position of the LCB in each sequence */
	std::vector< uint > left_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
	std::vector< uint > right_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
	double weight;		/**< The weight (or coverage) of this LCB */
	bool to_be_deleted;	/**< set to true if this LCB is about to be deleted, but the deletion hasn't yet been processed */
};

} // namespace mems


#endif  // __LCB_h__

