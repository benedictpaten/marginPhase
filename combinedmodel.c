#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <htslib/faidx.h>
#include "interface.h"
#include "vcftohapman.h"

// placeholder for actual read-partitioning hmm state
typedef struct hmmstate hmmstate;
typedef haplotypeStateNode hsn_t;

static size_t SIZEOF_PTR = sizeof(hmmstate*);
static size_t MAX_ALLELE_PAIRS = 25;
static size_t MAX_HMM_CHILDREN = 50;
// static size_t MAX_CHUNK_SIZE = ;
// static haplotypeManager* hap_manager;

// void chunk() {
//   // chunk procedure
// 
//   // ensure that this is ceil...
//   size_t n_chunks = n_old * MAX_HMM_CHILDREN * MAX_ALLELE_PAIRS / MAX_CHUNK_SIZE;
//   size_t chunk_prnt_size = ;
// 
// 
//   hmmstate*** new_x_chunks = malloc(n_chunks * SIZEOF_PTR);
//   hsn_t*** new_y_chunks = malloc(n_chunks * SIZEOF_PTR);
//   hsn_t*** new_z_chunks = malloc(n_chunks * SIZEOF_PTR);
//   size_t* chunk_sizes = malloc(n_chunks * sizeof(size_t));
// 
//   size_t last_bd = 0;
//   size_t this_chunk_prnt_size;
//   for(size_t i = 0; i < n_chunks; i++) {
//     if(i == n_chunks - 1) {
//       this_chunk_prnt_size = n_old - (n_chunks - 1) * chunk_prnt_size;
//     } else {
//       this_chunk_prnt_size = chunk_prnt_size;
//     }
//     hmmstate* temp_prnt_x = malloc(this_chunk_prnt_size * MAX_HMM_CHILDREN * SIZEOF_PTR);
//     size_t* hSN_prnt_state_offsets = malloc(this_chunk_prnt_size * sizeof(size_t));
//     
//     hmmstate** temp_new_x = 
//     hsn_t** temp_new_y =
//     hsn_t** temp_new_z =
//   }
//   size_t n_states;
//   for(size_t i = 0; i < n_chunks; i++) {
//     n_states += *chunk_sizes[i];
//   }
// 
//   return_x = malloc(n_states * SIZEOF_PTR);
//   return_y = malloc(n_states * SIZEOF_PTR);
//   return_z = malloc(n_states * SIZEOF_PTR);
// }
// 
// void build_first_combined_column(
//           hmmstate** x_first_column,        // first column hmmstates only
//           size_t n_x_first,
//           hmmstate** return_x,              // will write to this, pass in NULL
//           hsn_t** return_y,    // will write to this, pass in NULL
//           hsn_t** return_z,    // will write to this, pass in NULL
//           size_t* return_state_count,       // will overwrite this
//           double threshold) {               // to trim cross product states
//   hmmstate** new_x = malloc(n_x_first * MAX_ALLELE_PAIRS * SIZEOF_PTR);
//   hsn_t** new_y = malloc(n_x_first * MAX_ALLELE_PAIRS * SIZEOF_PTR);
//   hsn_t** new_z = malloc(n_x_first * MAX_ALLELE_PAIRS * SIZEOF_PTR);
// 
//   // hSN is first column iff child of root
//   hsn_t* root = haplotypeManager_get_root_node(hap_manager);
//   
//   // temp storage for interleaved allele pairs
//   char* next_alleles = malloc(MAX_ALLELE_PAIRS * 2);
//   int n_next_alleles;
// 
//   j = i = 0;
// 
//   for(i < n_x_first; i++) {
//   	get_next_alleles(x_opt[i], next_alleles, &n_next_alleles);
//   	for(size_t k = 0; k < n_next_alleles; k++) {
//   		new_x[j] = x_opt[i];
//   		new_y[j] = hsn_t_get_child(root, next_alleles[k * 2]);
//   		new_z[j] = hsn_t_get_child(root, next_alleles[k * 2 + 1]);
//   		if(new_y[j] != NULL && new_z[j] != NULL) {
//   			if(combined_likelihood(new_x[j], new_y[j], new_z[j]) > threshold) {
//   				j++;
//   			}
//         // otherwise overwrite index j on next iteration
//   		}
//   	}
//   }
//   *return_state_count = j;
//   return_x = malloc(return_state_count, SIZEOF_PTR);
//   return_y = malloc(return_state_count, SIZEOF_PTR);
//   return_z = malloc(return_state_count, SIZEOF_PTR);
//   memcpy(return_x, new_x, return_state_count * SIZEOF_PTR);
//   memcpy(return_y, new_y, return_state_count * SIZEOF_PTR);
//   memcpy(return_z, new_z, return_state_count * SIZEOF_PTR);
//   free(new_x);
//   free(new_y);
//   free(new_y);
// }

// TODO fill in how to extract probability from read partitioning state
double combined_likelihood(hmmstate* x, hsn_t* y, hsn_t* z) {
  return haplotypeStateNode_total_probability(y)
         * haplotypeStateNode_total_probability(z)
         // hmmstate_probability(x)
         ;
}

void append_hmm_children(hmmstate* x, hmmstate** write_to, size_t* n_writing, size_t start) {
  
}

void get_next_alleles(hmmstate* x, char* next_alleles, int* n_next_alleles) {
  
}

void build_next_combined_column(
          hmmstate** old_x,                 // last column, passed as 3 arrays; 
          hsn_t** old_y,       //    entries with same read hmm
          hsn_t** old_z,       //    state are contiguous entries
          size_t old_state_count,           
          hmmstate** return_x,              // will write to this, pass in NULL
          hsn_t** return_y,    // will write to this, pass in NULL
          hsn_t** return_z,    // will write to this, pass in NULL
          size_t* return_state_count,       // will overwrite this
          double threshold) {               // to trim cross product states
  
  // the ._opt pointer arrays correspond to the hmmstate* dimension of the 3D
  // matrix "column" we are building
  hmmstate** x_opt = malloc(old_state_count * MAX_HMM_CHILDREN * SIZEOF_PTR);
  hsn_t** y_opt_prnt = malloc(old_state_count * MAX_HMM_CHILDREN * SIZEOF_PTR);
  hsn_t** z_opt_prnt = malloc(old_state_count * MAX_HMM_CHILDREN * SIZEOF_PTR);

  size_t j, i;
  j = i = 0;
  for(i; i < old_state_count; i++) {
  	size_t n_hmm_children;
  	append_hmm_children(old_x[i],  // children of this hmm state
  			    x_opt,                 // array of all children
  		      &n_hmm_children,       // number being appended
  			    j);                    // index to start writing
    for(size_t k = 0; k < n_hmm_children; k++) {
  		y_opt_prnt[j + k] = old_y[i];
  		z_opt_prnt[j + k] = old_z[j];
  	}
  	j += n_hmm_children;
  }
  size_t n_x_opt = j;

  // expand the x_opt dimension by two hsn_t* dimensions
  hmmstate** new_x = malloc(n_x_opt * MAX_ALLELE_PAIRS * SIZEOF_PTR);
  hsn_t** new_y = malloc(n_x_opt * MAX_ALLELE_PAIRS * SIZEOF_PTR);
  hsn_t** new_z = malloc(n_x_opt * MAX_ALLELE_PAIRS * SIZEOF_PTR);

  // temp storage for interleaved allele pairs
  char* next_alleles = malloc(MAX_ALLELE_PAIRS * 2);
  int n_next_alleles;

  j = i = 0;

  for(i; i < n_x_opt; i++) {
  	get_next_alleles(x_opt[i], next_alleles, &n_next_alleles);
  	for(size_t k = 0; k < n_next_alleles; k++) {
  		new_x[j] = x_opt[i];
  		new_y[j] = haplotypeStateNode_get_child(y_opt_prnt[i], next_alleles[k * 2]);
  		new_z[j] = haplotypeStateNode_get_child(z_opt_prnt[i], next_alleles[k * 2 + 1]);
  		if(new_y[j] != NULL && new_z[j] != NULL) {
  			if(combined_likelihood(new_x[j], new_y[j], new_z[j]) > threshold) {
  				j++;
  			}
        // otherwise overwrite index j on next iteration
  		}
  	}
  }
  free(x_opt);
  free(y_opt_prnt);
  free(z_opt_prnt);
  
  *return_state_count = j;
  // return_x = new_x;
  // return_y = new_y;
  // return_z = new_z;
  return_x = malloc(*return_state_count * SIZEOF_PTR);
  return_y = malloc(*return_state_count * SIZEOF_PTR);
  return_z = malloc(*return_state_count * SIZEOF_PTR);
  memcpy(return_x, new_x, *return_state_count * SIZEOF_PTR);
  memcpy(return_y, new_y, *return_state_count * SIZEOF_PTR);
  memcpy(return_z, new_z, *return_state_count * SIZEOF_PTR);
  free(new_x);
  free(new_y);
  free(new_y);
}

int main(int argc, char* argv[]) {
  return 0;
}