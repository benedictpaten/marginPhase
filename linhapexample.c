#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <htslib/faidx.h>
#include "interface.h"
#include "vcftohapman.h"

// arguments:
// 1. reference sequence path
//      expects a FASTA
// 2. reference interval
//      expects format "chr:start-end"
// 3. vcf path
//      expects a vcf
// 4. recombination penalty
//      expects a log-scale value < 0
// 5. mutation penalty
//      expects a log-scale value < 0
// 6. trimming cutoff (log scale)
//      expects a log-scale value <= 0; 0 symbolizes no cutoff
int main(int argc, char* argv[]) {
	if(argc != 7) {
		printf("arguments are [ref path] [interval] [vcf path] [recomb penalty] [mutation penalty] [trimming cutoff] \n");
		return 1;
	}

  // Parse interval input to get start position, end position, ref sequence
  char* interval_str = argv[2];
	fprintf(stderr, "loading reference interval %s from %s\n", argv[2], argv[1]);
  int32_t ref_start;
  int32_t ref_end;
  get_interval_bounds(interval_str, &ref_start, &ref_end);
	if(ref_end != -1) {
		fprintf(stderr, "ref start is %d, ref end is %d\n", ref_start, ref_end);
	}
	faidx_t* ref_seq_fai = fai_load(argv[1]);
	char* reference_sequence = fai_fetch(ref_seq_fai, interval_str, &ref_end);
	fprintf(stderr, "loaded reference sequence of length %d\n", strlen(reference_sequence));
	if(ref_end == -1) {
		ref_end = strlen(reference_sequence);
		fprintf(stderr, "ref start is %d, ref end is %d\n", ref_start, ref_end);
	}
  
  // build index  
	linearReferenceStructure* reference = NULL;
	haplotypeCohort* cohort = NULL;
	int built_index = lh_indices_from_vcf(argv[3], ref_start, ref_end, &reference, &cohort);
	
	if(built_index == 0) {
		fprintf(stderr, "Input vcf is empty\n");
		return 1;
	} else {
		fprintf(stderr, "built haplotypeCohort from vcf %s\n", argv[3]);
	}
  
  // SIMULATE READ DP (for simplicity)
  // in general, the information needed is:
  //      1. reference start position
  //      2. number of sites
  //      3. position of sites
  //      4. the full sequence inferred for the read; can be unassigned at sites
  
  size_t read_DP_ref_start = ref_start;
  size_t* read_DP_sites = NULL;
  size_t n_read_DP_sites;
  char* read_DP_seq = NULL;
  
  double recombination_penalty = atof(argv[4]);
  double mutation_penalty = atof(argv[5]);
  double threshold = atof(argv[6]);
  double share_rate = atof(argv[7]);
  size_t cohort_size = haplotypeCohort_n_haplotypes(cohort);
  
  haplotypeCohort_sim_read_query(cohort,
                                 reference_sequence,
                                 mutation_penalty,
                                 recombination_penalty,
                                 cohort_size,
                                 share_rate,
                                 &read_DP_sites,
                                 &n_read_DP_sites,
                                 &read_DP_seq);

	haplotypeManager* hap_manager = haplotypeManager_build_int_from_index(
            reference_sequence,
            ref_end - ref_start,
            reference,
            cohort,
            mutation_penalty, 
            recombination_penalty,
            read_DP_ref_start,
            n_read_DP_sites,
            read_DP_sites,
            read_DP_seq, 
            threshold);

	haplotypeStateNode* n = haplotypeManager_get_root_node(hap_manager);
	haplotypeStateNode* options[5];
	// fills options-vector with children of n; options vector must be
	// a minimum of number of children
	haplotypeStateNode_get_next_options(n, options);

	for(int i = 0; i < haplotypeStateNode_number_of_children(n); i++) {
		n = options[i];
		// what allele does this node have?		
		char allele = haplotypeStateNode_allele(n);
		printf("%c %lf\n", allele, haplotypeStateNode_local_probability(n, hap_manager));
	}
	
	// print the whole thing
	haplotypeManager_print_transition_likelihoods(hap_manager);
	printf("\n");
	haplotypeManager_print_prefix_likelihoods(hap_manager);
		
	return 0;
}
