#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
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
// 7. share rate (what proportion of cohort sites are sites in the simulated read input)
int main(int argc, char* argv[]) {
	if(argc != 7) {
		printf("arguments are [ref path] [interval] [vcf path] [recomb penalty] [mutation penalty] [trimming cutoff]\n");
		return 1;
	}

  // Parse interval input to get start position, end position, ref sequence
  // Extract reference sequence from fasta
	char* interval_str = argv[2];
	fprintf(stderr, "loading reference interval %s from %s\n", argv[2], argv[1]);
  int32_t region_beg;
  int32_t region_end;
  get_interval_bounds(interval_str, &region_beg, &region_end);
	if(region_end != -1) {
		fprintf(stderr, "ref start is %d, ref end is %d\n", region_beg, region_end);
	}
	faidx_t* ref_seq_fai = fai_load(argv[1]);
	int32_t length = region_end;
	char* ref_seq = fai_fetch(ref_seq_fai, interval_str, &length);
	fprintf(stderr, "loaded reference sequence of length %lu\n", strlen(ref_seq));
	if(region_end == -1) {
		region_end = strlen(ref_seq);
		fprintf(stderr, "ref start is %d, ref end is %d\n", region_beg, region_end);
	}
  
  // build index  
	siteIndex* reference = NULL;
	haplotypeCohort* cohort = NULL;
	int built_index = lh_indices_from_vcf(argv[3], region_beg, region_end, &reference, &cohort);
	
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
  
  size_t read_region_beg = region_beg;
  size_t* read_sites = NULL;
  size_t n_read_sites = 0;
  char* read_seq = (char*)malloc(strlen(ref_seq) + 1);
	strcpy(read_seq, ref_seq);
	char* r_alleles_1 = NULL;
  char* r_alleles_2 = NULL;
  
  double recombination_penalty = atof(argv[4]);
  double mutation_penalty = atof(argv[5]);
  double threshold = atof(argv[6]);
  size_t cohort_size = haplotypeCohort_n_haplotypes(cohort);
  
	clock_t start, end;
  double cpu_time_used;

	haplotypeManager* hap_manager = haplotypeManager_build_from_idx(
            ref_seq,
            region_end - region_beg,
            reference,
            cohort,
            mutation_penalty, 
            recombination_penalty,
            read_region_beg,
            n_read_sites,
            read_sites,
            read_seq);
						
	haplotypeManager_init_opt_idx(hap_manager,
														  	r_alleles_1,
																r_alleles_2);
		
	start = clock();
	haplotypeManager_build_tree_interval(hap_manager, threshold);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// printf("%f sec to build tree\n", cpu_time_used);
	haplotypeManager_print_terminal_nodes(hap_manager);
	// haplotypeManager_print_prefix_likelihoods(hap_manager);
	
	haplotypeManager_delete(hap_manager);
	free(read_sites);
	free(read_seq);
	free(ref_seq);
	free(ref_seq_fai);

	return 0;
}
