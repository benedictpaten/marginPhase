#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <htslib/faidx.h>
#include "interface.h"
#include "vcftohapman.h"

int main(int argc, char* argv[]) {
	// -- input handling -------------------------------------------------------------------------------------------------
  if(argc != 5) {
		printf("arguments are [ref path] [interval] [vcf path] [cohort size] [number of trials]\n");
		return 1;
	}
  
  char* reference_path = argv[1];
  char* interval_str = argv[2];
  char* vcf_path = argv[3];
  size_t cohort_size = atoi(argv[4]);
  size_t n_trials = atoi(argv[5]);
  
  // static parameters
  double MUT_PEN = 2.303 * -9;
  double RECOMB_PEN = 2.303 * -6;
  size_t LINEAR_MAX_SAMPLES = 10000;
  size_t QUADRATIC_MAX_SAMPLES = 500;
  size_t RANDOM_HAPLOTYPE_GENERATIONS = 3;
  
  // Parse interval input to get start position, end position, ref sequence
	
	fprintf(stderr, "loading reference interval %s from %s\n", interval_str, reference_path);
  int32_t region_beg;
  int32_t region_end;
  get_interval_bounds(interval_str, &region_beg, &region_end);
	if(region_end != -1) {
		fprintf(stderr, "ref start is %d, ref end is %d\n", region_beg, region_end);
	}
  
  // Extract reference sequence from fasta
	faidx_t* ref_seq_fai = fai_load(reference_path);
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
	int built_index = lh_indices_from_vcf_subset(vcf_path, region_beg, region_end, &reference, &cohort, cohort_size);
	
	if(built_index == 0) {
		fprintf(stderr, "Input vcf is empty\n");
		return 1;
	} else {
		fprintf(stderr, "built haplotypeCohort from vcf %s\n", vcf_path);
	}
  
  // build penalty container
  double recombination_penalty = atof(argv[4]);
  double mutation_penalty = atof(argv[5]);
  penaltySet* penalties = penaltySet_build(recombination_penalty, mutation_penalty, cohort_size);
  
  for(size_t j = 0; j < n_trials; j++) {
  	inputHaplotype* query_ih = haplotypeCohort_random_haplo(cohort, reference, RANDOM_HAPLOTYPE_GENERATIONS, penalties, length);
    
    fprintf(stderr, "simulated read query using %d generations \n", RANDOM_HAPLOTYPE_GENERATIONS);
    
    double time_used_fast;
  	double time_used_quad;
  	double time_used_linear;
  	
  	for(size_t i = 0; i < 10; i++) {
      fastFwdAlgState* haplotype_matrix = fastFwdAlgState_initialize(reference, penalties, cohort);
    	slowFwdSolver* linear_fwd = slowFwd_initialize(reference, penalties, cohort);
    	slowFwdSolver* quadratic_fwd = slowFwd_initialize(reference, penalties, cohort);

    	struct timeval tv1, tv2, tv3, tv4;
    	gettimeofday(&tv1, NULL);
    	double result = fastFwdAlgState_score(haplotype_matrix, query_ih);
      gettimeofday(&tv2, NULL);
      time_used_fast += (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    	fprintf(stderr, "finished fast fwd alg replicate %d in %f sec\n", i, (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));
    	
      if(cohort_size > QUADRATIC_MAX_SAMPLES) {
        time_used_quad = 0;
      } else {
        gettimeofday(&tv2, NULL);
        double result_slowq = slowFwd_solve_quadratic(quadratic_fwd, query_ih);
        gettimeofday(&tv3, NULL);
        fprintf(stderr, "finished quadratic fwd alg replicate %d in %f sec\n", i, (double) (tv3.tv_usec - tv2.tv_usec) / 1000000 + (double) (tv3.tv_sec - tv2.tv_sec));
  	    time_used_quad += (double) (tv3.tv_usec - tv2.tv_usec) / 1000000 + (double) (tv3.tv_sec - tv2.tv_sec);
      }
      
      if(cohort_size > LINEAR_MAX_SAMPLES) {
        time_used_linear = 0;
      } else {
        gettimeofday(&tv3, NULL);
      	double result_slowl = slowFwd_solve_linear(linear_fwd, query_ih);
        gettimeofday(&tv4, NULL);
        time_used_linear += (double) (tv4.tv_usec - tv3.tv_usec) / 1000000 + (double) (tv4.tv_sec - tv3.tv_sec);
      	fprintf(stderr, "finished linear fwd alg replicate %d in %f sec\n", i, (double) (tv4.tv_usec - tv3.tv_usec) / 1000000 + (double) (tv4.tv_sec - tv3.tv_sec));
      }
      
      fastFwdAlgState_delete(haplotype_matrix);
    	slowFwdSolver_delete(quadratic_fwd);
    	slowFwdSolver_delete(linear_fwd);
    }

  	fprintf(stdout, "%f\t%f\t%f\t%d\t%d\t%d\n", time_used_fast/3, time_used_linear/3, time_used_quad/3, haplotypeCohort_sum_MACs(cohort), haplotypeCohort_n_sites(cohort), atoi(argv[4]), region_end - region_beg);
    
    inputHaplotype_delete(query_ih);
  }
	haplotypeCohort_delete(cohort);
	penaltySet_delete(penalties);
  return 0;
}