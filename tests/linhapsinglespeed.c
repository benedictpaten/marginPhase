#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <htslib/faidx.h>
#include "interface.h"
#include "vcftohapman.h"

int main(int argc, char* argv[]) {
	if(argc != 8) {
		printf("arguments are [ref path] [interval] [vcf path] [recomb penalty] [mutation penalty] [cohort size]\n");
		return 1;
	}

  // Parse interval input to get start position, end position, ref sequence
  // Extract reference sequence from fasta
	char* interval_str = argv[2];
	fprintf(stdout, "loading reference interval %s from %s\n", argv[2], argv[1]);
  int32_t region_beg;
  int32_t region_end;
  get_interval_bounds(interval_str, &region_beg, &region_end);
	if(region_end != -1) {
		fprintf(stdout, "ref start is %d, ref end is %d\n", region_beg, region_end);
	}
	faidx_t* ref_seq_fai = fai_load(argv[1]);
	int32_t length = region_end;
	char* ref_seq = fai_fetch(ref_seq_fai, interval_str, &length);
	fprintf(stdout, "loaded reference sequence of length %lu\n", strlen(ref_seq));
	if(region_end == -1) {
		region_end = strlen(ref_seq);
		fprintf(stdout, "ref start is %d, ref end is %d\n", region_beg, region_end);
	}
  
  // build index  
	linearReferenceStructure* reference = NULL;
	haplotypeCohort* cohort = NULL;
	int built_index = lh_indices_from_vcf_subset(argv[3], region_beg, region_end, &reference, &cohort, atoi(argv[7]));
	
	if(built_index == 0) {
		fprintf(stdout, "Input vcf is empty\n");
		return 1;
	} else {
		fprintf(stdout, "built haplotypeCohort from vcf %s\n", argv[3]);
	}

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
	fprintf(stdout, "%d cohort size\n", cohort_size);
  
  penaltySet* penalties = penaltySet_build(recombination_penalty,
                               mutation_penalty,
                               cohort_size);

  haplotypeCohort_sim_read_query_2(cohort,
                                 ref_seq,
                                 mutation_penalty,
                                 recombination_penalty,
                                 0,
                                 &read_sites,
                                 &n_read_sites,
                                 read_seq,
                                 &r_alleles_1,
                                 &r_alleles_2);
  for(size_t i = 0; i < n_read_sites; i++) {
    read_seq[read_sites[i]] = r_alleles_1[i];
  }
  printf("%s\n", read_seq);
  
  inputHaplotype* input_haplotype = 
              inputHaplotype_build(ref_seq, read_seq, reference, region_beg);
  haplotypeMatrix* haplotype_matrix = 
              haplotypeMatrix_initialize(reference, penalties, cohort);
							
								clock_t start, end;
							  double cpu_time_used;
									start = clock();

  double result = haplotypeMatrix_score(haplotype_matrix, input_haplotype);
  	end = clock();
			cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;


  printf("%f\n", result);
	
	fprintf(stderr, "%f %d %d\n", cpu_time_used, atoi(argv[7]), region_end - region_beg);
  
  inputHaplotype_delete(input_haplotype);
  haplotypeMatrix_delete(haplotype_matrix);
  penaltySet_delete(penalties);
  return 0;
}