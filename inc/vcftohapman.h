#ifndef LH_INPUT_H
#define LH_INPUT_H

#include <interface.h>

int lh_indices_from_vcf(char* vcf_path, size_t ref_start, size_t ref_end, linearReferenceStructure** return_lr, haplotypeCohort** return_cohort);
int lh_indices_from_vcf_subset(char* vcf_path, size_t ref_start, size_t ref_end, linearReferenceStructure** return_lr, haplotypeCohort** return_cohort, size_t number);

void get_interval_bounds(const char* str, int32_t* beg, int32_t* end);

#endif
