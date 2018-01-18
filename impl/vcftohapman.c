#include <htslib/vcf.h>
#include "interface.h"
#include <stdlib.h>

int lh_indices_from_vcf_subset(char* vcf_path, size_t ref_start, size_t ref_end, siteIndex** return_lr, haplotypeCohort** return_cohort, size_t number) {
  vcfFile* cohort_vcf = vcf_open(vcf_path, "r");
  if (cohort_vcf == NULL) {
    return 0;
  } else {
    fprintf(stderr, "loading vcf %s\n", vcf_path);
  }
  
  bcf_hdr_t* cohort_hdr = bcf_hdr_read(cohort_vcf);
  bcf1_t* record = bcf_init1();
  
  size_t number_of_haplotypes;
  size_t max_number = bcf_hdr_nsamples(cohort_hdr) * 2;
  if(number > max_number || number == 0) {
    fprintf(stderr, "requested invalid number of haplotypes from cohort: wanted %d; exist %d\n", number, max_number);
    number_of_haplotypes = max_number;
  } else {
    number_of_haplotypes = number;
  }
  size_t length = ref_end - ref_start;
  
  size_t* sample_ids = malloc(number_of_haplotypes * sizeof(size_t));
  n_random_uints(sample_ids, number_of_haplotypes, max_number);
  
  siteIndex* reference = siteIndex_init_empty(ref_start);
  haplotypeCohort* cohort = haplotypeCohort_init_empty(number_of_haplotypes, reference);
  int built_initial_span = 0;
  
  fprintf(stderr, "building haplotype cohort matrix, progress: [");
  fflush(stderr);
  int stepsize = (ref_end - ref_start)/50;
  int progress;
  int steps;
  int laststep = 0;
  int stepsmade;
  int sites_added = 0;
  while(bcf_read(cohort_vcf, cohort_hdr, record) == 0) {
    size_t site_position = record->pos;
    if(site_position >= ref_start && site_position <= ref_end) {
      //TODO handle non-SNPs
      if (bcf_is_snp(record) == 1) {
        bcf_unpack(record, BCF_UN_ALL);
        if(built_initial_span == 0) {
          siteIndex_set_initial_span(reference, site_position - ref_start);
          built_initial_span = 1;
        }
        
        int64_t site_index = siteIndex_add_site(reference, site_position);
        haplotypeCohort_add_record(cohort);
        
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(cohort_hdr, record, &gt_arr, &ngt_arr);
        if(site_index >= 0) {
          for(size_t i = 0; i < number_of_haplotypes; i++) {
            int allele_index = bcf_gt_allele(gt_arr[sample_ids[i]]);
            char allele_value = record->d.allele[allele_index][0];
            haplotypeCohort_set_sample_allele(cohort, site_index, i, allele_value);
          }
        }
        sites_added++;
        progress = site_position - ref_start;
        steps = progress/stepsize;
        if(steps > laststep) {
          stepsmade = steps - laststep;
          laststep = steps;
          for(size_t i = 0; i < stepsmade; i++) {
            fprintf(stderr, "=");
            fflush(stderr);
          }
        }
        free(gt_arr);
      }
    } else if(site_position > ref_end) {
      break;
    }
  }
  fprintf(stderr, "]\n");
  
  fprintf(stderr, "number of sites %d\n", sites_added);
  
  siteIndex_calc_spans(reference, ref_end - ref_start);
  haplotypeCohort_populate_counts(cohort);

  *return_lr = reference;
  *return_cohort = cohort;
  
  free(cohort_hdr);
  free(record);
  vcf_close(cohort_vcf);
  return 1;
}

void get_interval_bounds(const char* str, int32_t* beg, int32_t* end) {
  // adapted from an htslib function
  char *s, *ep;
  size_t i, k, length, colon_pos, hyphen_pos;
    
  *beg = *end = -1;
  length = strlen(str);
  s = (char*)malloc(length+1);
  // remove space and commas
  for (i = k = 0; i < length; ++i) if (str[i] != ' ' && str[i] != ',') s[k++] = str[i];
  s[k] = 0;
    
  colon_pos = length = k;
  for (i = length; i > 0; --i) if (s[i - 1] == ':') break; // look for colon from the end
  if (i > 0) colon_pos = i - 1;

  if(colon_pos < length) {
    if(s[colon_pos + 1] == '-') {
      *beg = 0;
      i = colon_pos + 2; 
    } else {
      *beg = strtol(s + colon_pos + 1, &ep, 10);
      for (i = ep - s; i < k;) if (s[i++] == '-') break;
    }
    if(i < k) *end = strtol(s + i, &ep, 10);
  }
  free(s);
}

int lh_indices_from_vcf(char* vcf_path, size_t ref_start, size_t ref_end, siteIndex** return_lr, haplotypeCohort** return_cohort) {
  return lh_indices_from_vcf_subset(vcf_path, ref_start, ref_end, return_lr, return_cohort, 0);
}

