#include <interface.h>
#include <htslib/vcf.h>

int lh_indices_from_vcf(char* vcf_path, size_t ref_start, size_t ref_end, linearReferenceStructure** return_lr, haplotypeCohort** return_cohort) {
  vcfFile* cohort_vcf = vcf_open(vcf_path, "r");
  if (cohort_vcf == NULL) {
      return 0;
  }
  
  bcf_hdr_t* cohort_hdr = bcf_hdr_read(cohort_vcf);
  bcf1_t* record = bcf_init1();
  
  size_t number_of_haplotypes = bcf_hdr_nsamples(cohort_hdr) * 2;
  size_t length = ref_end - ref_start;
  
  linearReferenceStructure* reference = linearReferenceStructure_init_empty(length);
  haplotypeCohort* cohort = haplotypeCohort_init_empty(number_of_haplotypes, reference);
  
  while(bcf_read(cohort_vcf, cohort_hdr, record) == 0) {
    bcf_unpack(record, BCF_UN_ALL);
    
    //TODO handle non-SNPs
    if (bcf_is_snp(record) != 0) {
      size_t site = record->pos;
      int64_t site_index = linearReferenceStructure_add_site(reference, site);
      haplotypeCohort_add_record(cohort, site);
      
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int ngt = bcf_get_genotypes(cohort_hdr, record, &gt_arr, &ngt_arr);
      // TODO: check that ngt is correct
      if(site_index >= 0) {
        for(size_t i = 0; i < ngt; i++) {
          int allele_index = bcf_gt_allele(gt_arr[i]);
          char allele_value = record->d.allele[allele_index][0];
          haplotypeCohort_set_sample_allele(cohort, i, site_index, allele_value);
        }
      } 
    }
  }
  
  linearReferenceStructure_calc_spans(reference, ref_end - ref_start);
  haplotypeCohort_populate_counts(cohort);

  *return_lr = reference;
  *return_cohort = cohort;
  
  vcf_close(cohort_vcf);
  return 1;
}
