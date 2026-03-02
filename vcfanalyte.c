/*
 * VcfAnalyte - Genetic Variant Analysis Tool
 * 
 * A C-based tool for parsing, annotating, filtering, and summarizing
 * genetic variants from VCF files. Supports VCF 4.2+ specification.
 * 
 * Features:
 * - VCF 4.2+ compliant parsing
 * - SNP, indel, CNV, SV detection
 * - Functional impact prediction
 * - Flexible filtering options
 * - Summary statistics
 * - Tab-delimited and JSON output
 * 
 * Compilation:
 *   gcc -O2 -o vcfanalyte vcfanalyte.c -lm
 *   or use provided Makefile
 * 
 * Usage:
 *   ./vcfanalyte input.vcf [options]
 *   cat input.vcf | ./vcfanalyte [options]
 * 
 * Author: Dany Mukesha
 * License: MIT
 */

#define _GNU_SOURCE
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>

/* ============================================================================
 * CONFIGURATION AND CONSTANTS
 * ============================================================================ */

#define MAX_LINE_LENGTH 100000
#define MAX_FIELDS 100
#define MAX_ALT_ALLELES 100
#define MAX_SAMPLES 1000
#define MAX_CHROM_LENGTH 50
#define MAX_ID_LENGTH 200
#define MAX_INFO_LENGTH 5000
#define MAX_GTF_FEATURES 50000

/* Variant type enumeration */
typedef enum {
    VAR_TYPE_SNP,           /* Single nucleotide polymorphism */
    VAR_TYPE_MNP,           /* Multi-nucleotide polymorphism */
    VAR_TYPE_INS,           /* Insertion */
    VAR_TYPE_DEL,           /* Deletion */
    VAR_TYPE_SV,            /* Structural variant */
    VAR_TYPE_CNV,           /* Copy number variant */
    VAR_TYPE_UNKNOWN        /* Unknown variant type */
} VariantType;

/* Functional impact classification */
typedef enum {
    IMPACT_SYNONYMOUS,      /* No amino acid change */
    IMPACT_MISSENSE,        /* Amino acid change */
    IMPACT_NONSENSE,        /* Premature stop codon */
    IMPACT_FRAMESHIFT,      /* Frameshift indel */
    IMPACT_IN_FRAME,        /* In-frame indel */
    IMPACT_SPLICE_SITE,     /* Splice region variant */
    IMPACT_STOP_LOST,       /* Stop codon lost */
    IMPACT_STOP_GAINED,     /* Stop codon gained (nonsense) */
    IMPACT_START_LOST,      /* Start codon lost */
    IMPACT_INTERGENIC,      /* Non-coding region */
    IMPACT_INTRON,          /* Intron variant */
    IMPACT_UTR5,            /* 5' UTR variant */
    IMPACT_UTR3,            /* 3' UTR variant */
    IMPACT_UNKNOWN          /* Unknown impact */
} FunctionalImpact;

/* Filter criteria structure */
typedef struct {
    double min_qual;                /* Minimum quality score */
    double max_qual;                /* Maximum quality score */
    int min_alt_freq;               /* Minimum alternate frequency (count) */
    double min_af;                  /* Minimum allele frequency */
    double max_af;                  /* Maximum allele frequency */
    int filter_snv;                 /* Include SNPs */
    int filter_indel;               /* Include indels */
    int filter_sv;                  /* Include SVs */
    int filter_cnv;                 /* Include CNVs */
    char chrom_filter[MAX_CHROM_LENGTH]; /* Chromosome to filter (empty = all) */
    int passed_only;                /* Only include variants that passed filters */
} FilterCriteria;

/* Sample genotype structure */
typedef struct {
    char *gt_string;                /* Genotype string (e.g., "0/1") */
    int *alleles;                   /* Array of allele indices */
    int num_alleles;                /* Number of alleles */
    double *phased;                 /* Phase information */
    int is_phased;                  /* Whether genotype is phased */
    int is_missing;                 /* Whether genotype is missing */
    double *gt_quals;               /* Genotype qualities (GQs) */
    double *ados;                   /* Allele depths */
} Genotype;

/* Variant record structure */
typedef struct {
    char chrom[MAX_CHROM_LENGTH];   /* Chromosome */
    int pos;                        /* Position (1-based) */
    int end;                        /* End position (for SVs) */
    char id[MAX_ID_LENGTH];         /* Variant ID (rs number) */
    char *ref;                      /* Reference allele */
    char *alt[MAX_ALT_ALLELES];     /* Alternate alleles */
    int num_alt;                    /* Number of alternate alleles */
    double qual;                    /* Quality score */
    char *filter;                   /* Filter status */
    char *info;                     /* INFO field string */
    char **info_keys;               /* Parsed INFO keys */
    char **info_values;             /* Parsed INFO values */
    int num_info;                   /* Number of INFO fields */
    
    /* Parsed variant properties */
    VariantType var_type;           /* Variant type */
    FunctionalImpact impact;        /* Functional impact */
    double *allele_freqs;           /* Allele frequencies */
    int *ac;                        /* Allele counts */
    int an;                         /* Total allele number */
    char sv_type[MAX_ID_LENGTH];    /* SV type (if applicable) */
    float cn;                       /* Copy number (for CNVs) */
    
    /* Genotype data */
    Genotype *samples;              /* Array of sample genotypes */
    int num_samples;                /* Number of samples */
    char **sample_names;            /* Sample names */
    
    /* Complexity scoring */
    float complexity_score;         /* Novel complexity metric */
    
    /* Annotation */
    char gene[MAX_ID_LENGTH];       /* Gene name (if annotated) */
    char consequence[MAX_ID_LENGTH]; /* Consequence type */
} Variant;

/* GTF/GFF feature for annotation */
typedef struct {
    char chrom[MAX_CHROM_LENGTH];
    char feature_type[MAX_ID_LENGTH];
    char gene_id[MAX_ID_LENGTH];
    char gene_name[MAX_ID_LENGTH];
    char transcript_id[MAX_ID_LENGTH];
    int start;
    int end;
    char strand;
} GtfFeature;

/* Statistics structure */
typedef struct {
    long total_variants;
    long snp_count;
    long mnp_count;
    long indel_count;
    long sv_count;
    long cnv_count;
    long unknown_count;
    
    long synonymous_count;
    long missense_count;
    long nonsense_count;
    long frameshift_count;
    long inframe_count;
    long splice_count;
    long intergenic_count;
    long intron_count;
    long utr_count;
    long unknown_impact_count;
    
    long passed_count;
    long filtered_count;
    
    double ti_tv_ratio;             /* Transition/transversion ratio */
    long ti_count;
    long tv_count;
    
    /* Per-chromosome counts */
    int chrom_counts[100];
    char chrom_names[100][MAX_CHROM_LENGTH];
    int num_chroms;
} VariantStats;

/* ============================================================================
 * MEMORY MANAGEMENT FUNCTIONS
 * ============================================================================ */

/*
 * Safe memory allocation with error checking
 */
void *xmalloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL && size > 0) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(1);
    }
    return ptr;
}

/*
 * Safe reallocation with error checking
 */
void *xrealloc(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (new_ptr == NULL && size > 0) {
        fprintf(stderr, "Error: Memory reallocation failed\n");
        exit(1);
    }
    return new_ptr;
}

/*
 * Safe string duplication
 */
char *xstrdup(const char *s) {
    if (s == NULL) return NULL;
    size_t len = strlen(s) + 1;
    char *copy = (char *)xmalloc(len);
    memcpy(copy, s, len);
    return copy;
}

/* ============================================================================
 * STRING UTILITY FUNCTIONS
 * ============================================================================ */

/*
 * Trim whitespace from both ends of a string
 */
char *trim(char *str) {
    if (str == NULL) return NULL;
    
    while (isspace((unsigned char)*str)) str++;
    
    if (*str == 0) return str;
    
    char *end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    
    *(end + 1) = '\0';
    return str;
}

/*
 * Count occurrences of a character in a string
 */
int count_char(const char *str, char c) {
    int count = 0;
    while (*str) {
        if (*str == c) count++;
        str++;
    }
    return count;
}

/*
 * Split string by delimiter into array of strings
 */
int split_string(const char *str, char delim, char ***result) {
    if (str == NULL || result == NULL) return 0;
    
    char *copy = xstrdup(str);
    int count = 1;
    
    for (char *p = copy; *p; p++) {
        if (*p == delim) count++;
    }
    
    *result = (char **)xmalloc(count * sizeof(char *));
    
    char *token = strtok(copy, (char[]){delim, 0});
    int i = 0;
    while (token != NULL) {
        (*result)[i++] = xstrdup(token);
        token = strtok(NULL, (char[]){delim, 0});
    }
    
    free(copy);
    return i;
}

/*
 * Free split string array
 */
void free_split_string(char **arr, int count) {
    if (arr == NULL) return;
    for (int i = 0; i < count; i++) {
        free(arr[i]);
    }
    free(arr);
}

/* ============================================================================
 * VCF PARSING FUNCTIONS
 * ============================================================================ */

/*
 * Check if line is a VCF header line
 */
int is_header_line(const char *line) {
    return line[0] == '#';
}

/*
 * Check if line is the VCF column header (starts with #CHROM)
 */
int is_column_header(const char *line) {
    return strncmp(line, "#CHROM", 6) == 0;
}

/*
 * Parse VCF header to extract sample names
 */
int parse_vcf_header(FILE *fp, char ***sample_names, int *num_samples) {
    char line[MAX_LINE_LENGTH];
    *num_samples = 0;
    *sample_names = NULL;
    
    while (fgets(line, sizeof(line), fp)) {
        trim(line);
        
        if (is_column_header(line)) {
            /* Parse sample names from header */
            char **fields = NULL;
            int num_fields = split_string(line, '\t', &fields);
            
            if (num_fields > 9) {
                *num_samples = num_fields - 9;
                *sample_names = (char **)xmalloc(*num_samples * sizeof(char *));
                for (int i = 0; i < *num_samples; i++) {
                    (*sample_names)[i] = xstrdup(fields[9 + i]);
                }
            }
            
            free_split_string(fields, num_fields);
            break;
        }
    }
    
    return 0;
}

/*
 * Parse INFO field string into key-value pairs
 */
int parse_info_field(const char *info_str, char ***keys, char ***values) {
    if (info_str == NULL || strlen(info_str) == 0) {
        *keys = NULL;
        *values = NULL;
        return 0;
    }
    
    char *info_copy = xstrdup(info_str);
    int count = count_char(info_copy, ';') + 1;
    
    *keys = (char **)xmalloc(count * sizeof(char *));
    *values = (char **)xmalloc(count * sizeof(char *));
    
    char *token = strtok(info_copy, ";");
    int i = 0;
    
    while (token != NULL && i < count) {
        char *eq = strchr(token, '=');
        if (eq) {
            *eq = '\0';
            (*keys)[i] = xstrdup(trim(token));
            (*values)[i] = xstrdup(trim(eq + 1));
        } else {
            (*keys)[i] = xstrdup(trim(token));
            (*values)[i] = xstrdup("");
        }
        token = strtok(NULL, ";");
        i++;
    }
    
    free(info_copy);
    return i;
}

/*
 * Get INFO field value by key
 */
const char *get_info_value(const Variant *var, const char *key) {
    for (int i = 0; i < var->num_info; i++) {
        if (strcmp(var->info_keys[i], key) == 0) {
            return var->info_values[i];
        }
    }
    return NULL;
}

/*
 * Parse genotype string from FORMAT field
 */
int parse_genotype(const char *gt_str, Genotype *gt) {
    /* Initialize all fields to NULL/0 to avoid garbage */
    gt->gt_string = NULL;
    gt->alleles = NULL;
    gt->num_alleles = 0;
    gt->phased = NULL;
    gt->is_phased = 0;
    gt->is_missing = 0;
    gt->gt_quals = NULL;
    gt->ados = NULL;
    
    if (gt_str == NULL || strcmp(gt_str, ".") == 0) {
        gt->is_missing = 1;
        gt->gt_string = xstrdup(".");
        gt->num_alleles = 0;
        gt->alleles = NULL;
        return 0;
    }
    
    gt->is_missing = 0;
    gt->gt_string = xstrdup(gt_str);
    
    /* Determine if phased (|) or unphased (/) */
    gt->is_phased = (strchr(gt_str, '|') != NULL);
    char delim = gt->is_phased ? '|' : '/';
    
    /* Count alleles */
    gt->num_alleles = count_char(gt_str, delim) + 1;
    gt->alleles = (int *)xmalloc(gt->num_alleles * sizeof(int));
    
    /* Parse allele indices */
    char *copy = xstrdup(gt_str);
    char *token = strtok(copy, (char[]){delim, 0});
    int i = 0;
    
    while (token != NULL && i < gt->num_alleles) {
        if (strcmp(token, ".") == 0) {
            gt->alleles[i] = -1;  /* Missing allele */
        } else {
            gt->alleles[i] = atoi(token);
        }
        token = strtok(NULL, (char[]){delim, 0});
        i++;
    }
    
    free(copy);
    return 0;
}

/*
 * Parse a single VCF record line into a Variant structure
 */
int parse_vcf_record(const char *line, Variant *var, char **sample_names, int num_samples) {
    if (line == NULL || strlen(line) == 0) return -1;
    
    char **fields = NULL;
    int num_fields = split_string(line, '\t', &fields);
    
    if (num_fields < 8) {
        fprintf(stderr, "Error: VCF record has fewer than 8 fields\n");
        free_split_string(fields, num_fields);
        return -1;
    }
    
    /* Parse mandatory fields */
    strncpy(var->chrom, fields[0], MAX_CHROM_LENGTH - 1);
    var->chrom[MAX_CHROM_LENGTH - 1] = '\0';
    var->pos = atoi(fields[1]);
    strncpy(var->id, fields[2], MAX_ID_LENGTH - 1);
    var->id[MAX_ID_LENGTH - 1] = '\0';
    
    var->ref = xstrdup(fields[3]);
    
    /* Parse alternate alleles */
    {
        char **alt_arr = NULL;
        var->num_alt = split_string(fields[4], ',', &alt_arr);
        for (int i = 0; i < var->num_alt && i < MAX_ALT_ALLELES; i++) {
            var->alt[i] = xstrdup(alt_arr[i]);
        }
        free_split_string(alt_arr, var->num_alt);
    }
    
    /* Quality score */
    if (strcmp(fields[5], ".") == 0) {
        var->qual = 0.0;
    } else {
        var->qual = atof(fields[5]);
    }
    
    /* Filter status */
    var->filter = xstrdup(fields[6]);
    
    /* INFO field */
    var->info = xstrdup(fields[7]);
    var->num_info = parse_info_field(var->info, &var->info_keys, &var->info_values);
    
    /* Parse additional INFO fields for SV/CNV */
    const char *svtype = get_info_value(var, "SVTYPE");
    if (svtype) {
        strncpy(var->sv_type, svtype, MAX_ID_LENGTH - 1);
    } else {
        var->sv_type[0] = '\0';
    }
    
    const char *cn_str = get_info_value(var, "CN");
    if (cn_str) {
        var->cn = atof(cn_str);
    } else {
        var->cn = 2.0;  /* Default diploid */
    }
    
    /* END position for SVs */
    const char *end_str = get_info_value(var, "END");
    if (end_str) {
        var->end = atoi(end_str);
    } else {
        var->end = var->pos + strlen(var->ref) - 1;
    }
    
    /* Parse allele frequencies from INFO */
    const char *af_str = get_info_value(var, "AF");
    var->allele_freqs = NULL;
    if (af_str && strcmp(af_str, ".") != 0) {
        var->allele_freqs = (double *)xmalloc(var->num_alt * sizeof(double));
        char *af_copy = xstrdup(af_str);
        char *token = strtok(af_copy, ",");
        int i = 0;
        while (token != NULL && i < var->num_alt) {
            var->allele_freqs[i++] = atof(token);
            token = strtok(NULL, ",");
        }
        free(af_copy);
    }
    
    /* Parse allele counts */
    const char *ac_str = get_info_value(var, "AC");
    var->ac = NULL;
    if (ac_str && strcmp(ac_str, ".") != 0) {
        var->ac = (int *)xmalloc(var->num_alt * sizeof(int));
        char *ac_copy = xstrdup(ac_str);
        char *token = strtok(ac_copy, ",");
        int i = 0;
        while (token != NULL && i < var->num_alt) {
            var->ac[i++] = atoi(token);
            token = strtok(NULL, ",");
        }
        free(ac_copy);
    }
    
    /* AN (total allele number) */
    const char *an_str = get_info_value(var, "AN");
    var->an = an_str ? atoi(an_str) : 0;
    
    /* Parse sample genotypes if present */
    var->num_samples = 0;
    var->samples = NULL;
    var->sample_names = NULL;
    
    if (num_samples > 0 && num_fields > 8) {
        var->num_samples = num_samples;
        var->samples = (Genotype *)xmalloc(num_samples * sizeof(Genotype));
        var->sample_names = (char **)xmalloc(num_samples * sizeof(char *));
        
        /* Parse FORMAT field if present (reserved for future use) */
        if (num_fields > 8 && strcmp(fields[8], ".") != 0) {
            /* FORMAT field: fields[8] */
        }
        
        for (int i = 0; i < num_samples; i++) {
            var->sample_names[i] = xstrdup(sample_names[i]);
            
            if (num_fields > 9 + i && strcmp(fields[9 + i], ".") != 0) {
                parse_genotype(fields[9 + i], &var->samples[i]);
            } else {
                parse_genotype(".", &var->samples[i]);
            }
        }
    }
    
    /* Initialize gene and consequence */
    var->gene[0] = '\0';
    var->consequence[0] = '\0';
    
    /* Try to get gene annotation from INFO */
    const char *gene = get_info_value(var, "GENE");
    if (gene) {
        strncpy(var->gene, gene, MAX_ID_LENGTH - 1);
    }
    
    const char *conseq = get_info_value(var, "CSQ");
    if (conseq) {
        strncpy(var->consequence, conseq, MAX_ID_LENGTH - 1);
    } else {
        const char *conseq_alt = get_info_value(var, "ANN");
        if (conseq_alt) {
            strncpy(var->consequence, conseq_alt, MAX_ID_LENGTH - 1);
        }
    }
    
    free_split_string(fields, num_fields);
    return 0;
}

/* ============================================================================
 * VARIANT TYPE AND IMPACT PREDICTION
 * ============================================================================ */

/*
 * Determine variant type based on ref and alt alleles
 * 
 * Biological rationale:
 * - SNPs are single nucleotide changes
 * - MNPs are multiple nucleotide changes of same length
 * - Insertions have alt longer than ref
 * - Deletions have alt shorter than ref
 * - SVs are identified by SVTYPE INFO field or large variants
 * - CNVs have copy number information
 */
VariantType determine_variant_type(const Variant *var) {
    if (var->num_alt == 0) return VAR_TYPE_UNKNOWN;
    
    const char *ref = var->ref;
    const char *alt = var->alt[0];  /* Use first alt for type determination */
    
    int ref_len = strlen(ref);
    int alt_len = strlen(alt);
    
    /* Check for SV/CNV from INFO field first */
    if (strlen(var->sv_type) > 0) {
        if (strcmp(var->sv_type, "CNV") == 0 || strcmp(var->sv_type, "COPY_NUMBER") == 0) {
            return VAR_TYPE_CNV;
        }
        return VAR_TYPE_SV;
    }
    
    /* Check for large structural variants (> 50bp) */
    if (abs(alt_len - ref_len) > 50 || var->end - var->pos > 50) {
        if (alt_len > ref_len) return VAR_TYPE_INS;
        if (alt_len < ref_len) return VAR_TYPE_DEL;
        return VAR_TYPE_SV;
    }
    
    /* Small variants */
    if (ref_len == 1 && alt_len == 1) {
        return VAR_TYPE_SNP;  /* Single nucleotide substitution */
    }
    
    if (ref_len == alt_len) {
        return VAR_TYPE_MNP;  /* Multi-nucleotide polymorphism */
    }
    
    if (alt_len > ref_len) {
        return VAR_TYPE_INS;  /* Insertion */
    }
    
    return VAR_TYPE_DEL;  /* Deletion */
}

/*
 * Determine functional impact based on variant type and position
 * 
 * This is a simplified impact prediction based on:
 * - Variant type (SNP, indel)
 * - Position relative to gene annotations (would need GTF for accuracy)
 * - For known genes, we classify based on codon changes
 * 
 * In production, this would use transcript sequences and variant positions
 */
FunctionalImpact predict_functional_impact(const Variant *var) {
    VariantType vtype = var->var_type;
    
    /* If we have gene annotation, use it */
    if (strlen(var->gene) > 0) {
        if (strlen(var->consequence) > 0) {
            if (strstr(var->consequence, "synonymous")) return IMPACT_SYNONYMOUS;
            if (strstr(var->consequence, "missense")) return IMPACT_MISSENSE;
            if (strstr(var->consequence, "nonsense") || strstr(var->consequence, "stop_gained")) {
                return IMPACT_NONSENSE;
            }
            if (strstr(var->consequence, "frameshift")) return IMPACT_FRAMESHIFT;
            if (strstr(var->consequence, "inframe")) return IMPACT_IN_FRAME;
            if (strstr(var->consequence, "splice")) return IMPACT_SPLICE_SITE;
            if (strstr(var->consequence, "stop_lost")) return IMPACT_STOP_LOST;
            if (strstr(var->consequence, "start_lost")) return IMPACT_START_LOST;
            if (strstr(var->consequence, "UTR")) return IMPACT_UTR5;
            if (strstr(var->consequence, "intron")) return IMPACT_INTRON;
            if (strstr(var->consequence, "intergenic")) return IMPACT_INTERGENIC;
        }
    }
    
    /* Simplified classification based on variant type */
    switch (vtype) {
        case VAR_TYPE_SNP:
            /* For SNPs in coding regions without annotation, 
               check if they cause amino acid change */
            {
                const char *ref = var->ref;
                const char *alt = var->alt[0];
                if (strlen(ref) == 1 && strlen(alt) == 1) {
                    /* Simple heuristic: check if it's a transition or transversion
                       and whether it affects coding sequence */
                    return IMPACT_MISSENSE;  /* Conservative assumption */
                }
            }
            return IMPACT_UNKNOWN;
            
        case VAR_TYPE_MNP:
            return IMPACT_MISSENSE;  /* Usually missense */
            
        case VAR_TYPE_INS:
        case VAR_TYPE_DEL:
            {
                int ref_len = strlen(var->ref);
                int alt_len = strlen(var->alt[0]);
                int len_diff = abs(alt_len - ref_len);
                
                /* Frameshift if indel length not divisible by 3 */
                if (len_diff % 3 != 0) {
                    return IMPACT_FRAMESHIFT;
                } else {
                    return IMPACT_IN_FRAME;
                }
            }
            
        case VAR_TYPE_SV:
            return IMPACT_INTERGENIC;  /* Typically large rearrangements */
            
        case VAR_TYPE_CNV:
            return IMPACT_INTERGENIC;  /* Copy number changes */
            
        default:
            return IMPACT_UNKNOWN;
    }
}

/*
 * Compute variant complexity score
 * 
 * Novel feature: This scoring system helps identify complex genomic regions
 * that may warrant closer investigation. Higher scores indicate more complex
 * variants that could be difficult to genotype or may have regulatory effects.
 * 
 * Factors:
 * - Number of alternate alleles (multi-allelic = more complex)
 * - Length of indel/SV (longer = more complex)
 * - Presence of structural variant annotations
 * - Number of samples with the variant
 * - Allele frequency (rare variants score higher)
 */
float compute_complexity_score(const Variant *var) {
    float score = 0.0;
    
    /* Base complexity from alternate allele count */
    score += (var->num_alt - 1) * 2.0;  /* Each additional alt adds 2 points */
    
    /* Indel/SV length complexity */
    if (var->num_alt > 0) {
        int ref_len = strlen(var->ref);
        int alt_len = strlen(var->alt[0]);
        int len_diff = abs(alt_len - ref_len);
        
        if (len_diff > 0) {
            /* Logarithmic scaling for indel length */
            score += log1p(len_diff) * 1.5;
        }
    }
    
    /* SV/CNV complexity */
    if (var->var_type == VAR_TYPE_SV || var->var_type == VAR_TYPE_CNV) {
        score += 5.0;  /* Base score for structural variants */
        
        /* Additional complexity for SV type */
        if (strcmp(var->sv_type, "INVERSION") == 0) score += 3.0;
        if (strcmp(var->sv_type, "TRANSLOCATION") == 0) score += 4.0;
        if (var->var_type == VAR_TYPE_CNV) {
            /* CNV complexity based on copy number deviation from diploid */
            score += fabs(var->cn - 2.0) * 2.0;
        }
    }
    
    /* Sample heterogeneity complexity */
    if (var->num_samples > 0) {
        int het_count = 0;
        for (int i = 0; i < var->num_samples; i++) {
            Genotype *gt = &var->samples[i];
            if (!gt->is_missing && gt->num_alleles >= 2) {
                if (gt->alleles[0] != gt->alleles[1]) {
                    het_count++;
                }
            }
        }
        /* Higher heterozygosity across samples adds complexity */
        score += (het_count / (float)var->num_samples) * 3.0;
    }
    
    /* Rarity bonus (rare variants may be more complex) */
    if (var->allele_freqs && var->allele_freqs[0] > 0) {
        if (var->allele_freqs[0] < 0.01) {
            score += 2.0;  /* Rare variant bonus */
        } else if (var->allele_freqs[0] < 0.05) {
            score += 1.0;
        }
    }
    
    return score;
}

/*
 * Calculate transition/transversion ratio for SNPs
 * 
 * Transitions (Ti): A<->G, C<->T
 * Transversions (Tv): A<->C, A<->T, G<->C, G<->T
 * 
 * This ratio is important in population genetics as it typically
 * exceeds 2.0 in most genomes due to methylation-mediated deamination
 */
void classify_snp(const Variant *var, int *is_transition) {
    if (var->var_type != VAR_TYPE_SNP || var->num_alt < 1) {
        *is_transition = -1;
        return;
    }
    
    char ref = toupper(var->ref[0]);
    char alt = toupper(var->alt[0][0]);
    
    /* Transitions */
    if ((ref == 'A' && alt == 'G') || (ref == 'G' && alt == 'A') ||
        (ref == 'C' && alt == 'T') || (ref == 'T' && alt == 'C')) {
        *is_transition = 1;
        return;
    }
    
    /* Transversions */
    *is_transition = 0;
}

/* ============================================================================
 * FILTERING FUNCTIONS
 * ============================================================================ */

/*
 * Check if variant passes filter criteria
 */
int variant_passes_filter(const Variant *var, const FilterCriteria *filter) {
    /* Quality filter */
    if (var->qual < filter->min_qual) return 0;
    if (filter->max_qual > 0 && var->qual > filter->max_qual) return 0;
    
    /* Passed filter only */
    if (filter->passed_only && strcmp(var->filter, "PASS") != 0) return 0;
    
    /* Chromosome filter */
    if (strlen(filter->chrom_filter) > 0 && 
        strcmp(var->chrom, filter->chrom_filter) != 0) {
        return 0;
    }
    
    /* Variant type filter */
    switch (var->var_type) {
        case VAR_TYPE_SNP:
        case VAR_TYPE_MNP:
            if (!filter->filter_snv) return 0;
            break;
        case VAR_TYPE_INS:
        case VAR_TYPE_DEL:
            if (!filter->filter_indel) return 0;
            break;
        case VAR_TYPE_SV:
            if (!filter->filter_sv) return 0;
            break;
        case VAR_TYPE_CNV:
            if (!filter->filter_cnv) return 0;
            break;
        default:
            break;
    }
    
    /* Allele frequency filter */
    if (var->allele_freqs && var->allele_freqs[0] >= 0) {
        if (filter->min_af > 0 && var->allele_freqs[0] < filter->min_af) return 0;
        if (filter->max_af < 1.0 && var->allele_freqs[0] > filter->max_af) return 0;
    }
    
    return 1;
}

/* ============================================================================
 * STATISTICS FUNCTIONS
 * ============================================================================ */

/*
 * Initialize statistics structure
 */
void init_stats(VariantStats *stats) {
    memset(stats, 0, sizeof(VariantStats));
    stats->ti_tv_ratio = 0.0;
    stats->num_chroms = 0;
}

/*
 * Update statistics with a variant
 */
void update_stats(VariantStats *stats, const Variant *var, int passed) {
    stats->total_variants++;
    
    if (passed) {
        stats->passed_count++;
        
        /* Count by variant type */
        switch (var->var_type) {
            case VAR_TYPE_SNP:
                stats->snp_count++;
                {
                    int is_transition;
                    classify_snp(var, &is_transition);
                    if (is_transition == 1) {
                        stats->ti_count++;
                    } else if (is_transition == 0) {
                        stats->tv_count++;
                    }
                }
                break;
            case VAR_TYPE_MNP:
                stats->mnp_count++;
                break;
            case VAR_TYPE_INS:
            case VAR_TYPE_DEL:
                stats->indel_count++;
                break;
            case VAR_TYPE_SV:
                stats->sv_count++;
                break;
            case VAR_TYPE_CNV:
                stats->cnv_count++;
                break;
            default:
                stats->unknown_count++;
                break;
        }
        
        /* Count by impact */
        switch (var->impact) {
            case IMPACT_SYNONYMOUS:
                stats->synonymous_count++;
                break;
            case IMPACT_MISSENSE:
                stats->missense_count++;
                break;
            case IMPACT_NONSENSE:
            case IMPACT_STOP_GAINED:
                stats->nonsense_count++;
                break;
            case IMPACT_FRAMESHIFT:
                stats->frameshift_count++;
                break;
            case IMPACT_IN_FRAME:
                stats->inframe_count++;
                break;
            case IMPACT_SPLICE_SITE:
                stats->splice_count++;
                break;
            case IMPACT_INTERGENIC:
                stats->intergenic_count++;
                break;
            case IMPACT_INTRON:
                stats->intron_count++;
                break;
            case IMPACT_UTR5:
            case IMPACT_UTR3:
                stats->utr_count++;
                break;
            default:
                stats->unknown_impact_count++;
                break;
        }
    } else {
        stats->filtered_count++;
    }
    
    /* Chromosome counts */
    int found = 0;
    for (int i = 0; i < stats->num_chroms; i++) {
        if (strcmp(stats->chrom_names[i], var->chrom) == 0) {
            stats->chrom_counts[i]++;
            found = 1;
            break;
        }
    }
    if (!found && stats->num_chroms < 100) {
        strncpy(stats->chrom_names[stats->num_chroms], var->chrom, MAX_CHROM_LENGTH - 1);
        stats->chrom_counts[stats->num_chroms] = 1;
        stats->num_chroms++;
    }
    
    /* Calculate Ti/Tv ratio */
    if (stats->tv_count > 0) {
        stats->ti_tv_ratio = (double)stats->ti_count / stats->tv_count;
    }
}

/*
 * Print statistics summary
 */
void print_stats(const VariantStats *stats, FILE *out) {
    fprintf(out, "\n========== VARIANT STATISTICS ==========\n\n");
    fprintf(out, "Total variants processed: %ld\n", stats->total_variants);
    fprintf(out, "Passed filters:           %ld\n", stats->passed_count);
    fprintf(out, "Filtered out:             %ld\n\n", stats->filtered_count);
    
    fprintf(out, "--- Variant Types ---\n");
    fprintf(out, "SNPs:            %ld\n", stats->snp_count);
    fprintf(out, "MNPs:            %ld\n", stats->mnp_count);
    fprintf(out, "Indels:          %ld\n", stats->indel_count);
    fprintf(out, "SVs:             %ld\n", stats->sv_count);
    fprintf(out, "CNVs:            %ld\n", stats->cnv_count);
    fprintf(out, "Unknown:         %ld\n\n", stats->unknown_count);
    
    fprintf(out, "--- Functional Impact ---\n");
    fprintf(out, "Synonymous:      %ld\n", stats->synonymous_count);
    fprintf(out, "Missense:        %ld\n", stats->missense_count);
    fprintf(out, "Nonsense:        %ld\n", stats->nonsense_count);
    fprintf(out, "Frameshift:      %ld\n", stats->frameshift_count);
    fprintf(out, "In-frame:        %ld\n", stats->inframe_count);
    fprintf(out, "Splice site:     %ld\n", stats->splice_count);
    fprintf(out, "Intergenic:      %ld\n", stats->intergenic_count);
    fprintf(out, "Intron:          %ld\n", stats->intron_count);
    fprintf(out, "UTR:             %ld\n", stats->utr_count);
    fprintf(out, "Unknown:         %ld\n\n", stats->unknown_impact_count);
    
    fprintf(out, "--- SNP Statistics ---\n");
    fprintf(out, "Transitions (Ti): %ld\n", stats->ti_count);
    fprintf(out, "Transversions(Tv): %ld\n", stats->tv_count);
    fprintf(out, "Ti/Tv ratio:     %.4f\n\n", stats->ti_tv_ratio);
    
    fprintf(out, "--- Per-Chromosome Counts ---\n");
    for (int i = 0; i < stats->num_chroms; i++) {
        fprintf(out, "chr%s: %d\n", stats->chrom_names[i], stats->chrom_counts[i]);
    }
    fprintf(out, "===========================================\n");
}

/* ============================================================================
 * OUTPUT FUNCTIONS
 * ============================================================================ */

/*
 * Get string representation of variant type
 */
const char *variant_type_str(VariantType vt) {
    switch (vt) {
        case VAR_TYPE_SNP: return "SNP";
        case VAR_TYPE_MNP: return "MNP";
        case VAR_TYPE_INS: return "INSERTION";
        case VAR_TYPE_DEL: return "DELETION";
        case VAR_TYPE_SV: return "SV";
        case VAR_TYPE_CNV: return "CNV";
        default: return "UNKNOWN";
    }
}

/*
 * Get string representation of functional impact
 */
const char *impact_str(FunctionalImpact impact) {
    switch (impact) {
        case IMPACT_SYNONYMOUS: return "synonymous";
        case IMPACT_MISSENSE: return "missense";
        case IMPACT_NONSENSE: return "nonsense";
        case IMPACT_FRAMESHIFT: return "frameshift";
        case IMPACT_IN_FRAME: return "in_frame";
        case IMPACT_SPLICE_SITE: return "splice_site";
        case IMPACT_STOP_LOST: return "stop_lost";
        case IMPACT_STOP_GAINED: return "stop_gained";
        case IMPACT_START_LOST: return "start_lost";
        case IMPACT_INTERGENIC: return "intergenic";
        case IMPACT_INTRON: return "intron";
        case IMPACT_UTR5: return "UTR5";
        case IMPACT_UTR3: return "UTR3";
        default: return "unknown";
    }
}

/*
 * Print variant in tab-delimited format
 */
void print_variant_tsv(const Variant *var, FILE *out) {
    /* Build alt allele string */
    char alt_str[MAX_INFO_LENGTH] = {0};
    for (int i = 0; i < var->num_alt; i++) {
        if (i > 0) strcat(alt_str, ",");
        strncat(alt_str, var->alt[i], MAX_INFO_LENGTH - strlen(alt_str) - 1);
    }
    
    /* Build allele frequency string */
    char af_str[500] = {0};
    if (var->allele_freqs) {
        for (int i = 0; i < var->num_alt; i++) {
            if (i > 0) strcat(af_str, ",");
            char tmp[50];
            snprintf(tmp, sizeof(tmp), "%.4f", var->allele_freqs[i]);
            strncat(af_str, tmp, 500 - strlen(af_str) - 1);
        }
    } else {
        strcpy(af_str, ".");
    }
    
    fprintf(out, "%s\t%d\t%s\t%s\t%s\t%.1f\t%s\t%s\t%.2f\t%s\t%s\n",
            var->chrom,
            var->pos,
            var->id,
            var->ref,
            alt_str,
            var->qual,
            var->filter,
            variant_type_str(var->var_type),
            var->complexity_score,
            impact_str(var->impact),
            af_str);
}

/*
 * Print header for TSV output
 */
void print_tsv_header(FILE *out) {
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tTYPE\tCOMPLEXITY\tIMPACT\tAF\n");
}

/*
 * Print variant in JSON format
 */
void print_variant_json(const Variant *var, FILE *out, int is_first) {
    if (!is_first) fprintf(out, ",\n");
    
    fprintf(out, "  {\n");
    fprintf(out, "    \"chrom\": \"%s\",\n", var->chrom);
    fprintf(out, "    \"pos\": %d,\n", var->pos);
    fprintf(out, "    \"end\": %d,\n", var->end);
    fprintf(out, "    \"id\": \"%s\",\n", var->id);
    fprintf(out, "    \"ref\": \"%s\",\n", var->ref);
    
    fprintf(out, "    \"alt\": [");
    for (int i = 0; i < var->num_alt; i++) {
        if (i > 0) fprintf(out, ", ");
        fprintf(out, "\"%s\"", var->alt[i]);
    }
    fprintf(out, "],\n");
    
    fprintf(out, "    \"qual\": %.1f,\n", var->qual);
    fprintf(out, "    \"filter\": \"%s\",\n", var->filter);
    fprintf(out, "    \"variant_type\": \"%s\",\n", variant_type_str(var->var_type));
    fprintf(out, "    \"functional_impact\": \"%s\",\n", impact_str(var->impact));
    fprintf(out, "    \"complexity_score\": %.4f,\n", var->complexity_score);
    
    if (strlen(var->sv_type) > 0) {
        fprintf(out, "    \"sv_type\": \"%s\",\n", var->sv_type);
    }
    if (var->cn != 2.0) {
        fprintf(out, "    \"copy_number\": %.1f,\n", var->cn);
    }
    if (strlen(var->gene) > 0) {
        fprintf(out, "    \"gene\": \"%s\",\n", var->gene);
    }
    if (strlen(var->consequence) > 0) {
        fprintf(out, "    \"consequence\": \"%s\",\n", var->consequence);
    }
    
    fprintf(out, "    \"allele_frequencies\": [");
    if (var->allele_freqs) {
        for (int i = 0; i < var->num_alt; i++) {
            if (i > 0) fprintf(out, ", ");
            fprintf(out, "%.4f", var->allele_freqs[i]);
        }
    }
    fprintf(out, "],\n");
    
    fprintf(out, "    \"num_alt_alleles\": %d,\n", var->num_alt);
    fprintf(out, "    \"num_samples\": %d\n", var->num_samples);
    
    fprintf(out, "  }");
}

/* ============================================================================
 * VARIANT STORAGE AND CLEANUP
 * ============================================================================ */

/*
 * Free variant structure
 */
void free_variant(Variant *var) {
    if (var == NULL) return;
    
    free(var->ref);
    
    for (int i = 0; i < var->num_alt; i++) {
        free(var->alt[i]);
    }
    
    free(var->filter);
    free(var->info);
    
    for (int i = 0; i < var->num_info; i++) {
        free(var->info_keys[i]);
        free(var->info_values[i]);
    }
    free(var->info_keys);
    free(var->info_values);
    
    if (var->allele_freqs) free(var->allele_freqs);
    if (var->ac) free(var->ac);
    
    /* Free sample genotypes */
    if (var->samples) {
        for (int i = 0; i < var->num_samples; i++) {
            if (var->samples[i].gt_string) free(var->samples[i].gt_string);
            if (var->samples[i].alleles) free(var->samples[i].alleles);
            if (var->samples[i].gt_quals) free(var->samples[i].gt_quals);
            if (var->samples[i].ados) free(var->samples[i].ados);
        }
        free(var->samples);
    }
    
    if (var->sample_names) {
        for (int i = 0; i < var->num_samples; i++) {
            if (var->sample_names[i]) free(var->sample_names[i]);
        }
        free(var->sample_names);
    }
}

/* ============================================================================
 * MAIN PROGRAM
 * ============================================================================ */

/*
 * Print usage information
 */
void print_usage(const char *prog) {
    fprintf(stderr, "VcfAnalyte - Genetic Variant Analysis Tool\n\n");
    fprintf(stderr, "Usage: %s [options] [input.vcf]\n\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help              Show this help message\n");
    fprintf(stderr, "  -o, --output FILE       Output file (default: stdout)\n");
    fprintf(stderr, "  -j, --json               Output in JSON format\n");
    fprintf(stderr, "  -s, --stats              Print summary statistics\n");
    fprintf(stderr, "  -H, --no-header          Don't print output header\n");
    fprintf(stderr, "\nFilter Options:\n");
    fprintf(stderr, "  -q, --min-qual FLOAT     Minimum quality score (default: 0)\n");
    fprintf(stderr, "  -Q, --max-qual FLOAT     Maximum quality score (default: no limit)\n");
    fprintf(stderr, "  -t, --type TYPE          Variant type to include (SNP, INDEL, SV, CNV, all)\n");
    fprintf(stderr, "  -f, --min-af FLOAT       Minimum allele frequency (default: 0)\n");
    fprintf(stderr, "  -F, --max-af FLOAT       Maximum allele frequency (default: 1.0)\n");
    fprintf(stderr, "  -c, --chrom CHROM        Filter by chromosome\n");
    fprintf(stderr, "  -p, --passed-only        Only include variants that passed filters\n");
    fprintf(stderr, "\nInput Options:\n");
    fprintf(stderr, "  -i, --stdin              Read from stdin\n");
    fprintf(stderr, "\nExample:\n");
    fprintf(stderr, "  %s input.vcf -s -o results.txt\n", prog);
    fprintf(stderr, "  cat input.vcf | %s -s --json > results.json\n", prog);
    fprintf(stderr, "  %s input.vcf -t SNP -q 30 -f 0.01\n", prog);
}

/*
 * Parse command line arguments
 */
int parse_arguments(int argc, char **argv, 
                   char **input_file, 
                   char **output_file,
                   int *output_json,
                   int *print_stats_flag,
                   int *print_header,
                   FilterCriteria *filter) {
    
    /* Set defaults */
    *input_file = NULL;
    *output_file = NULL;
    *output_json = 0;
    *print_stats_flag = 0;
    *print_header = 1;
    
    /* Default filter criteria */
    filter->min_qual = 0.0;
    filter->max_qual = 0.0;
    filter->min_alt_freq = 0;
    filter->min_af = 0.0;
    filter->max_af = 1.0;
    filter->filter_snv = 1;
    filter->filter_indel = 1;
    filter->filter_sv = 1;
    filter->filter_cnv = 1;
    filter->chrom_filter[0] = '\0';
    filter->passed_only = 0;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            *output_file = argv[++i];
        }
        else if (strcmp(argv[i], "-j") == 0 || strcmp(argv[i], "--json") == 0) {
            *output_json = 1;
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--stats") == 0) {
            *print_stats_flag = 1;
        }
        else if (strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "--no-header") == 0) {
            *print_header = 0;
        }
        else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--min-qual") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            filter->min_qual = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-Q") == 0 || strcmp(argv[i], "--max-qual") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            filter->max_qual = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            char *type = argv[++i];
            filter->filter_snv = 0;
            filter->filter_indel = 0;
            filter->filter_sv = 0;
            filter->filter_cnv = 0;
            
            if (strcmp(type, "SNP") == 0) filter->filter_snv = 1;
            else if (strcmp(type, "INDEL") == 0) filter->filter_indel = 1;
            else if (strcmp(type, "SV") == 0) filter->filter_sv = 1;
            else if (strcmp(type, "CNV") == 0) filter->filter_cnv = 1;
            else if (strcmp(type, "all") == 0) {
                filter->filter_snv = 1;
                filter->filter_indel = 1;
                filter->filter_sv = 1;
                filter->filter_cnv = 1;
            }
        }
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--min-af") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            filter->min_af = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-F") == 0 || strcmp(argv[i], "--max-af") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            filter->max_af = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--chrom") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            strncpy(filter->chrom_filter, argv[++i], MAX_CHROM_LENGTH - 1);
        }
        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--passed-only") == 0) {
            filter->passed_only = 1;
        }
        else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--stdin") == 0) {
            *input_file = "-";  /* Special value for stdin */
        }
        else if (argv[i][0] != '-') {
            *input_file = argv[i];
        }
        else {
            fprintf(stderr, "Error: Unknown option: %s\n", argv[i]);
            return -1;
        }
    }
    
    return 0;
}

/*
 * Main program entry point
 */
int main(int argc, char **argv) {
    /* Set stdout to unbuffered mode for immediate output */
    setbuf(stdout, NULL);
    
    char *input_file = NULL;
    char *output_file = NULL;
    int output_json = 0;
    int print_stats_flag = 0;
    int print_header = 1;
    FilterCriteria filter;
    
    /* Parse command line arguments */
    if (parse_arguments(argc, argv, &input_file, &output_file, 
                       &output_json, &print_stats_flag, &print_header,
                       &filter) != 0) {
        print_usage(argv[0]);
        return 1;
    }
    
    /* Open input file */
    FILE *fin = stdin;
    int close_fin = 0;
    
    if (input_file != NULL && strcmp(input_file, "-") != 0) {
        fin = fopen(input_file, "r");
        if (fin == NULL) {
            fprintf(stderr, "Error: Cannot open input file '%s': %s\n", 
                    input_file, strerror(errno));
            return 1;
        }
        close_fin = 1;
    }
    
    /* Open output file */
    FILE *fout = stdout;
    int close_fout = 0;
    
    if (output_file != NULL) {
        fout = fopen(output_file, "w");
        if (fout == NULL) {
            fprintf(stderr, "Error: Cannot open output file '%s': %s\n",
                    output_file, strerror(errno));
            if (close_fin) fclose(fin);
            return 1;
        }
        close_fout = 1;
    }
    
    /* Parse VCF header to get sample names */
    char **sample_names = NULL;
    int num_samples = 0;
    
    rewind(fin);
    parse_vcf_header(fin, &sample_names, &num_samples);
    
    /* Print JSON opening if needed */
    if (output_json) {
        fprintf(fout, "{\n");
        fprintf(fout, "  \"variants\": [\n");
    }
    
    /* Process VCF records */
    char line[MAX_LINE_LENGTH];
    VariantStats stats;
    init_stats(&stats);
    
    int first_variant = 1;
    int variant_count = 0;
    
    while (fgets(line, sizeof(line), fin)) {
        /* Skip header lines */
        if (is_header_line(line) && !is_column_header(line)) {
            continue;
        }
        
        /* Skip column header */
        if (is_column_header(line)) {
            continue;
        }
        
        /* Trim and skip empty lines */
        trim(line);
        if (strlen(line) == 0) continue;
        
        /* Parse variant record */
        Variant var;
        memset(&var, 0, sizeof(Variant));
        
        if (parse_vcf_record(line, &var, sample_names, num_samples) != 0) {
            continue;
        }
        
        /* Determine variant type and impact */
        var.var_type = determine_variant_type(&var);
        var.impact = predict_functional_impact(&var);
        var.complexity_score = compute_complexity_score(&var);
        
        /* Check if variant passes filter */
        int passed = variant_passes_filter(&var, &filter);
        
        /* Update statistics */
        update_stats(&stats, &var, passed);
        
        /* Output variant if it passes filter */
        if (passed) {
            if (output_json) {
                print_variant_json(&var, fout, first_variant);
                first_variant = 0;
            } else {
                if (print_header) {
                    print_tsv_header(fout);
                    print_header = 0;
                }
                print_variant_tsv(&var, fout);
            }
            fflush(fout);
            variant_count++;
        }
        
        /* Clean up variant */
        free_variant(&var);
    }
    
    /* Print JSON closing */
    if (output_json) {
        fprintf(fout, "\n  ],\n");
        fprintf(fout, "  \"total_variants\": %d,\n", variant_count);
        fprintf(fout, "  \"statistics\": {\n");
        fprintf(fout, "    \"snp_count\": %ld,\n", stats.snp_count);
        fprintf(fout, "    \"indel_count\": %ld,\n", stats.indel_count);
        fprintf(fout, "    \"sv_count\": %ld,\n", stats.sv_count);
        fprintf(fout, "    \"cnv_count\": %ld,\n", stats.cnv_count);
        fprintf(fout, "    \"ti_tv_ratio\": %.4f,\n", stats.ti_tv_ratio);
        fprintf(fout, "    \"synonymous\": %ld,\n", stats.synonymous_count);
        fprintf(fout, "    \"missense\": %ld,\n", stats.missense_count);
        fprintf(fout, "    \"nonsense\": %ld\n", stats.nonsense_count);
        fprintf(fout, "  }\n");
        fprintf(fout, "}\n");
    }
    
    /* Print statistics if requested */
    if (print_stats_flag && !output_json) {
        print_stats(&stats, fout);
    }
    
    /* Flush output before cleanup */
    fflush(fout);
    
    /* Clean up */
    if (close_fin) fclose(fin);
    if (close_fout) fclose(fout);
    
    /* Free sample names from header */
    if (sample_names) {
        for (int i = 0; i < num_samples; i++) {
            free(sample_names[i]);
        }
        free(sample_names);
    }
    
    return 0;
}
