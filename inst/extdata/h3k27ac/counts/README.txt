Per Doug's request, files were formatted in this way:

  - The first 6 columns were BED formatted
    - (i.e. chrom, chromStart, chromEnd, name, score, strand)
    - The rest of the columns are the count columns

  - H3K27Ac columns were renamed as follows:

	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0000_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0000_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0030_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0030_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0060_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0060_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0090_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0090_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0120_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0120_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0240_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0240_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0360_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_0360_2.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_1440_1.1.1_filter_sorted.bam"
	"LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_1440_2.1.1_filter_sorted.bam"

     to

	"k27_0000_r1"  
        "k27_0000_r2"  
        "k27_0030_r1"   
        "k27_0030_r2"  
        "k27_0060_r1"  
        "k27_0060_r2"  
        "k27_0090_r1"  
        "k27_0090_r2"  
        "k27_0120_r1"  
        "k27_0120_r2"  
        "k27_0240_r1"  
        "k27_0240_r2"  
        "k27_0360_r1"  
        "k27_0360_r2"  
        "k27_1440_r1"  
        "k27_1440_r2"  