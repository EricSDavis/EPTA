Per Doug's request, files were formatted in this way:

  - The first 6 columns were BED formatted
    - (i.e. chrom, chromStart, chromEnd, name, score, strand)
    - The rest of the columns are the count columns

  - ATAC columns were renamed as follows:

    	"LIMA_ATAC_THP1_WT_LPIF_S_0000_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0000_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0030_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0030_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0060_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0060_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0090_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0090_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0120_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0120_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0240_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0240_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0360_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_0360_3.1.2"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_1440_2.1.1"  
     	"LIMA_ATAC_THP1_WT_LPIF_S_1440_3.1.2"  

     to

     	"atac_0000_r1"  
     	"atac_0000_r2"  
     	"atac_0030_r1"  
     	"atac_0030_r2"  
     	"atac_0060_r1"  
     	"atac_0060_r2"  
     	"atac_0090_r1"  
     	"atac_0090_r2"  
     	"atac_0120_r1"  
     	"atac_0120_r2"  
     	"atac_0240_r1"  
     	"atac_0240_r2"  
     	"atac_0360_r1"  
     	"atac_0360_r2"  
     	"atac_1440_r1"  
     	"atac_1440_r2"  