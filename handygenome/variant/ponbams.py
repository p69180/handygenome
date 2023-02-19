import os
import re

import pysam

import handygenome.common as common
import handygenome.bameditor as bameditor


PON_BAM_PATHS_WITHOUT_NAMES = common.RefverDict({
    'GRCh37': {
        'BGI': [
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG01.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG02.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG03.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG04.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG05.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG06.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG07.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG08.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG09.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG10.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG11.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG12.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG13.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG14.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG15.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG16.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG17.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG18.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG19.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG20.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG21.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG22.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG23.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/BGI_BAM/BGI-WG24.normal.rmBDBI.cram',
        ],
        'PCAWG': [
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4389.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4395.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4396.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4397.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4398.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-4420.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-05-5429.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-38-4628.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-44-2659.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-44-6148.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-49-4486.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-49-4512.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-49-6742.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-50-5930.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-50-5932.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-50-6591.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-50-6597.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-6972.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-6982.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-6984.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-6986.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-7281.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-55-8299.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-64-1678.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-64-1680.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-67-3771.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-67-6215.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-73-4659.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-73-4666.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-75-5147.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-75-6203.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-75-7030.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-78-7158.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-78-7535.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-91-6840.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-91-6847.normal.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA_BAM/TCGA-97-8171.normal.cram',
        ],
        'SNULUNG': [
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-14.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-4.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-42.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-44.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-51.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-57.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-6.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-69.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-78.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-87.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-89.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F13.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F16.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F2.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F29.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F31.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F37.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-F6.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF1.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF104.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF107.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF115.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF13.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF20.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF21.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF23.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF24.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF3.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF34.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF37.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF4.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF43.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF53.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF55.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF56.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF57.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF58.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF62.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF67.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF71.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF76.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF77.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF78.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF79.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF85.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-FF86.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-SC126.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-SC81.normal.rmBDBI.cram',
            '/home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/SNUH_BAM/LU-SC97.normal.rmBDBI.cram',
        ],
        'NEBlow_FFPE_LCM': [
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-1-A.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-1-B.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-11.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-5-A.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-6-B.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/1015-6.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FA-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FB-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FC-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FD-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FE-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FF-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FH-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/FI-N-1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/Thy1-B.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/Thy2-A.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/08_Thyroid_LCM/01_aligned_bam/Thy2-B.al.st.dedup.realign.bam',
        ],
        'NEBlow_FreshFrozen_LCM': [
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_C1.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_C2.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_C3.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_C4.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_C5.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_L1.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_L2.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_L3.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_L4.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_L5.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R1.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R10_lane1.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R11.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R12.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R2.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R3.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R4.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R5.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R6.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R7.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R8.fmarked.realigned.bam',
            '/home/users/team_projects/Lineage_tracing/DB3/02-2_BAM_liverLCM/DB3_liver_R9.fmarked.realigned.bam',
        ],
        'NEBlow_PAXgene_organoid': [
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_AIR_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_AIR_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_AIR_BAS_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_P1_P.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_P2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_P2_P.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_AT_P3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_C_B_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_P_B_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_P_B_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_P_B_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES1_P_B_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_AT_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_1_BAS_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_AT_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_2_BAS_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_3_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_3_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_3_BAS_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/ES6_3_BAS_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/NS2_1_AT_P1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/NS2_1_AT_P1_P.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/NS2_1_AT_P2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/NS2_1_AT_P2_P.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/RSPO1_merged_P1_P4.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/S0_1_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/S0_1_AT_3.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/S0_2_AT_1.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/S0_2_AT_2.al.st.dedup.realign.bam',
            '/home/users/twkim/00_data/10_smoking_in_lung/01_aligned_bam/DNA/S0_2_AT_3.al.st.dedup.realign.bam',
        ],
    },
})


def _get_pon_bam_paths():
    result = dict()
    for refver, subdict in PON_BAM_PATHS_WITHOUT_NAMES.items():
        result[refver] = dict()
        for cohort, bam_path_list in subdict.items():
            result[refver][cohort] = dict()
            for bam_path in bam_path_list:
                sampleid = 'PON_' + os.path.basename(bam_path).split('.')[0]
                result[refver][cohort][sampleid] = bam_path

    return common.RefverDict(result)


PON_BAM_PATHS = _get_pon_bam_paths()


def get_pon_sample_names(cohort_names, refver):
    """Args:
        cohort_names: list or tuple of cohort names
    """
    result = list()
    pon_bam_paths = PON_BAM_PATHS[refver]
    for cohort in cohort_names:
        if cohort not in pon_bam_paths.keys():
            raise Exception(
                f'Available PON cohort name for reference version {refver} is: '
                f'{tuple(pon_bam_paths.keys())}'
            )
        result.extend(sorted(pon_bam_paths[cohort].keys()))

    return result



