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


def get_pon_bam_paths():
    result = dict()
    for refver, subdict in PON_BAM_PATHS_WITHOUT_NAMES.items():
        result[refver] = dict()
        for cohort, bam_path_list in subdict.items():
            result[refver][cohort] = dict()
            for bam_path in bam_path_list:
                sampleid = 'PON_' + os.path.basename(bam_path).split('.')[0]
                result[refver][cohort][sampleid] = bam_path

    return common.RefverDict(result)


PON_BAM_PATHS = get_pon_bam_paths()


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



#def get_ponfilter(refver, cohortname, **kwargs):
#    return PonFilter(
#        samples=sorted(PON_BAM_PATHS[refver][cohortname].keys()), 
#        **kwargs,
#    )
#
#
#class PonFilter(SamplewiseFilter):
#    @common.get_deco_arg_choices({"mode": ("mean", "max", "median")})
#    def __init__(
#        self,
#        samples,
#
#        valid_sample_num_cutoff=20,
#        min_germline_vaf=0.2,
#        germline_sample_ratio_cutoff=0.1,
#        # max_noise_vaf=0.1,
#        lowest_subset_fraction=0.7,
#        lowest_subset_num_cutoff=10,
#        lowest_subset_snr_cutoff=5,
#        nearby_ratio=1.3,
#        #nearby_subset_fraction_cutoff=0.3,
#        nearby_subset_num_cutoff=5,
#        nearby_subset_snr_cutoff=3,
#        mode="mean",
#        verbose=False,
#    ):
#        # set logger
#        self.logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)
#
#        # sanity check
#        if nearby_ratio < 1:
#            raise Exception(
#                f'"nearby_ratio" argument must be equal to or greater than 1.'
#            )
#
#        # main
#        self.samples = samples
#        self.params = {
#            "valid_sample_num_cutoff": valid_sample_num_cutoff,
#            "min_germline_vaf": min_germline_vaf,
#            "germline_sample_ratio_cutoff": germline_sample_ratio_cutoff,
#            # "max_noise_vaf": max_noise_vaf,
#            "lowest_subset_fraction": lowest_subset_fraction,
#            "lowest_subset_num_cutoff": lowest_subset_num_cutoff,
#            "lowest_subset_snr_cutoff": lowest_subset_snr_cutoff,
#            "nearby_ratio": nearby_ratio,
#            #"nearby_subset_fraction_cutoff": nearby_subset_fraction_cutoff,
#            "nearby_subset_num_cutoff": nearby_subset_num_cutoff,
#            "nearby_subset_snr_cutoff": nearby_subset_snr_cutoff,
#            "mode": mode,
#        }
#
#    def __repr__(self):
#        infostr = ", ".join(
#            [
#                f"samples: {self.samples}",
#                f"params: {self.params}",
#            ]
#        )
#        return f"<PonFilter ({infostr})>"
#
#    @staticmethod
#    def get_vaf_array(allele_index, readstats_dict):
#        """Returns:
#        np.ndarray composed of VAFs of PON samples, where samples in
#            "excluded_pon_samples" and samples with nan vafs are excluded.
#        Result is sorted in ascending order.
#        """
#        pon_vafs = list()
#        for sampleid, readstats in readstats_dict.items():
#            readstats = readstats_dict[sampleid]
#            denom = readstats.get_total_rppcount()
#            if denom == 0:
#                vaf = np.nan
#            else:
#                numer = readstats["rppcounts"][allele_index]
#                vaf = numer / denom
#
#            if not np.isnan(vaf):
#                pon_vafs.append(vaf)
#
#        pon_vafs.sort()
#        pon_vaf_array = np.array(pon_vafs)
#
#        return pon_vaf_array
#
#    def get_rppcount_array(self, allele_index, readstats_dict):
#        return np.array([
#            readstats["rppcounts"][allele_index]
#            for readstats in readstats_dict.values()
#        ])
#
#    def get_lowest_summary_count(self, pon_rppcount_array):
#        rppcounts_nonzero = pon_rppcount_array[pon_rppcount_array > 0]
#        subset_num = int(len(rppcounts_nonzero) * self.params["subset_fraction"])
#
#        if subset_num >= self.params["subset_num_cutoff"]:
#            subset = sorted(rppcounts_nonzero)[:subset_num]
#            noise_summary = getattr(np, self.params["mode"])(subset)
#        else:
#            noise_summary = None
#
#        return noise_summary
#
#    def check_greater_than_lowest_count(self, query_rppcount, pon_rppcount_array):
#        noise_summary = self.get_lowest_summary_count(pon_rppcount_array)
#        if noise_summary is None:
#            return True
#        else:
#            snr = query_rppcount / noise_summary
#            return snr >= self.params["lowest_subset_snr_cutoff"]
#
#    def get_lowest_summary_vaf(self, pon_vaf_array):
#        subset_num = int(pon_vaf_array.shape[0] * self.params["lowest_subset_fraction"])
#        if subset_num == 0:
#            lowest_summary = None
#        else:
#            subset = pon_vaf_array[:subset_num]
#            lowest_summary = getattr(np, self.params["mode"])(subset)
#
#        return lowest_summary
#
#    def check_greater_than_lowest_vaf(self, query_vaf, pon_vaf_array):
#        if np.isnan(query_vaf):
#            return True
#        else:
#            lowest_summary = self.get_lowest_summary_vaf(pon_vaf_array)
#            self.logger.debug(f"lowest_summary {lowest_summary}")
#
#            if (lowest_summary is None) or (lowest_summary == 0):
#                return True
#            else:
#                snr = query_vaf / lowest_summary
#                self.logger.debug(f"snr {snr}")
#
#                return snr >= self.params["lowest_subset_snr_cutoff"]
#
#    def check_greater_than_nearby(self, query_vaf, pon_vaf_array):
#        if np.isnan(query_vaf):
#            return True
#        else:
#            upper_bound = query_vaf * self.params["nearby_ratio"]
#            lower_bound = query_vaf / self.params["nearby_ratio"]
#            vafs_nearby_query = pon_vaf_array[
#                np.logical_and(pon_vaf_array <= upper_bound, pon_vaf_array >= lower_bound)
#            ]
#            self.logger.debug(f"upper_bound, {upper_bound}")
#            self.logger.debug(f"lower_bound, {lower_bound}")
#            self.logger.debug(f"vafs_nearby_query {vafs_nearby_query}")
#
#            if len(vafs_nearby_query) == 0:
#                return True
#            else:
#                #if len(vafs_nearby_query) >= len(pon_vafs) * self.params["nearby_subset_fraction_cutoff"]:
#                if len(vafs_nearby_query) >= self.params["nearby_subset_num_cutoff"]:
#                    nearby_summary = getattr(np, self.params["mode"])(vafs_nearby_query)
#                    self.logger.debug(f"nearby_summary {nearby_summary}")
#
#                    if nearby_summary == 0:
#                        return True
#                    else:
#                        snr = query_vaf / nearby_summary
#                        self.logger.debug(f'snr {snr}')
#                        return snr >= self.params["nearby_subset_snr_cutoff"]
#                else:
#                    return True
#
#    # @functools.cache
#    def check_germline_pattern(self, pon_vaf_array):
#        n_germline_vaf = (pon_vaf_array >= self.params["min_germline_vaf"]).sum()
#        germline_sample_ratio = n_germline_vaf / pon_vaf_array.shape[0]
#        self.logger.debug(f"germline_sample_ratio {germline_sample_ratio}")
#
#        return germline_sample_ratio >= self.params["germline_sample_ratio_cutoff"]
#
#    def check(
#        self,
#        vp,
#        sampleid,
#        allele_index=1,
#        excluded_pon_samples=None,
#        exclude_target=True,
#        do_germlinepat=True,
#        do_lowest=True,
#        do_nearby=True,
#    ):
#        # set excluded_pon_samples
#        if excluded_pon_samples is None:
#            excluded_pon_samples = set()
#        else:
#            excluded_pon_samples = set(excluded_pon_samples)
#
#        if exclude_target:
#            excluded_pon_samples.add(sampleid)
#
#        used_pon_samples = set(self.samples).difference(excluded_pon_samples)
#        if len(used_pon_samples) < self.params["valid_sample_num_cutoff"]:
#            self.logger.warning(
#                f'The number of valid PON samples ({len(used_pon_samples)}) is '
#                f'less than the parameter "valid_sample_num_cutoff" ({self.params["valid_sample_num_cutoff"]}).'
#            )
#
#        # sanity check
#        if not used_pon_samples.issubset(set(vp.readstats_dict.keys())):
#            raise Exception(f'Not all of PON sample IDs are included in the VariantPlus object.')
#    
#        # set params
#        readstats_dict = {
#            sampleid: vp.readstats_dict[sampleid] for sampleid in used_pon_samples
#        }
#        #pon_rppcount_array = self.get_rppcount_array(allele_index, readstats_dict)
#        pon_vaf_array = self.get_vaf_array(allele_index, readstats_dict)  # nan is removed from "pon_vaf_array"
#        query_vaf = vp.get_vaf(sampleid, allele_index)
#        #query_rppcount = vp.readstats_dict[sampleid]["rppcounts"][allele_index]
#
#        # main
#        subtests = list()
#        if do_germlinepat:
#            subtests.append(
#                not self.check_germline_pattern(pon_vaf_array)
#            )
#        if do_lowest:
#            subtests.append(
#                # is_gt_noise = self.check_greater_than_lowest_count(query_rppcount, pon_rppcount_array)
#                self.check_greater_than_lowest_vaf(query_vaf, pon_vaf_array)
#            )
#        if do_nearby:
#            subtests.append(self.check_greater_than_nearby(query_vaf, pon_vaf_array))
#
#        return all(subtests)


# unused function versions

#@common.get_deco_arg_choices({"mode": ("mean", "max", "median")})
#def filter_pon(
#    vp, sampleid, allele_index, pon_samples,
#    exclude_target=True,
#    do_germlinepat=True,
#    do_lowest=True,
#    do_nearby=True,
#    verbose=True,
#
#    mode="mean",
#    min_germline_vaf=0.2,
#    germline_sample_ratio_cutoff=0.1,
#    lowest_subset_fraction=0.7,
#    lowest_subset_num_cutoff=10,
#    lowest_subset_snr_cutoff=5,
#    nearby_ratio=1.3,
#    nearby_subset_num_cutoff=5,
#    nearby_subset_snr_cutoff=3,
#):
#    # set logger
#    logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)
#
#    # set params
#    used_pon_samples = set(pon_samples)
#    if exclude_target:
#        used_pon_samples.difference_update({sampleid})
#
#    if len(used_pon_samples) < SUFFICIENT_PON_SAMPLE_NUM:
#        logger.warning(
#            f'The number of valid PON samples ({len(used_pon_samples)}) is '
#            f'less than the parameter "SUFFICIENT_PON_SAMPLE_NUM" '
#            f'({SUFFICIENT_PON_SAMPLE_NUM}).'
#        )
#
#    # sanity check
#    if not used_pon_samples.issubset(set(vp.readstats_dict.keys())):
#        raise Exception(f'Not all of PON sample IDs are included in the VariantPlus object.')
#
#    readstats_dict = {
#        sampleid: vp.readstats_dict[sampleid] for sampleid in used_pon_samples
#    }
#    #pon_rppcount_array = self.get_rppcount_array(allele_index, readstats_dict)
#    pon_vaf_array = get_vaf_array(allele_index, readstats_dict)  # nan is removed from "pon_vaf_array"
#
#    #query_rppcount = vp.readstats_dict[sampleid]["rppcounts"][allele_index]
#    query_vaf = vp.get_vaf(sampleid, allele_index)
#
#    # main
#    subtests = list()
#    if do_germlinepat:
#        subtests.append(
#            not check_germline_pattern(
#                pon_vaf_array, verbose, min_germline_vaf, germline_sample_ratio_cutoff, logger,
#            )
#        )
#    if do_lowest:
#        subtests.append(
#            # is_gt_noise = self.check_greater_than_lowest_count(query_rppcount, pon_rppcount_array)
#            check_greater_than_lowest_vaf(query_vaf, pon_vaf_array, logger, lowest_subset_fraction, mode, lowest_subset_snr_cutoff)
#        )
#    if do_nearby:
#        subtests.append(
#            check_greater_than_nearby(query_vaf, pon_vaf_array, logger, nearby_ratio, nearby_subset_num_cutoff, mode, nearby_subset_snr_cutoff)
#        )
#
#    return all(subtests)
#
#
###################
## helper functions
#
#
#def get_vaf_array(allele_index, readstats_dict):
#    """Returns:
#    A np.ndarray composed of VAFs of PON samples. nan is removed.
#    Result is sorted in ascending order.
#    """
#    pon_vafs = list()
#    for readstats in readstats_dict.values():
#        denom = readstats.get_total_rppcount()
#        if denom == 0:
#            continue
#
#        numer = readstats["rppcounts"][allele_index]
#        vaf = numer / denom
#        pon_vafs.append(vaf)
#
#    pon_vafs.sort()
#    pon_vaf_array = np.array(pon_vafs)
#
#    return pon_vaf_array
#
#
#def get_rppcount_array(allele_index, readstats_dict):
#    return np.array([
#        readstats["rppcounts"][allele_index]
#        for readstats in readstats_dict.values()
#    ])
#
#
#######
## check_germline_pattern
#
#def check_germline_pattern(pon_vaf_array, verbose, min_germline_vaf, germline_sample_ratio_cutoff, logger):
#    n_germline_vaf = (pon_vaf_array >= min_germline_vaf).sum()
#    germline_sample_ratio = n_germline_vaf / pon_vaf_array.shape[0]
#    logger.debug(f"germline_sample_ratio {germline_sample_ratio}")
#
#    return germline_sample_ratio >= germline_sample_ratio_cutoff
#
#
#######
## check_greater_than_lowest_vaf
#
#def check_greater_than_lowest_vaf(query_vaf, pon_vaf_array, logger, lowest_subset_fraction, mode, lowest_subset_snr_cutoff):
#    if np.isnan(query_vaf):
#        return True
#    else:
#        lowest_summary = get_lowest_summary_vaf(pon_vaf_array, lowest_subset_fraction, mode)
#        logger.debug(f"lowest_summary {lowest_summary}")
#
#        if (lowest_summary is None) or (lowest_summary == 0):
#            return True
#        else:
#            snr = query_vaf / lowest_summary
#            logger.debug(f"snr {snr}")
#
#            return snr >= lowest_subset_snr_cutoff
#
#
#def get_lowest_summary_vaf(pon_vaf_array, lowest_subset_fraction, mode):
#    subset_num = int(pon_vaf_array.shape[0] * lowest_subset_fraction)
#    if subset_num == 0:
#        lowest_summary = None
#    else:
#        subset = pon_vaf_array[:subset_num]
#        lowest_summary = getattr(np, mode)(subset)
#
#    return lowest_summary
#
#
#######
## check_greater_than_nearby
#
#def check_greater_than_nearby(query_vaf, pon_vaf_array, logger, nearby_ratio, nearby_subset_num_cutoff, mode, nearby_subset_snr_cutoff):
#    if np.isnan(query_vaf):
#        return True
#    else:
#        upper_bound = query_vaf * nearby_ratio
#        lower_bound = query_vaf / nearby_ratio
#        vafs_nearby_query = pon_vaf_array[
#            np.logical_and(
#                pon_vaf_array <= upper_bound, 
#                pon_vaf_array >= lower_bound
#            )
#        ]
#        logger.debug(f'upper_bound {upper_bound}')
#        logger.debug(f'lower_bound {lower_bound}')
#        logger.debug(f'vafs_nearby_query {vafs_nearby_query}')
#
#        if len(vafs_nearby_query) == 0:
#            return True
#        else:
#            if len(vafs_nearby_query) >= nearby_subset_num_cutoff:
#                nearby_summary = getattr(np, mode)(vafs_nearby_query)
#                logger.debug(f'nearby_summary {nearby_summary}')
#
#                if nearby_summary == 0:
#                    return True
#                else:
#                    snr = query_vaf / nearby_summary
#                    logger.debug(f'snr {snr}')
#                    return snr >= nearby_subset_snr_cutoff
#            else:
#                return True
#
#
#######
## other unused but kept ones
#
#def get_lowest_summary_count(pon_rppcount_array, subset_fraction, subset_num_cutoff, mode):
#    rppcounts_nonzero = pon_rppcount_array[pon_rppcount_array > 0]
#    subset_num = int(len(rppcounts_nonzero) * subset_fraction)
#
#    if subset_num >= subset_num_cutoff:
#        subset = sorted(rppcounts_nonzero)[:subset_num]
#        noise_summary = getattr(np, mode)(subset)
#    else:
#        noise_summary = None
#
#    return noise_summary
#
#
#def check_greater_than_lowest_count(query_rppcount, pon_rppcount_array, subset_fraction, subset_num_cutoff, mode, lowest_subset_snr_cutoff):
#    noise_summary = get_lowest_summary_count(pon_rppcount_array, subset_fraction, subset_num_cutoff, mode)
#    if noise_summary is None:
#        return True
#    else:
#        snr = query_rppcount / noise_summary
#        return snr >= lowest_subset_snr_cutoff

###########################


