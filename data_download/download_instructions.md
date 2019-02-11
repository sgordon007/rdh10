directory: /190111_100PE_HS4K2A

:: RUN INFORMATION ::
100PE HiSeq4000 (minimum 2% PhiX control run on all libraries unless otherwise stated, up to 25% for diversity)

:: LANE ASSIGNMENTS ::
Lane 1: Alison Hanh Nguyen and Bachtrog Lab
Lane 2: Joshua Tworig and Feller Lab
Lane 3: Marta Vuckovic and Napoli Lab
Lane 4: Marta Vuckovic and Napoli Lab
Lane 5: Jen Quick-Cleveland and Ares Lab
Lane 6: Steve Eacker
Lane 7: Steve Eacker
Lane 8: Beeke Wienert and Corn Lab

:: LANE SUMMARY ::
Lane 1:  DBAN005 + 2% PhiX, 3nM, INDEX (369M Reads, 1.2% PhiX Aligned)
Lane 2:  MFJT01 + 2% PhiX, 3nM, INDEX (325M Reads, 2.3% PhiX Aligned)
Lane 3:  JNVM01,02,02 + 2% PhiX, 3nM, INDEX (394M Reads, 1.7% PhiX Aligned)
Lane 4:  JNVM01,02,02 + 2% PhiX, 3nM, INDEX (394M Reads, 1.6% PhiX Aligned)
Lane 5:  MJO_001 + 2% PhiX, 3nM, INDEX (361M Reads, 1.8% PhiX Aligned)
Lane 6:  ILSE37 + 0% PhiX, 3nM, INDEX (372M Reads, 0% PhiX Aligned)
Lane 7:  ILSE39 + 0% PhiX, 3nM, INDEX (378M Reads, 0% PhiX Aligned)
Lane 8:  JCBW_20190106 + 2% PhiX, 3nM, INDEX (341M Reads, 1.9% PhiX Aligned)

GSL uses the Illumina supported CASAVA bcl2fastq (v2.19) program to demultiplex data. While a set file naming system has not been determined please refer to either the directory associated with your lane assignment noted above or refer to the directory associated with the submitting user or laboratory PI. However, file naming is typically based on the library/sample name and the lane number. Within these directories, compressed FASTQ files will be listed with a single other file (laneBarcode.html) which displays the demultiplexing statistics for all sequencing lanes that were demultiplexed together. (This file is not always generated for MiSeq runs if the demultiplexing occurs onboard the instrument). Scroll through this file until you see your lane and sample(s).


cat ftp_file_list.txt | awk '{print $5}' > ftp_file_sizecol.txt


rsync -av user@server:$BSCRATCH/mus_2019_01_29/ /Volumes/Seagate\ Backup\ Plus\ Drive/mus_2019_01_29/

rsync -av ~/Downloads/kallisto-0.45.0/test user@server:$BSCRATCH/mus_2019_01_29/test_data/

