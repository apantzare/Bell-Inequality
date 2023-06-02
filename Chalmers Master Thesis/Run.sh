#!/usr/bin/env bash
#SBATCH -A C3SE2023-1-7
#SBATCH -t 2-10:00:00
#SBATCH -n 128
#SBATCH -C ICELAKE
#SBATCH --job-name=Axels_test

module purge
module load PYTHIA/8.309-foss-2022b
module load ROOT/6.26.10-foss-2022b

# Copy data to TMPDIR
cd $TMPDIR
cp /cephyr/users/paaxel/Vera/main91 .

# Run your code
./main91
hadd New_Data.root New_Data_0.root New_Data_1.root New_Data_2.root New_Data_3.root New_Data_4.root New_Data_5.root New_Data_6.root New_Data_7.root New_Data_8.root New_Data_9.root New_Data_10.root New_Data_11.root New_Data_12.root New_Data_13.root New_Data_14.root New_Data_15.root New_Data_16.root New_Data_17.root New_Data_18.root New_Data_19.root New_Data_20.root New_Data_21.root New_Data_22.root New_Data_23.root New_Data_24.root New_Data_25.root New_Data_26.root New_Data_27.root New_Data_28.root New_Data_29.root New_Data_30.root New_Data_31.root New_Data_32.root New_Data_33.root New_Data_34.root New_Data_35.root New_Data_36.root New_Data_7.root New_Data_38.root New_Data_39.root New_Data_40.root New_Data_41.root New_Data_42.root New_Data_43.root New_Data_44.root New_Data_45.root New_Data_46.root New_Data_47.root New_Data_48.root New_Data_49.root New_Data_50.root New_Data_51.root New_Data_52.root New_Data_53.root New_Data_54.root New_Data_55.root New_Data_56.root New_Data_57.root New_Data_58.root New_Data_59.root New_Data_60.root New_Data_61.root New_Data_62.root New_Data_63.root New_Data_64.root New_Data_65.root New_Data_66.root New_Data_67.root New_Data_68.root New_Data_69.root New_Data_70.root New_Data_71.root New_Data_72.root New_Data_73.root New_Data_74.root New_Data_75.root New_Data_76.root New_Data_77.root New_Data_78.root New_Data_79.root New_Data_80.root New_Data_81.root New_Data_82.root New_Data_83.root New_Data_84.root New_Data_85.root New_Data_86.root New_Data_87.root New_Data_88.root New_Data_89.root New_Data_90.root New_Data_91.root New_Data_92.root New_Data_93.root New_Data_94.root New_Data_95.root New_Data_96.root New_Data_97.root New_Data_98.root New_Data_99.root New_Data_100.root New_Data_101.root New_Data_102.root New_Data_103.root New_Data_104.root New_Data_105.root New_Data_106.root New_Data_107.root New_Data_108.root New_Data_109.root New_Data_110.root New_Data_111.root New_Data_112.root New_Data_113.root New_Data_114.root New_Data_115.root New_Data_116.root New_Data_117.root New_Data_118.root New_Data_119.root New_Data_120.root New_Data_121.root New_Data_122.root New_Data_123.root New_Data_124.root New_Data_125.root New_Data_126.root New_Data_127.root
# Retrieve results
cp New_Data.root $SLURM_SUBMIT_DIR
