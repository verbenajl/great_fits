#!/bin/bash
# Example: How to use combine_fits command

# Step 1: Create a list file with your FITS files
cat > my_fits_list.txt << 'EOF'
# My FITS files to combine
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part1.fits
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part2.fits
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part3.fits
EOF

# Step 2: Run the combine_fits command
combine_fits --input my_fits_list.txt --output /home/diskB/data_soft/tests/SKYCHOPDIFF/combined_output.fits

# Step 3: Verify the output
print_oifits_info --fits /home/diskB/data_soft/tests/SKYCHOPDIFF/combined_output.fits
