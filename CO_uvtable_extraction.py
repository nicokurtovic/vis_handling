#!/usr/bin/env python
# -*- coding: utf-8 -*-
# /home/kurtovic/Documents/CASA/casa-pipeline-release-5.6.2-6.el7/bin/casa

################################################################################
#                             GAS UV TABLE EXTRACTION                          #
################################################################################

'''
Max Planck Institute for Extraterrestrial Physics
L-G Group

Nicolas Kurtovic
Contact: kurtovic at mpe.mpg.de

The following code extracts the visibilities from a measurement set, and
stores them in txt, to be later read in python or other code. 

If you use this file, please cite Kurtovic & Pinilla 2024.
https://ui.adsabs.harvard.edu/abs/2024A%26A...687A.188K/abstract
'''

# Import the analysis scripts, available in:
# https://casaguides.nrao.edu/index.php/Analysis_Utilities
sys.path.append('/path_to_analysis_scripts/analysis_scripts/')
import analysisUtils as au

# Execute the supporting functions
execfile('./CO_to_ascii.py')


################################################################################
#                                                                              #
################################################################################

# Prefix for naming files
prefix = 'PDS111'

# Paths to measurement sets
work_dir = os.getcwd() + '/'
msfiles_dir  = work_dir + 'msfiles/'
uvtables_dir = work_dir + 'uvtables/'

# Measurement set
gas_msfile = msfiles_dir + 'PDS111_SB_12CO_selfcal.ms.contsub.cvel'

# Line information
line = '12CO'
freq_line = 230.5380000e9 # 12CO J=2-1, Hz
freq_GHz  = str(freq_line*1e-9)+'GHz'

# Initialize weight column
initweights(vis=gas_msfile, wtmode='weight', dowtsp=True)

# Extract flagged visibilities? False for no, True for yes
keepflags = False


################################################################################
#                              UVTABLE EXTRACTION                              #
################################################################################

# From here onwards, you should not need to modify anything

# Array to save frequencies
freqs = []
# Extract each channel, one by one, and save the uvtable
for i in range(au.getNChanFromCaltable(gas_msfile)[0]):
    # Print channel number
    print(i)
    # Name of ascii file to save
    ascii_file = uvtables_dir+prefix+'_'+line+'_chan'+str(i)+'.txt'
    # Temporary measurement set for extracting visibilities
    ms_file    = msfiles_dir + prefix + '_channel.ms'

    # Split the channel
    os.system('rm -rf ' + ms_file)
    split(vis=gas_msfile,      \
          outputvis=ms_file,   \
          spw='0:'+str(i),     \
          keepflags=keepflags, \
          datacolumn='DATA')

    # Get frequency of the channel, and save
    freq_chan = get_freqchan(ms_file)
    freqs.append(freq_chan)

    # Remove uvtable if it already exists, to avoid overwriting
    os.system('rm -rf '+ascii_file)
    # Write the uvtable in txt. 
    ms_to_ascii(ms_file, ascii_file, with_flags=True)
    # Delete temporary msfile
    os.system('rm -rf ' + ms_file)

# Save frequency array
np.save(uvtables_dir+prefix+'_frequency_'+line, np.squeeze(freqs))
# Read frequency array
freqs = np.load(uvtables_dir+prefix+'_frequency_'+line+'.npy')

# Light speed
light_speed = 299792458.0 # m/s
# Calculate velocity relative to line rest frequency
vels = -(((freqs / freq_line) * light_speed) - light_speed).astype(int)
# Save velocity array
np.save(uvtables_dir+prefix+'_velocity_'+line, np.squeeze(vels))
# Read velocity array
vels = np.load(uvtables_dir+prefix+'_velocity_'+line+'.npy')

# Check weight of the visibilities
os.system('du -sch '+uvtables_dir)
print (' ')

# Save for your own reference
# 173M	./uvtables/
# 173M	total

