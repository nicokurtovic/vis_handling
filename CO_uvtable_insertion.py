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

The following code insert model or residual visibilities written in txt
into a measurement set. This code does not erase or replace the original
measurement set.

If you use this file, please cite Kurtovic & Pinilla 2024.
https://ui.adsabs.harvard.edu/abs/2024A%26A...687A.188K/abstract
'''

# Import the analysis scripts, available in:
# https://casaguides.nrao.edu/index.php/Analysis_Utilities
sys.path.append('./analysis_scripts/')
import analysisUtils as au

# Execute the supporting functions
execfile('./CO_to_ascii.py')


################################################################################
#                                                                              #
################################################################################

# Prefix for naming files
prefix = 'PDS111'

# Paths
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
# Should be the same as the extracted visibilities
keepflags = False

# Prefix for the model visibilities. They must be in the uvtables_dir
# In this case, we will assume the model visibilities are called
# PDS111_12CO_mod{chan_number}.txt, so that channel
# 0 is called PDS111_12CO_mod0.txt, channel 1 is PDS111_12CO_mod1.txt, 
# channel 14 is PDS111_12CO_mod14.txt, and so on.
prefix_mod = prefix + '_' + line + '_mod'

# Name of the model measurement set to be created. 
# It will be created in the msfiles_dir
# In this example, the name will be PDS111_12CO_mod.ms
new_mod_ms = prefix_mod + '.ms'


################################################################################
#                               UVTABLE INSERTION                              #
################################################################################

# From here onwards, you should not need to modify anything

# Names
for i in range(au.getNChanFromCaltable(gas_msfile)[0]):
    # Print channel
    print(i)
    # Name of ascii file and ms files
    ascii_file  = uvtables_dir + prefix_mod+str(i) + '.txt'
    ms_file     = msfiles_dir  + prefix + '_channel.ms'
    ms_file_new = msfiles_dir  + prefix + '_auxmod.ms'

    # Split the channel
    os.system('rm -rf ' + ms_file)
    os.system('rm -rf ' + ms_file_new)
    split(vis=gas_msfile, \
          outputvis=ms_file, \
          spw='0:'+str(i), \
          keepflags=keepflags, datacolumn='DATA')

    # Load and overwrite visibilities
    ascii_to_ms(ascii_file, ms_file, ms_file_new)

    # If it is the first channel, then rename ms file.
    if i == 0:
        # Change name
        os.system('rm -rf ' + msfiles_dir + new_mod_ms)
        os.system('mv ' + ms_file_new + ' ' + msfiles_dir + new_mod_ms)
    # If not the first channel, then concatenate
    else:
        concat(vis=[ms_file_new, \
                    msfiles_dir + new_mod_ms], \
               concatvis=msfiles_dir + new_mod_ms, \
               freqtol='1Hz')
    # Remove temporal files
    os.system('rm -rf ' + ms_file)
    os.system('rm -rf ' + ms_file_new)


