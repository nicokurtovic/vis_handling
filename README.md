---
title: Vishandling
intro: Collection of codes to extract, handle, and insert visibilities from gas measurement sets to txt files, in CASA. 
---

## Extract the visibilities of your gas emission

The gas emission in a measurement set (from now on, ms files) is contained in channels, which are grouped in spectral windows. For this example, we will assume that only *one spectral window* is present in the ms file. For multiple spectral windows, please check what to do if I have multiple spectral windows.

The necessary functions to extract the visibilities of each channel are in the file *CO_to_ascii.py*, which you should not need to modify. The first step in extracting the visibilities is to execute this file, as well as importing the analysis utilities of CASA (CASA Team et al., 2022; Hunter et al., 2023). 

<code>\`\`\`
testing code
<code>

```
testing code
```
