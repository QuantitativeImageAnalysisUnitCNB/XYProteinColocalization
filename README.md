# XYProteinColocalization
This is a Java-based script to measure colocalization between X protein and Y protein within specific cell regions or whole cell.
## Overview of Procedure

- The basic workflow initially consists of typing both the `path to the input directory` (in which the raw images to be analyzed are located) and the `path to the output directory` (in which the results generated from colocalization analysis will be stored). 
- This Java-based script which works within ImageJ or Fiji ecosystem performs pixel-by-pixel intensity correlation methods in batch-mode in two different proteins by automatically selecting a specific region of interest (ROI) within cell or overthe whole cell, depending on what was predefined by user.
- Once done, this tool is able to automatically extract each image from the input directory file as a multichannel z-stack TIFF format to be analyzed. 
- Then the command `“Splitting multi-channel images”` is called to split each z-stack volume into their respective channel components then being each resulting stack manipulated for extracting each slice to be analyzed independently. 
- To isolate either ROIs or whole cells as specific areas for colocalization analysis, `“Auto-Threshold”` plugin is called in order to binarize images using the `“Default”` method subsequently outlining those areas, then `“Create Selection”` command is called in order to keep them as ROIs thus computing from them their corresponding features (morphological, geometrical and descriptive statistics). 
- Whether DAPI signal is available on image, ROIs whose contour’s centroid are contained in areas belonging to DAPI signal, they will be considered as relevant, conversely, whether there is no DAPI signal on images, all ROIs will be considered as relevant to be analyzed.
- Then the `“Combine”` command from `“ROI Manager”` is called combining all ROIs by applying the union operator.
- Once done, to perform colocalization analysis between X protein and Y protein and in order to keep only those pixels within ROI having pixel-intensity values among [media-std] and [median+std], filtering actions will be applied over the two channels. At this point, mean intensity values for each protein together with colocalization indicators are computed.


The first one is the Pearson’s coefficient is a statistic which provides an estimate of how strong is the association of pixel-by-pixel intensity signals among X protein and Y protein. This pixel-by-pixel intensity relationship is described by the slope of linear regression. Then the slope measures the overlapping degree of the signal levels of two proteins ranging from 1 (fluorescence intensities from X and Y proteins are perfectly and linearly related) to -1 (fluorescence intensities from X and Y proteins are perfectly but inversely related) and 0 standing for no correlation. The main advantage of this coefficient as a evaluator for colocalization is its simplicity to subtracts the mean intensity pixel-by-pixel being thus independent of signal levels and signal from background.
## Pearson's Coefficient

$$r_\mathrm{p} =  \frac {\sum ((X_i - X_{\mathrm{avg}}) (Y_i - Y_{\mathrm{avg}}))} {\sqrt{\sum (X_i-X_{\mathrm{avg}})^2 \sum (Y_i - Y_{\mathrm{avg}})^2}}$$



with $X_i$ and  $Y_i$ being the pixel values with the pixel index and Xavg and Yavg the averages of the Xi and  Yi protein signals respectively and the summations with indexi over all the image pixels.
## Manders’ overlap coefficients
Secondly, Manders’ overlap coefficients are based on Pearson’s correlation coefficient, overcoming the difficulty of understanding negative Pearson’s correlation coefficients by eliminating the subtraction of mean signals from the equation. These coefficient’s values range from 0 (non-overlapping) to 1 (total overlapping) measuring the fraction of total X protein fluorescence that colocalizes with the fluorescence of Y protein. In this context, M1 is the ratio of the “summed intensities of pixels from the X protein for which the intensity in the Y protein is above zero” to the total intensity in the X protein being M2 defined conversely for Y protein.
$$M_1 = \frac {\sum X_{\mathrm{coloc},i}}  {\sum X_i}$$
where $X_i,coloc$ = $X_i$ if $Y_i$ >0 and $X_i,coloc$ = 0 if $Y_i$ =0
$$M_2 = \frac {\sum Y_{\mathrm{coloc},i}}  {\sum Y_i}$$
where $Y_i,coloc$ = $Y_i$ if $X_i$ >0 and $Y_i,coloc$ = 0 if $X_i$ =0
## Two-sample Kolmogorov Smirnov test
Finally, the p-value from the Two-sample Kolmogorov Smirnov test (statistical significance indicated by *p < 0.05.) is computed to determine if the X and Y protein fluorescence signals follow the same spatial distribution or not.  Therefore this metric provides information which goes well beyond the colocalization analysis as X and Y proteins may not colocalize but still be spatially-dependent on each other. This is done by comparing the distribution of values from the first product among the pixel-by-pixel intensities in which the mean intensity is subtracted from each pixel’s intensity value from each protein and the second product among the normalized pixel-by-pixel intensities from X protein (as detailed above) and its flipped version from Y protein.


