
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JTable;
import javax.swing.table.TableModel;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

public class XYProteinColocalization_ implements PlugIn {

	Thread thread;
	KolmogorovSmirnovTest testKS;
	PearsonsCorrelation pc;
	ImagePlus imps;

	public XYProteinColocalization_() {

	}

	public void run(String arg0) {
		// Select the directory in which the images to be analyzed are stored as is
		// indicated below
		File imageFolder = new File("/home/acayuela/data/cris_tirso_agnes/imagesToAnalyze");
		File[] listOfFiles = imageFolder.listFiles();

		thread = new Thread(new Runnable() {
			public void run() {
				// Iterating through the files belonging the directory
				for (int x = 0; x < listOfFiles.length; x++) {
					imps = new ImagePlus(imageFolder.getAbsolutePath() + File.separator + listOfFiles[x].getName());
					// Split channels for each image
					ImagePlus[] channels = ChannelSplitter.split(imps);

					ImagePlus chXProtein = channels[0].duplicate();
					ImagePlus chYProtein = channels[1].duplicate();
					ImagePlus chDapi = channels[2].duplicate();
					/// Select this channel in case you want to isolate the analysis to a specific
					/// region of the cell
					ImagePlus chCinet = channels[3].duplicate();
					// Get stack from each channel
					ImagePlus chXProteinSlices[] = stack2images(chXProtein);
					ImagePlus chYProteinSlices[] = stack2images(chYProtein);
					ImagePlus chDapiSlices[] = stack2images(chDapi);
					ImagePlus chCinetSlices[] = stack2images(chCinet);
					// Define the list to store the colocalization indicators
					List<String> pvaluePerSlice = new ArrayList<String>();
					List<String> pCorrelationPerSlice = new ArrayList<String>();
					List<String> m1PerSlice = new ArrayList<String>();
					List<String> m2PerSlice = new ArrayList<String>();
					List<Double> pixelXProteinYProteinMultMirrorTotal = new ArrayList<Double>();
					List<Double> pixelXProteinYProteinMultTotal = new ArrayList<Double>();
					List<Double> pixelXProteinTotal = new ArrayList<Double>();
					List<Double> pixelYProteinTotal = new ArrayList<Double>();
					List<Double> pixelValuesYProteinM1Total = new ArrayList<Double>();
					List<Double> pixelValuesXProteinM2Total = new ArrayList<Double>();
					List<String> meanXProteinTotal = new ArrayList<String>();
					List<String> meanYProteinTotal = new ArrayList<String>();

					// IS cinet included in Dapi?
					// We only keep those specific regions within DAPI to ensure only signal from
					// nuclei
					// To do that we iterate through all slices of stack
					for (int i = 0; i < chXProteinSlices.length; i++) {

						IJ.resetMinAndMax(chXProteinSlices[i]);
						// ImagePlus impToSegment = chXProteinSlices[i].duplicate();
						// We define here the list "pixelXProteinYProteinMult" to store the product
						// among XY pixel normalized intensities
						List<Double> pixelXProteinYProteinMult = new ArrayList<Double>();
						// We define here the list "pixelXProteinYProteinMultMirror" to store the
						// product among XY pixel normalized intensities using the mirror version
						List<Double> pixelXProteinYProteinMultMirror = new ArrayList<Double>();
						IJ.resetMinAndMax(chDapiSlices[i]);
						// Segment DAPI slices using Max Entropy mehtod of Auto Threhold Global methods
						// from MorpholibJ
						IJ.run(chDapiSlices[i], "Auto Threshold", "method=MaxEntropy ignore_black white");
						// Then a selection is create to isolate DAPI nuclei bg
						IJ.run(chDapiSlices[i], "Create Selection", "");
						// Get background from DAPI rois
						Roi roiBg = chDapiSlices[i].getRoi();
						chDapiSlices[i].setRoi(roiBg);
						// Measure the background on DAPI
						double bgMeasure = chDapiSlices[i].getStatistics().mean;
						// The we do the inverse of roiBg to isolate DAPI nuclei
						IJ.run(chDapiSlices[i], "Make Inverse", "");
						// Get DAPI rois not bg
						Roi roiDapi = chDapiSlices[i].getRoi();
						IJ.resetMinAndMax(chCinetSlices[i]);
						// Segmet cinet channel (containing the specific areas for analysis) with
						// Default method from MorpholibJ
						IJ.run(chCinetSlices[i], "Auto Threshold", "method=Default ignore_black white stack");
						// We invert LUT to appy "Fill Holes"
						IJ.run(chCinetSlices[i], "Invert LUT", "");
						IJ.run(chCinetSlices[i], "Fill Holes", "");
						IJ.run(chCinetSlices[i], "Invert LUT", "");
						// We call once again "Invert LUT" to go on with analysis
						// Apply a morphological operation as erosion to isolate specific regions
//						chCinetSlices[i] = new ImagePlus("",
//								Morphology.erosion(chCinetSlices[i].getProcessor(), Strel.Shape.DISK.fromRadius((10))));
						// Apply a median filter to keep only relevant cinet areas
						IJ.run(chCinetSlices[i], "Median...", "radius=4");
						IJ.run(chCinetSlices[i], "Create Selection", "");
						IJ.run(chCinetSlices[i], "Make Inverse", "");
						// Get specific areas for measuring as rois
						Roi roiCinet = chCinetSlices[i].getRoi();
						// Keep them as separated rois
						Roi[] roisCinet = new ShapeRoi(roiCinet).getRois();
						List<Roi> roisCinetDef = new ArrayList<Roi>();

						Roi roiToAdd = null;
						RoiManager rm = RoiManager.getInstance();
						if (null == rm)
							rm = new RoiManager();
						for (int j = 0; j < roisCinet.length; j++) {
							// Check that specific rois are containd into DAPI signal
							for (int z = 0; z < roisCinet[j].getContainedPoints().length; z++) {
								if (roisCinet[j].getContainedPoints().length > 4
										&& new ShapeRoi(roiDapi).contains(roisCinet[j].getContainedPoints()[x].x,
												roisCinet[j].getContainedPoints()[x].y) == true) {
									roiToAdd = roisCinet[j];
								}
							}
							// Add rois which are contained into DAPI rois
							if (roiToAdd != null && roisCinetDef.contains(roiToAdd) == false) {
								roisCinetDef.add(roiToAdd);
								rm.addRoi(roiToAdd);

							}

						}
						// Initializing list to store pixel intensities accomplishing threholds
						List<Double> pixelValuesXProtein0 = new ArrayList<Double>();
						List<Double> pixelValuesYProtein0 = new ArrayList<Double>();
						// Initializing list to store pixel intensities accomplishing threholds
						// (normalized version: subtracting the mean Intensity value)
						List<Double> pixelValuesXProtein1 = new ArrayList<Double>();
						List<Double> pixelValuesYProtein1 = new ArrayList<Double>();
						List<Double> pixelValuesYProteinM1 = new ArrayList<Double>();
						List<Double> pixelValuesXProteinM2 = new ArrayList<Double>();
						List<Double> meanXProtein = new ArrayList<Double>();
						List<Double> meanYProtein = new ArrayList<Double>();
						Roi roiToMeasure = null;

						if (null == rm)
							rm = new RoiManager();
						if (roisCinetDef.size() != 0) {
							for (int j = 0; j < roisCinetDef.size(); j++)
								rm.addRoi(roisCinetDef.get(j));
							int[] indexes = new int[rm.getCount()];
							for (int j = 0; j < rm.getCount(); j++)
								indexes[j] = j;
							rm.setSelectedIndexes(indexes);
							if (indexes.length > 1) {
								// Combine them as an individual roi
								rm.runCommand(chCinetSlices[i], "Combine");
								roiToMeasure = chCinetSlices[i].getRoi();
							}
							if (indexes.length == 1)
								roiToMeasure = rm.getRoi(0);

							rm.reset();
							// Now it is time to measure mean intensity into specific areas (not whole cell)
							// for XProtein
							IJ.resetMinAndMax(chXProteinSlices[i]);
							chXProteinSlices[i].setRoi(roiToMeasure);
							meanXProtein.add(chXProteinSlices[i].getStatistics().mean);
							meanXProteinTotal.add(String.valueOf(chXProteinSlices[i].getStatistics().mean));
							chXProteinSlices[i].setRoi(roiToMeasure);
							// Lower threshold for XProtein
							double th1XProtein = chXProteinSlices[i].getStatistics().mean
									- chXProteinSlices[i].getStatistics().stdDev;
							// Upper threshold for YProtein
							double th2XProtein = chXProteinSlices[i].getStatistics().mean
									+ chXProteinSlices[i].getStatistics().stdDev;
							// Now it is time to measure mean intensity into specific areas (not whole cell)
							// for YProtein
							IJ.resetMinAndMax(chYProteinSlices[i]);
							chYProteinSlices[i].setRoi(roiToMeasure);
							meanYProtein.add(chYProteinSlices[i].getStatistics().mean);
							meanYProteinTotal.add(String.valueOf(chYProteinSlices[i].getStatistics().mean));
							chYProteinSlices[i].setRoi(roiToMeasure);
							// Lower threshold for XProtein
							double th1YProtein = chYProteinSlices[i].getStatistics().mean
									- chYProteinSlices[i].getStatistics().stdDev;
							// Upper threshold for XProtein
							double th2YProtein = chYProteinSlices[i].getStatistics().mean
									+ chYProteinSlices[i].getStatistics().stdDev;

							for (int z = 0; z < roiToMeasure.getContainedPoints().length; z++) {
								// Keeping only those pixels accomplishing threshold conditions for both
								// X-Yproteins
								if ((chXProteinSlices[i].getProcessor().getPixelValue(
										roiToMeasure.getContainedPoints()[z].x,
										roiToMeasure.getContainedPoints()[z].y) >= (Double) th1XProtein
										&& chXProteinSlices[i].getProcessor().getPixelValue(
												roiToMeasure.getContainedPoints()[z].x,
												roiToMeasure.getContainedPoints()[z].y) <= (Double) th2XProtein)
										&& (chYProteinSlices[i].getProcessor().getPixelValue(
												roiToMeasure.getContainedPoints()[z].x,
												roiToMeasure.getContainedPoints()[z].y) >= (Double) th1YProtein
												&& chYProteinSlices[i].getProcessor().getPixelValue(
														roiToMeasure.getContainedPoints()[z].x, roiToMeasure
																.getContainedPoints()[z].y) <= (Double) th2YProtein)) {

									pixelValuesXProtein0.add((double) chXProteinSlices[i].getProcessor().getPixelValue(
											roiToMeasure.getContainedPoints()[z].x,
											roiToMeasure.getContainedPoints()[z].y));
									pixelValuesYProtein0.add((double) chYProteinSlices[i].getProcessor().getPixelValue(

											roiToMeasure.getContainedPoints()[z].x,
											roiToMeasure.getContainedPoints()[z].y));

								}

							}
							// Here we calculate the mean intensity for each X,Y protein
							Double avgXProtein = (Double) pixelValuesXProtein0.stream().mapToDouble(d -> d).average()
									.orElse(0.0);
							Double avgYProtein = (Double) pixelValuesYProtein0.stream().mapToDouble(d -> d).average()
									.orElse(0.0);
							// Here we store the normalized intensity values in pixelValuesXProtein!
							for (int z = 0; z < pixelValuesXProtein0.size(); z++) {
								pixelValuesXProtein1.add(pixelValuesXProtein0.get(z) - avgXProtein);
								pixelValuesYProtein1.add(pixelValuesYProtein0.get(z) - avgYProtein);
							}
							// Keeping pixels accomplishing conditions established by Manders
							for (int z = 0; z < pixelValuesYProtein0.size(); z++) {
								if ((pixelValuesXProtein0.get(z) > 0 && pixelValuesYProtein0.get(z) > 0)
										|| (pixelValuesXProtein0.get(z) == 0 && pixelValuesYProtein0.get(z) == 0)) {
									pixelValuesYProteinM1.add(pixelValuesYProtein0.get(z));
									pixelValuesYProteinM1Total.add(pixelValuesYProtein0.get(z).doubleValue());
								}

								if ((pixelValuesYProtein0.get(z) > 0 && pixelValuesXProtein0.get(z) > 0)
										|| (pixelValuesYProtein0.get(z) == 0 && pixelValuesXProtein0.get(z) == 0)) {
									pixelValuesXProteinM2.add(pixelValuesXProtein0.get(z));
									pixelValuesXProteinM2Total.add(pixelValuesXProtein0.get(z).doubleValue());
								}

							}
							// Calculating Manders (M1 and M2) coefficients for X and YProtein
							double m1 = (pixelValuesYProteinM1.stream().mapToDouble(d -> d).sum())
									/ (pixelValuesYProtein0.stream().mapToDouble(d -> d).sum());
							double m2 = (pixelValuesXProteinM2.stream().mapToDouble(d -> d).sum())
									/ (pixelValuesXProtein0.stream().mapToDouble(d -> d).sum());

							for (int z = 0; z < pixelValuesXProtein1.size(); z++)
								pixelXProteinYProteinMult.add(pixelValuesXProtein1.get(z).doubleValue()
										* pixelValuesYProtein1.get(z).doubleValue());
							// Now XProtein channel is flipped to calculate Kolmogorov Smirnov with mirror
							// version on XProtein
							IJ.run(chXProteinSlices[i], "Flip Horizontally", "");
							IJ.run(chXProteinSlices[i], "Flip Vertically", "");
							// Here we create the list to store the pixel intensity values for mirror
							// version within
							// specific area for both proteins
							List<Double> pixelValuesXProteinMirror0 = new ArrayList<Double>();
							List<Double> pixelValuesYProteinMirror0 = new ArrayList<Double>();
							List<Double> pixelValuesXProteinMirror1 = new ArrayList<Double>();
							List<Double> pixelValuesYProteinMirror1 = new ArrayList<Double>();

							// Here we iterate through pixels in specific mask area to get pixel values in
							// mirror version
							for (int z = 0; z < roiToMeasure.getContainedPoints().length; z++) {

								pixelValuesXProteinMirror0.add((double) chXProteinSlices[i].getProcessor()
										.getPixelValue(roiToMeasure.getContainedPoints()[z].x,
												roiToMeasure.getContainedPoints()[z].y));
								pixelValuesYProteinMirror0.add((double) chYProteinSlices[i].getProcessor()
										.getPixelValue(roiToMeasure.getContainedPoints()[z].x,
												roiToMeasure.getContainedPoints()[z].y));

							}
							// Then mean value for intensity is calculated in specific area in mirror
							// version
							Double avgXProteinMirror = (Double) pixelValuesXProteinMirror0.stream().mapToDouble(d -> d)
									.average().orElse(0.0);
							Double avgYProteinMirror = (Double) pixelValuesYProteinMirror0.stream().mapToDouble(d -> d)
									.average().orElse(0.0);
							// We normalize pixel intensity values for mirror version
							for (int z = 0; z < pixelValuesXProteinMirror0.size(); z++) {
								pixelValuesXProteinMirror1.add(pixelValuesXProteinMirror0.get(z) - avgXProteinMirror);

							}
							for (int z = 0; z < pixelValuesYProteinMirror0.size(); z++) {
								pixelValuesYProteinMirror1.add(pixelValuesYProteinMirror0.get(z) - avgYProteinMirror);

							}
							// Here we multiply both pixel intensity populations
							for (int z = 0; z < pixelValuesXProteinMirror1.size(); z++)
								pixelXProteinYProteinMultMirror.add(pixelValuesXProteinMirror1.get(z).doubleValue()
										* pixelValuesYProteinMirror1.get(z).doubleValue());
							// Arrays were created to store values from the pixelXProteinYProteinMul list
							// and the
							// corresponding mirror version
							double[] pixelXProteinYProteinMultMirrorArray = new double[pixelXProteinYProteinMultMirror
									.size()];
							double[] pixelXProteinYProteinMultArray = new double[pixelXProteinYProteinMult.size()];

							for (int c = 0; c < pixelXProteinYProteinMultMirrorArray.length; c++) {

								pixelXProteinYProteinMultMirrorArray[c] = pixelXProteinYProteinMultMirror.get(c)
										.doubleValue();
								pixelXProteinYProteinMultMirrorTotal
										.add(pixelXProteinYProteinMultMirror.get(c).doubleValue());
							}
							for (int c = 0; c < pixelXProteinYProteinMultArray.length; c++) {
								pixelXProteinYProteinMultArray[c] = pixelXProteinYProteinMult.get(c).doubleValue();
								pixelXProteinYProteinMultTotal.add(pixelXProteinYProteinMult.get(c).doubleValue());

							}
							// Arrays were created to store initial pixel intensity (not normalized) values
							// from list
							double[] pixelXProteinArray = new double[pixelValuesXProtein0.size()];
							double[] pixelYProteinArray = new double[pixelValuesYProtein0.size()];

							for (int z = 0; z < pixelValuesXProtein0.size(); z++) {
								pixelXProteinArray[z] = pixelValuesXProtein0.get(z);
								pixelXProteinTotal.add(pixelValuesXProtein0.get(z).doubleValue());

							}

							for (int z = 0; z < pixelValuesYProtein0.size(); z++) {
								pixelYProteinArray[z] = pixelValuesYProtein0.get(z);
								pixelYProteinTotal.add(pixelValuesYProtein0.get(z).doubleValue());
							}
							// KolmogorovSmirnovTest class is defined
							testKS = new KolmogorovSmirnovTest();
							// Pearson Correlation class is defined
							pc = new PearsonsCorrelation();
							// p-value from Kolmogorov Smirnow is calculated
							double pValueFil = testKS.kolmogorovSmirnovTest(pixelXProteinYProteinMultArray,
									pixelXProteinYProteinMultMirrorArray);
							double pcCoefficient = pc.correlation(pixelXProteinArray, pixelYProteinArray);

							// Then indicator per slice are stored
							pvaluePerSlice.add(String.valueOf(pValueFil));
							pCorrelationPerSlice.add(String.valueOf(pcCoefficient));
							m1PerSlice.add(String.valueOf(m1));
							m2PerSlice.add(String.valueOf(m2));
						}
						if (roisCinetDef.size() == 0) {
							meanXProteinTotal.add("NaN Value");
							meanYProteinTotal.add("NaN Value");
							pvaluePerSlice.add("NaN Value");
							pCorrelationPerSlice.add("NaN Value");
							m1PerSlice.add("NaN Value");
							m2PerSlice.add("NaN Value");
						}
					}
					// Column names per table were defined
					String[] columnNames = new String[] { "Slice", "XProtein-Mean Intensity", "YProtein-Mean Intensity",
							"Kolmogorov pValue", "Pearson Correlation Coef.", "Manders m1", "Manders m2" };
					// Data per slice is stored on a 2D array
					String[][] rowData = new String[chXProteinSlices.length + 1][columnNames.length];
					for (int z = 0; z < rowData.length - 1; z++) {
						rowData[z][0] = String.valueOf((z + 1));
						rowData[z][1] = String.valueOf(meanXProteinTotal.get(z));
						rowData[z][2] = String.valueOf(meanYProteinTotal.get(z));
						rowData[z][3] = String.valueOf(pvaluePerSlice.get(z));
						rowData[z][4] = String.valueOf(pCorrelationPerSlice.get(z));
						rowData[z][5] = String.valueOf(m1PerSlice.get(z));
						rowData[z][6] = String.valueOf(m2PerSlice.get(z));

					}

//					double[] pixelXProteinYProteinMultMirrorTotalArray = new double[pixelXProteinYProteinMultMirrorTotal.size()];
//					double[] pixelXProteinYProteinMultTotalArray = new double[pixelXProteinYProteinMultTotal.size()];
//
//					for (int c = 0; c < pixelXProteinYProteinMultMirrorTotalArray.length; c++) {
//						pixelXProteinYProteinMultMirrorTotalArray[c] = pixelXProteinYProteinMultMirrorTotal.get(c);
//						pixelXProteinYProteinMultTotalArray[c] = pixelXProteinYProteinMultTotal.get(c);
//					}
//					testKS = new KolmogorovSmirnovTest();
//					pc = new PearsonsCorrelation();
//					double pvalueFilFinal = testKS.kolmogorovSmirnovTest(pixelXProteinYProteinMultMirrorTotalArray,
//							pixelXProteinYProteinMultTotalArray);
//					double[] pixelXProteinArrayTotal = new double[pixelXProteinTotal.size()];
//					double[] pixelYProteinArrayTotal = new double[pixelYProteinTotal.size()];
//
//					for (int z = 0; z < pixelXProteinTotal.size(); z++) {
//						pixelXProteinArrayTotal[z] = pixelXProteinTotal.get(z);
//						pixelYProteinArrayTotal[z] = pixelYProteinTotal.get(z);
//
//					}
//
//					testKS = new KolmogorovSmirnovTest();
//					pc = new PearsonsCorrelation();
//					double pcCoefficientTotal = pc.correlation(pixelXProteinArrayTotal, pixelYProteinArrayTotal);
//					double m1 = (pixelValuesYProteinM1Total.stream().mapToDouble(d -> d).sum())
//							/ (pixelYProteinTotal.stream().mapToDouble(d -> d).sum());
//					double m2 = (pixelValuesXProteinM2Total.stream().mapToDouble(d -> d).sum())
//							/ (pixelXProteinTotal.stream().mapToDouble(d -> d).sum());
//					IJ.log("m1:     " + pixelValuesYProteinM1Total.size() + "---" + pixelYProteinTotal.size() + "------m2:   "
//							+ pixelValuesXProteinM2Total.size() + "----" + pixelXProteinTotal.size());
//
//					List<Double> meanXProteinTotalList = new ArrayList<Double>();
//					List<Double> meanYProteinTotalList = new ArrayList<Double>();
//					for (int y = 0; y < meanXProteinTotal.size(); y++) {
//						if (meanXProteinTotal.get(y) != "NaN Value")
//							meanXProteinTotalList.add(Double.valueOf(meanXProteinTotal.get(y)));
//						if (meanYProteinTotal.get(y) != "NaN Value")
//							meanYProteinTotalList.add(Double.valueOf(meanYProteinTotal.get(y)));
//					}
//
//					double meanXProteinTotalDef = meanXProteinTotalList.stream().mapToDouble(d -> d).average().getAsDouble();
//					double meanYProteinTotalDef = meanYProteinTotalList.stream().mapToDouble(d -> d).average().getAsDouble();
//				
//					rowData[rowData.length - 1][0] = " ";
//					rowData[rowData.length - 1][1] = String.valueOf(meanXProteinTotalDef);
//					rowData[rowData.length - 1][2] = String.valueOf(meanYProteinTotalDef);
//					rowData[rowData.length - 1][3] = String.valueOf(pvalueFilFinal);
//					rowData[rowData.length - 1][4] = String.valueOf(pcCoefficientTotal);
//					rowData[rowData.length - 1][5] = String.valueOf(m1);
//					rowData[rowData.length - 1][6] = String.valueOf(m2);
//				
					ResultsTable results = new ResultsTable(rowData.length);
					for (int i = 0; i < rowData.length; i++)
						for (int j = 0; j < rowData[i].length; j++)
							results.setValue(columnNames[j], i, rowData[i][j]);
					//Results per slice as stores for each image available in the directory.
					//(Here you can type the absolute path in which you want to save each table per image analyzed)
					exportToCSV(new JTable(rowData, columnNames),
							new File("/home/acayuela/data/cris_tirso_agnes/results/results2" + File.separator
									+ "Results_of_" + imps.getShortTitle() + ".csv"));

				}
			}
		});
		thread.start();
	}

	public static ImagePlus[] stack2images(ImagePlus imp) {
		String sLabel = imp.getTitle();
		String sImLabel = "";
		ImageStack stack = imp.getStack();

		int sz = stack.getSize();
		int currentSlice = imp.getCurrentSlice(); // to reset ***

		DecimalFormat df = new DecimalFormat("0000"); // for title
		ImagePlus[] arrayOfImages = new ImagePlus[imp.getStack().getSize()];
		for (int n = 1; n <= sz; ++n) {
			imp.setSlice(n); // activate next slice ***
			ImageProcessor ip = imp.getProcessor();
			ImageProcessor newip = ip.createProcessor(ip.getWidth(), ip.getHeight());
			newip.setPixels(ip.getPixelsCopy());
			sImLabel = imp.getStack().getSliceLabel(n);
			if (sImLabel == null || sImLabel.length() < 1) {
				sImLabel = "slice" + df.format(n) + "_" + sLabel;
			}
			ImagePlus im = new ImagePlus(sImLabel, newip);
			im.setCalibration(imp.getCalibration());
			arrayOfImages[n - 1] = im;

		}
		imp.setSlice(currentSlice);
		if (imp.isProcessor()) {
			ImageProcessor ip = imp.getProcessor();
			ip.setPixels(ip.getPixels());
		}
		imp.setSlice(currentSlice);
		return arrayOfImages;
	}

	public void exportToCSV(JTable table, File file) {
		try {

			TableModel modelo = table.getModel();
			FileWriter excel = new FileWriter(file);
			for (int i = 0; i < modelo.getColumnCount(); i++) {
				excel.write(modelo.getColumnName(i) + ",");
			}
			excel.write("\n");
			for (int i = 0; i < modelo.getRowCount(); i++) {
				for (int j = 0; j < modelo.getColumnCount(); j++) {
					String data = (String) modelo.getValueAt(i, j);
					if (data == "null") {
						data = "";
					}
					excel.write(data + ",");
				}
				excel.write("\n");
			}

			excel.close();

		} catch (IOException ex) {

		}
	}

}
