// This macro automates the analysis pipeline described in Panichnantakul et al (2021) using Comdet and Ezcolocalization for analysis.
// References:
// 1. Panichnantakul et al (2021) DNA Repair 105, 103156
// 2. Comdet. https://github.com/UU-cellbiology/ComDet
// 3. Stauffer, Weston, Huanjie Sheng, and Han N. Lim. "EzColocalization: An ImageJ plugin for visualizing and measuring colocalization in cells and organisms." Scientific reports 8.1 (2018): 15764.


// --- Variables --- //
backgroundImageDapi = "backgroundImage_dapi.tif"
filterRadius = 2
minNucleousSize = 75
maxNucleousSize = 1000
greenChannel = 2 // channel for detecting focci
redChannel = 3 // other channel
focciSize = 6
focciIntensity = 15

// --- Initialization --- //
run("Options...", "iterations=1 count=1 black");
run("Bio-Formats Macro Extensions"); 
setBatchMode(false);
run("ROI Manager..."); // open the ROI manager

// --- Select the Image --- //
#@ File (label = "Select a file", style = "file") file
imageFolder = File.getDirectory(file);
outputFile = "Analysis of " + File.getNameWithoutExtension(file);
outputFileFullPath = imageFolder + File.separator + outputFile;
backgroundImageDapi = imageFolder + backgroundImageDapi;
Ext.setId(file);

// Determine the number of series.
Ext.getSeriesCount(seriesCount);

// Create output tables
outputTable_Nuclei = "output_table_nuclei";
Table.create(outputTable_Nuclei);
var outputLine_nuclei = 0;

outputTableFocciGreenChannel = "output_table_focci_green_channel";
Table.create(outputTableFocciGreenChannel);
var outputLineFocciGreenChannel = 0;

outputTableFocciRedChannel = "output_table_focci_red_channel";
Table.create(outputTableFocciRedChannel);
var outputLineFocciRedChannel = 0;

//seriesNo = 101; 
//analyseImage(seriesNo); // remove for batch
for (seriesNo=1; seriesNo <= seriesCount; seriesNo++) { //uncomment for batch
	print("Analysing Series No "+ seriesNo + "/" + seriesCount);
	analyseImage(seriesNo);
}

// -- Save Output and Exit --- // 
selectWindow(outputTable_Nuclei);
saveAs("Results", outputFileFullPath + "_nuclei.csv");

selectWindow(outputTableFocciGreenChannel);
saveAs("Results", outputFileFullPath + "_focci_green_channel.csv");

selectWindow(outputTableFocciRedChannel);
saveAs("Results", outputFileFullPath + "_focci_red_channel.csv");

setBatchMode(false);

// --- Analyse each image --- //
function analyseImage(seriesNo) { 
	// Open the image and duplicate it
	run("Bio-Formats Importer", "open=" + file + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesNo);
	originalImage = getTitle();
	originalDuplicate = "original_duplicate";
	run("Duplicate...", "title=[" + originalDuplicate + "] duplicate");
	run("Gaussian Blur...", "sigma=1 stack"); //denoising filter
	
	// Starting
	nucleousImage = "nucleous_image";
	run("Duplicate...", "title=" + nucleousImage + " duplicate channels=1");
	
	// run a background subtraction if the background image exists
	if (File.exists(backgroundImageDapi)) {
		subtractBackgroundImage(nucleousImage,backgroundImageDapi);
	}
	run("Gaussian Blur...", "sigma=5");
	setAutoThreshold("Li dark");
	run("Convert to Mask");
	for (i = 0; i < 2; i++) { //expand the selection to crop a larger area
		run("Dilate");
	}
	clearRoiManager();
	run("Analyze Particles...", "size=" + minNucleousSize +"-" + maxNucleousSize + "Infinity exclude clear add");
	
	// Open nuclei images
	nNucleous = roiManager('count');
	for (nucleousId = 0; nucleousId < nNucleous; nucleousId++) {
		selectWindow(originalDuplicate);
	    roiManager('select', nucleousId);
	    run("Duplicate...", "title=nucleous_" + nucleousId + " duplicate");
	}
	
	// Analyse nuclei
	for (nucleousId = 0; nucleousId < nNucleous; nucleousId++) {
		print("Analysing nucleous: "+ nucleousId);
		nFocciGreen = analyseNucleous("nucleous_" + nucleousId, greenChannel, redChannel);
		print("Green channel analysed. Focci found = " + nFocciGreen);
		nfocciRed = analyseNucleous("nucleous_" + nucleousId, redChannel, greenChannel);
		print("Red channel analysed. Focci found = " + nfocciRed);
		
		selectWindow(outputTable_Nuclei);
		Table.set("Image_name", outputLine_nuclei, originalImage);
		Table.set("Cell_no", outputLine_nuclei, nucleousId);
		Table.set("Green_Focci", outputLine_nuclei, nFocciGreen);
		Table.set("Red_Focci", outputLine_nuclei, nfocciRed);
		outputLine_nuclei = outputLine_nuclei + 1;
		
		closeWindow("nucleous_" + nucleousId);
	}
	
	
	// Close windows
	selectWindow(outputTable_Nuclei); //save progress
	saveAs("Results", outputFileFullPath + "_nuclei.csv");
	Table.rename(outputFile + "_nuclei.csv", outputTable_Nuclei);	
	
	selectWindow(outputTableFocciGreenChannel); //save progress
	saveAs("Results", outputFileFullPath + "_focci_green_channel.csv");
	Table.rename(outputFile + "_focci_green_channel.csv", outputTableFocciGreenChannel);
	
	selectWindow(outputTableFocciRedChannel); //save progress
	saveAs("Results", outputFileFullPath + "_focci_red_channel.csv");
	Table.rename(outputFile + "_focci_red_channel.csv", outputTableFocciRedChannel);
	
	closeWindow(nucleousImage);
	closeWindow(originalDuplicate);
	closeWindow(originalImage);
}

// --- Analyse each Nucleous --- //
function analyseNucleous(image, focciChannel, comparisonChannel) { 
	selectWindow(image);
	run("Select None");
	focciChannelWindow = "focci_channel";
	run("Duplicate...", "title=" + focciChannelWindow + " duplicate channels=" + focciChannel);
	clearRoiManager();
	
	// Run ComDet
	run("Detect Particles", "ch1i ch1a=" + focciSize + " ch1s=" + focciIntensity + " rois=Rectangles add=[All detections] summary=Reset");
	closeWindow(focciChannelWindow);
	
	// loop through the focci
	nFocci = roiManager('count');
	for (focciId = 0; focciId < nFocci; focciId++) {
		analyseFocci(image, focciId, focciChannel, comparisonChannel);
	}
	return nFocci;
}

// --- Analyse each focci --- //
function analyseFocci(image, focciId, focciChannel, comparisonChannel) {
		selectWindow(image);
		roiManager("select", focciId);
		focciImage = "focci_" + focciId + "_main";
		run("Duplicate...", "title=" + focciImage + " duplicate channels=" + focciChannel);
		selectWindow(image);
		roiManager("select", focciId);
		comparisonImage = "focci_" + focciId + "_comparison";
		run("Duplicate...", "title=" + comparisonImage + " duplicate channels=" + comparisonChannel);
		run("EzColocalization ", "reporter_1_(ch.1)=" + focciImage + " reporter_2_(ch.2)=" + comparisonImage + " alignthold4=default tos metricthold1=costes' allft-c1-1=10 allft-c2-1=10 pcc metricthold2=all allft-c1-2=10 allft-c2-2=10 srcc metricthold3=all allft-c1-3=10 allft-c2-3=10 icq metricthold4=all allft-c1-4=10 allft-c2-4=10 mcc metricthold5=costes' allft-c1-5=10 allft-c2-5=10");
		analysisTable = "Metric(s) of " + comparisonImage;
		
		// Read the analysis table
		selectWindow(analysisTable);
		tosLinear = Table.get("TOS(linear)", 0);
		tosLog2 = Table.get("TOS(log2)", 0);
		pcc = Table.get("PCC", 0);
		srcc = Table.get("SRCC", 0);
		icq = Table.get("ICQ", 0);
		m1 = Table.get("M1", 0);
		m2 = Table.get("M2", 0);
		
		// Fill the output table
		if (focciChannel == greenChannel) { //Set variables to the green channel
			outputTable = outputTableFocciGreenChannel;
			outputLine = outputLineFocciGreenChannel;
		}
		if (focciChannel == redChannel) { //Set variables to the red channel
			outputTable = outputTableFocciRedChannel;
			outputLine = outputLineFocciRedChannel;
		}
		
		selectWindow(outputTable);
		Table.set("Image_name", outputLine, originalImage);
		Table.set("Cell_no", outputLine, nucleousId);
		Table.set("Focci_no", outputLine, focciId);
		Table.set("TOS(linear)", outputLine, tosLinear);
		Table.set("TOS(log2)", outputLine, tosLog2);
		Table.set("PCC", outputLine, pcc);
		Table.set("SRCC", outputLine, srcc);
		Table.set("ICQ", outputLine, icq);
		Table.set("M1", outputLine, m1);
		Table.set("M2", outputLine, m2);
		
		if (focciChannel == greenChannel) {
			outputLineFocciGreenChannel = outputLine + 1;
		}
		if (focciChannel == redChannel) {
			outputLineFocciRedChannel = outputLine + 1;
		}
		
		// Close windows
		closeWindow(focciImage);
		closeWindow(comparisonImage);
		closeTable(analysisTable);
		closeTable("Summary");
		
}

function subtractBackgroundImage(image,backgroundImagePath) {
	open(backgroundImagePath);
	backgroundImage = getTitle();
	imageCalculator("Divide create 32-bit", image, backgroundImage);
	rename("temp_image");
	closeWindow(image);
	closeWindow(backgroundImage);
	selectWindow("temp_image");
	rename(image);
}

function closeWindow(image){
	selectWindow(image);
	wait(50);
	close();
	wait(50);
}

function closeTable(table){
	if (isOpen(table)) {
		selectWindow(table);
		run("Close");
	}
}

function clearRoiManager() {
	if (roiManager("count") > 0) {
		roiManager("deselect");
		roiManager("delete");
	}
}
