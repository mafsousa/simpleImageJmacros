/* Title: Spectral FRET analysis
 * Version:
*
* Short description: A specific macro to analyze particular spectral fret images. Inputs are Donnor only, Acceptor only and FRET images.  
* Application of "Gustina and Trudeau 10.1073/pnas.0900180106" formulas
* 
* This macro should NOT be redistributed without author's permission. 
* Explicit acknowledgement to the ALM facility should be done in case of published articles (approved in C.E. 7/17/2017):     
* 
* "The authors acknowledge the support of i3S Scientific Platform Advanced Light Microscopy, 
* member of the national infrastructure PPBI-Portuguese Platform of BioImaging (supported by POCI-01-0145-FEDER-022122)."
* 
* Date: November/2019
* Author: Mafalda Sousa, mafsousa@ibmc.up.pt 
* Advanced Ligth Microscopy, I3S 
* PPBI-Portuguese Platform of BioImaging
*/



#@ File (label = "Input Donnor only directory", style = "directory") input_donnor
#@ File (label = "Input Acceptor only directory", style = "directory") input_acceptor
#@ File (label = "Input FRET directory", style = "directory") input_fret
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ Double(label = "Band length (um)", value = 0.5) band

run("Clear Results");
run("Close All");
roiManager("reset");

processFolder(input_donnor);

list = getFileList(input_donnor);
donnor_average = newArray(19);
for (aux = 0; aux < 19; aux++) {
	donnor_average[aux] = (getResult("A0", aux) +  getResult("A1", aux) +getResult("A2", aux) +getResult("A3", aux))/4;
}
print("Donnor average");
Array.print(donnor_average);
run("Clear Results");
roiManager("reset");

Acc_Ratio = processFolder2(input_acceptor); 
print("Acceptor ratio");
Array.print(Acc_Ratio);

processFolder3(input_fret,donnor_average);
selectWindow("Log");
saveAs("text", output + File.separator + "Summary results_" + band + ".txt");

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	donnor_index = 0;
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)){
			A = processFile(input, output, list[i],donnor_index);
			if (donnor_index == 0)	A0 = A;
			if (donnor_index == 1)	A1 = A;
			if (donnor_index == 2)	A2 = A;
			if (donnor_index == 3)	A3 = A;
			//print in table each column
		
			}			
			donnor_index = donnor_index + 1;			
		}
	Array.show("Table (indexes)", A0, A1, A2, A3);
	selectWindow("Table");
	saveAs("text", output + File.separator + "Table.csv");
	close();

	run("Clear Results");
	table_file = File.openAsString(output + File.separator + "Table.csv");
	ImportResultsTable(table_file);
}

function processFile(input, output, file,donnor_index) {
	roiManager("reset");
	A = process_donnor_only(input,file,donnor_index);
	return A;
}


// function to scan folders/subfolders/files to find files with correct suffix
function processFolder2(input) {
	list = getFileList(input);
	list = Array.sort(list);
	id = 0;
	
	Acceptor_ratio = newArray(list.length/2);
	for (i = 0; i < list.length; i=i+2) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder2(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)){
			Acceptor_ratio[id] = processFile2(input, output, list[i]);	
			id = id +1;		
		}			
	}
	return Acceptor_ratio;
}

function processFile2(input, output, file) {
	open(input + File.separator + file);
	job1 = getTitle();
	job2 = replace(input + File.separator + file, "Job1", "Job2");
	open(job2);
	job2 = getTitle();
	//create Results table to store the measurements
	
	A = process_acceptor_only(job1,job2);
	return A;

}

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder3(input,donnor_average) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i=i+2) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder3(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)){
			processFile3(input, list[i], donnor_average);	
			}			
	}
}

function processFile3(input, file,donnor_average) {
	open(input + File.separator + file);
	job1 = getTitle();
	job2 = replace(input + File.separator + file, "Job1", "Job2");
	open(job2);
	job2 = getTitle();
	print(job1);
	FRET_RatioA(job1,job2, donnor_average);

}

//process_acceptor_only(acceptor_only);
//donnor_image=getTitle();
//process_donnor_only(donnor_image);

function process_donnor_only(input,input,idx){
	
	open(input + File.separator + file);
	run("Clear Results");
	run("Duplicate...", " duplicate");
	original = getTitle();
	getDimensions(width, height, channels, slices, frames);
	
	
	run("Enhance Contrast", "saturated=0.35");
	/*----------------Cell segmentation---------------------------*/
	run("Z Project...", "projection=[Max Intensity]");
	run("Enhance Contrast", "saturated=0.35");
	//run("Extended Min & Max", "operation=[Extended Minima] dynamic=16 connectivity=4");
	//remove background
	//run("AutoCutB", "method=Mode extend=1 subtract");
	setAutoThreshold("Default dark");
	run("Threshold...");
	//waitForUser("Select best threshold");
	
	run("Options...", "iterations=1 count=1 black do=Close");
	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	//setThreshold(65, 255);
	run("Convert to Mask");
	
	run("Analyze Particles...", "size=50-infinity show=Outlines add");
	//TODO: confirm one cell only
	selectWindow(original);
	
	index = createband();

	roiManager("select", index);		
	run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
	roiManager("Multi Measure append");
	roiManager("select", 0);
	run("Clear Outside", "stack");
	roiManager("Select", index);
	run("Add Selection...");
	save(output + File.separator + original + "roi.roi");
	run("royal");
	run("Flatten", "stack");
	saveAs("Tiff", output + File.separator + "Images/"+ original +"_"+ band + "_band");	

	donnor_only = newArray(frames);
	donnor_only_norm = newArray(frames);
	for (r = 0;r<nResults;r++){
		donnor_only[r]=getResult("Mean1", r);
	}

	Array.getStatistics(donnor_only, min, donnor_max, donnor_averg, stdDev);
	
	for (r = 0;r<19;r++){
		donnor_only_norm[r]=donnor_only[r]/donnor_max;
	}
	
	
	return donnor_only_norm;
	
	
}

function process_acceptor_only(job1,job2){
	
	roiManager("reset");
	
	selectWindow(job2);
	run("Duplicate...", " duplicate");
	job2_proj = getTitle();
	getDimensions(width, height, channels, slices, frames);	
	run("Enhance Contrast", "saturated=0.35");
	/*----------------Cell segmentation---------------------------*/
	run("Z Project...", "projection=[Max Intensity]");
	run("Enhance Contrast", "saturated=0.35");
	setAutoThreshold("Li dark");
	run("Threshold...");
	//waitForUser("Select best threshold");
	
	run("Options...", "iterations=1 count=1 black do=Close");
	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Li dark");
	//run("Threshold...");
	//setThreshold(65, 255);
	run("Convert to Mask");
	
	run("Analyze Particles...", "size=50-infinity show=Outlines add");
	//TODO: confirm one cell only
	selectWindow(job1);
	index = createband();
	roiManager("select", index);		
	run("Set Measurements...", " mean redirect=None decimal=3");
	roiManager("Multi Measure ");
	//sort select last with max value and divide all mean values
	
	acceptor_j1 = newArray(19);
	acceptor_j2 = newArray(14);
	for (r = 0;r<nResults;r++){
		acceptor_j1[r] = getResult("Mean1", r);
	}

	run("Clear Results");
	selectWindow(job2);
	
	roiManager("Select", index);
	run("Set Measurements...", " mean redirect=None decimal=3");
	roiManager("Multi Measure ");	
	roiManager("select", 0);
	run("Clear Outside", "stack");
	roiManager("Select", index);
	run("Add Selection...");
	save(output + File.separator + job1 + "roi.roi");
	run("royal");
	run("Flatten", "stack");
	
	saveAs("Tiff", output + File.separator + "Images/"+ job2 +"_"+ band + "_band");	

	
	run("Close All");
	for (r = 0;r<14;r++){
		acceptor_j2[r] = getResult("Mean1", r);
	}

	RatioA0 = newArray(5);
	for (i = 0; i < 5; i++) {
		RatioA0[i] = acceptor_j1[9+i]/acceptor_j2[4+i];
	}

	Array.getStatistics(RatioA0, min, max, Ratio_A0, Ratio_A0_std);
	return Ratio_A0;	
}

function FRET_RatioA(job1,job2, donnor_average){
	roiManager("reset");

	selectWindow(job2);
	run("Duplicate...", " duplicate");
	job2_proj = getTitle();
	getDimensions(width, height, channels, slices, frames);	
	run("Enhance Contrast", "saturated=0.35");
	/*----------------Cell segmentation---------------------------*/
	run("Z Project...", "projection=[Max Intensity]");
	run("Enhance Contrast", "saturated=0.35");
	//run("Extended Min & Max", "operation=[Extended Minima] dynamic=16 connectivity=4");
	//remove background
	//run("AutoCutB", "method=Mode extend=1 subtract");
	setAutoThreshold("Li dark");
	run("Threshold...");
	
	run("Options...", "iterations=1 count=1 black do=Close");
	//waitForUser(" ");

	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Li dark");
	run("Convert to Mask");
	
	run("Analyze Particles...", "size=50-infinity show=Outlines add");
	//TODO: confirm one cell only
	selectWindow(job1);
	
	index = createband();
	roiManager("select", index);		
	run("Set Measurements...", " mean redirect=None decimal=3");
	roiManager("Multi Measure ");
	//sort select last with max value and divide all mean values
	
	fret_j1 = newArray(19);
	fret_j2 = newArray(14);
	for (r = 0;r<nResults;r++){
		fret_j1[r] = getResult("Mean1", r);
	}
	
	Array.getStatistics(fret_j1, min, FRET1_max, avg, std);
	run("Clear Results");
	print("Maximum CFP458_FRET: ", FRET1_max);

	selectWindow(job2);	
	roiManager("select", index);

	run("Set Measurements...", " mean redirect=None decimal=3");
	roiManager("Multi Measure ");	

	//save roi and image
	roiManager("select", 0);
	run("Clear Outside", "stack");
	roiManager("Select", index);
	run("Add Selection...");
	save(output + File.separator + job1 + "roi.roi");
	run("royal");
	run("Flatten", "stack");
	saveAs("Tiff", output + File.separator + "Images/"+ job2 + "_"+ band + "_band");	
	run("Close All");
	roiManager("reset");

	// calculate RatioA
	for (r = 0;r<14;r++){
		fret_j2[r] = getResult("Mean1", r);
	}
	
	FRET_norm = newArray(19);
	F458_tot = newArray(19);

	for (i = 0; i < 19; i++) {
		FRET_norm[i] = FRET1_max*donnor_average[i];
		F458_tot[i] = fret_j1[i]-FRET_norm[i];
	}
	print("CFP458_FRET norm:");
	Array.print(FRET_norm);
	print("F458_total:");
	Array.print(F458_tot);
	RatioA = newArray(5);
	for (i = 0; i < 5; i++) {
		RatioA[i] = F458_tot[9+i]/fret_j2[4+i];
		
	}
	print("Ratio A vector:");
	Array.print(RatioA);
	Array.getStatistics(RatioA, min, max, Ratio_A, Ratio_A_std);

	print("Ratio A of " + job1 + " = ",Ratio_A);
	
}

function ImportResultsTable(table_file){
     requires("1.35r");
     lineseparator = "\n";
     cellseparator = ",\t";

     // copies the whole RT to an array of lines
     lines=split(table_file, lineseparator);

     // recreates the columns headers
     labels=split(lines[0], cellseparator);
     if (labels[0]==" ")
        k=1; // it is an ImageJ Results table, skip first column
     else
        k=0; // it is not a Results table, load all columns
     for (j=k; j<labels.length; j++)
        setResult(labels[j],0,0);

     // dispatches the data into the new RT
     run("Clear Results");
     for (i=1; i<lines.length; i++) {
        items=split(lines[i], cellseparator);
        for (j=k; j<items.length; j++)
           setResult(labels[j],i-1,items[j]);
     }
     updateResults();
 }

 function createband(){
 	
 	roiManager("Select", 0);
 	if (band == 0){
 		indx = 0;
 		return 	indx;
 	}
 	else{
		run("Enlarge...", "enlarge=-" + band);
		roiManager("add");
		roiManager("Select", newArray(0,1));
		roiManager("XOR");
		roiManager("add");
		indx = 2;
		return indx;
 	}
	
 }