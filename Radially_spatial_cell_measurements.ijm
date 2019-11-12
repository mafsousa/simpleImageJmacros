/* Title: Radially spatial cell measurements
 * Version:v0.2
*
* Short description: Input stack files with CY3 in the first channel (result from macro FromInCellToHyperstack) 
* and the corresponding cell contours, nucleus contour and endossomes contours (result from Cellprofiler). 
* The idea is to create 2 bands from the original contour with "band_with" specified by the user. 
* Measure area and intensity inside each band, for each cell (over the endossomes region) and save it as an excel file
* Requirements: install IJPB-plugin (https://github.com/ijpb/MorphoLibJ) for morphological watershed
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

#@ File (label = "Select input original images directory", style = "directory") input
#@ File (label = "Select overlay image ", style = "directory") overlays
#@ File (label = "Select output directory", style = "directory") output
#@ Integer (label = "Band width", value = "4") band_width
#@ String (label = "File suffix", value = ".tif") suffix


setBatchMode(false);
run("Close All");
roiManager("reset");

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	//process for each file
	print("Processing: " + input + File.separator + file);
	roiManager("reset");
	
	//open original stack (1st channel is cy3)
	open(input + File.separator + file);
	original_stack = getTitle();
	//get only Cy3 channel
	run("Duplicate...", "title=Cy3");	
	run("8-bit");

	//get overlay name from the original name
	dotIndex = lastIndexOf( original_stack, "wv");
	title_begin = substring( original_stack , 0, dotIndex);

	//Note: predifined file names from cellprofiler output
	// edit the following lines with particular file names
	cell_overlay = title_begin + "wv FITC - FITC)CellsOutlines.tiff";	
	nuc_overlay = title_begin + "wv FITC - FITC)NucsOutlines.tiff";	
	endo_overlay = title_begin + "wv FITC - FITC)EndoOutlines.tiff";	
		
	//close stack
	selectWindow(original_stack);
	close();

	//open overlay images
	print(overlays + File.separator + cell_overlay);
	open(overlays + File.separator + cell_overlay);
	cells = getTitle();
	open(overlays + File.separator + nuc_overlay);
	nucs = getTitle();
	open(overlays + File.separator + endo_overlay);
	endo = getTitle();

	//identify correct contours based on cell nucleus
	selectWindow(nucs);
	run("Invert");
	run("Fill Holes");
	run("Find Maxima...", "prominence=1 strict exclude light output=[Single Points]");
	rename("marker");
	run("Dilate");
	selectWindow(cells);
	run("Invert");
	run("Fill Holes");
	//watershed MORPHOLIBJ
	run("Marker-controlled Watershed", "input=" + cells +" marker=marker mask=None binary calculate use");
	rename("watershed_lines");
	run("8-bit");
	setAutoThreshold("Default dark");
	setThreshold(0, 0);
	run("Convert to Mask");
	selectWindow(cells);
	run("Invert");
	imageCalculator("Subtract create", cells,"watershed_lines");

	//Invert since overlay image is 0 value for bands
	run("Invert");

	//get each contour, exclude on the edge of the image and add to RoiManager
	run("Analyze Particles...", "size=3000-Infinity show=Masks exclude add display");
	saveAs("Tiff", output + File.separator + title_begin + "_cells");

	//Perform roi manipulations to get the 3 bands (enlarge -, XOR, twice)
	number_cells = roiManager("count"); 
	print(number_cells);
	roi1 = newArray(number_cells);
	roi2 = newArray(number_cells);
	roi3 = newArray(number_cells);
	
	for (r = 0; r < number_cells; r++) {
		indx = 4*r +number_cells;		
		roiManager("select", r);
		run("Enlarge...", "enlarge=-" + band_width); //reduce roi		
		IJ.redirectErrorMessages();	
		roiManager("Add");	
		wait(200);				
		if (findSameRoi(r,indx) == true ) {
			print("1 - Unable to shrink roi " + r );
			roiManager("select", r);
			roiManager("Add");  //had just to keep indexes formula
			roiManager("Add");  	
			roiManager("Add");  				
  			roi1[r] = r;  //outside band
  			roi2[r] = -1;	//intermediate band
			roi3[r] = -1;	//inside band	
			continue;
		} else {
  			roiManager("select", newArray(r,indx));	
			roiManager("XOR");
			roiManager("Add");
			roi1[r] = indx+1;  //outside band			
		}			
		roiManager("select", indx);
		run("Enlarge...", "enlarge=-" + band_width);
		IJ.redirectErrorMessages(); //to avoid stoping in case there is no ROI to add
		roiManager("Add");		
		wait(200);	
		if (findSameRoi(indx,indx+2) == true ) {
			print("2 -Unable to shrink roi " + r );
			roiManager("select", indx+1);
			roiManager("Add");//had just to keep indexes formula
  			roi2[r] = -1;	//intermediate band
			roi3[r] = indx + 2;	//inside band	
			continue;	
		} else {
  			roiManager("select", newArray(indx,indx+2));	
			roiManager("XOR");
			roiManager("Add");
			roi2[r] = indx+3;	//intermediate band
			roi3[r] = indx+2; //inside band	
		}		
	}

	selectWindow(endo);
	//preprocess endossome countours
	run("Multiply...","value=255");
	run("Invert");
	run("Divide...","value=255");
	//Get Cy3 only inside endossomes
	imageCalculator("Multiply create",endo,"Cy3");
	result_endo_in_cy3 =getTitle();
	
	//create Results table to store the measurements
	title1 = "Results table";
	title2 = "["+title1+"]";
	f=title2;
	run("New... ", "name="+title2+" type=Table");
	print(f,"\\Headings:Cell\tArea_Roi1\tMean_Roi1\tArea_Roi2\tMean_Roi2\tArea_Roi3\tMean_Roi3");

	//Measure the three bands for each cell
	selectWindow(result_endo_in_cy3);

	for (s = 0; s < number_cells; s++) {
		if(roi2[s]==-1 && roi3[s]==-1){
			roiManager("select", newArray(roi1[s]));
		}
		else if(roi2[s]==-1 && roi3[s] !=-1){
			print(roi1[s],roi3[s]);
			roiManager("select", newArray(roi1[s],roi3[s]));
		}
		else if(roi2[s]!=-1 && roi3[s] !=-1){
			print(roi1[s],roi2[s],roi3[s]);
			roiManager("select", newArray(roi1[s],roi2[s],roi3[s]));	
		}		
		roiManager("Multi Measure");
		headings = split(String.getResultsHeadings);
		line = "";
		for (a=0; a<lengthOf(headings); a++)
    		line = line + getResult(headings[a],0) + "\t";
		
		array_line = split(line, "\t");		
		
		//change the order of rois to be: outside, intermediate, inside 
		if(roi2[s]==-1 && roi3[s]==-1){
			print(f,s+"\t"+array_line[0]+"\t"+array_line[1]+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "); //inside 	
		}
		else if(roi2[s]==-1 && roi3[s] !=-1){
			print(f,s+"\t"+array_line[0]+"\t"+array_line[1]+"\t"+array_line[2]+"\t"+array_line[3]+"\t"+" "+"\t"+" "); //change the order of rois to be: outside,  inside 	
		}
		else if(roi2[s]!=-1 && roi3[s] !=-1){
			print(f,s+"\t"+array_line[0]+"\t"+array_line[1]+"\t"+array_line[4]+"\t"+array_line[5]+"\t"+array_line[2]+"\t"+array_line[3]); //change the order of rois to be: outside, intermediate, inside 	
		}
	}

	//save output image with bands and results table
	roiManager("show all");
	saveAs("tiff",output + File.separator + title_begin );	
	run("Close All");
	selectWindow("Results table");
	saveAs("Text", output + File.separator + title_begin + ".xls");
	close("Results table");
}

function findSameRoi(ind1, ind2) {  
 	roiManager("Select", ind1);
 	Roi.getBounds(x1, y1, width, height);
	roiManager("Select", ind2);
 	Roi.getBounds(x2, y2, width, height);
	if (x1==x2 && y1==y2) { 
			return true; 
	} 
	return false; 
} 
