/* Lif to Tiff series
 * v0.2
*
* Read Lif file from input folder and save each serie (Position) as tiff file in the output folder (Batch mode)
* 
* This macro should NOT be redistributed without author's permission. 
* Explicit acknowledgement to the ALM facility should be done in case of published articles (approved in C.E. 7/17/2017):     
* 
* "The authors acknowledge the support of i3S Scientific Platform Advanced Light Microscopy, 
* member of the national infrastructure PPBI-Portuguese Platform of BioImaging (supported by POCI-01-0145-FEDER-022122)."
* 
* Date: July/2018
* Author: Mafalda Sousa, mafsousa@ibmc.up.pt 
* Advanced Ligth Microscopy, I3S 
* PPBI-Portuguese Platform of BioImaging
* */
 
#@ File (label = "Select input original images directory", style = "directory") inDir
#@ File (label = "Select overlay image ", style = "directory") overlays
#@ File (label = "Select uutput directory", style = "directory") outDir
#@ Integer (label = "Band width", value = "4") band_width
#@ String (label = "File suffix", value = ".lif") ext


filelist = getFileList(inDir); //load array of all files inside input directory
filename = filelist[0];
print("Processing " + filename + " conversion");     
setBatchMode(true);
for (i=1; i<= filelist.length; i++) {
 
     end = lastIndexOf(filename,".");
     series_name = "series_" + i;
     print(series_name);
    
 	if (i < 10){  
        tt = "Position00" + i;
 	}
 	else if (i>=100) {
 		 tt = "Position" + i;
 	}
 	else {
 		 tt = "Position0" + i;
 	}

     // open file, requires LOCI tools (aka Bio-Formats)
     run("Bio-Formats Importer", "open="+ inDir + filename +  " color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack " + series_name);
     print(i);
     saveAs("Tiff", outDir + tt + ".tiff");
     
}

setBatchMode(false);
print("-- Done --");