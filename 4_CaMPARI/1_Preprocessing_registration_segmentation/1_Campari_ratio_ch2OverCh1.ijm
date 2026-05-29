
//path = getDirectory("Choose a directory");
//setBatchMode(true);

path = "Path_to_images";


processFolder(path);

// Function to scan folders/subfolders/files to find files with correct suffix
function processFolder(path) {
    list = getFileList(path);
    print("List length is " + list.length);
    for (j = 0; j < list.length; j++) {
        // Construct full path for the current item
        fullPath = path + "/" + list[j];
        
        // Check if the current item is a directory
        if (File.isDirectory(fullPath)) {
            print("Processing folder: " + fullPath);
            processFolder(fullPath); // Recursive call to process subfolder
        } else if (endsWith(list[j], ".czi") || endsWith(list[j], ".lsm") || endsWith(list[j], ".tiff") || 
                   endsWith(list[j], ".nd2") || endsWith(list[j], ".vsi")) {
            // Process files with specified extensions
            processFile(path, list[j]);
        }
    }
}
function processFile(path, listname) {
filename=listname;
	    // Construct the full path to the file
	    fullFilePath = path + listname;
	
	    // Construct the output file path
	    outputFilePath = path + "processed_" + listname;
	
	    // Check if the output file already exists
	    if (File.exists(outputFilePath)) {
	        print("File already exists, skipping: " + outputFilePath);
	        return; // Skip this file
	    }
	
	    // Proceed with processing the file
	    print("Processing file: " + fullFilePath);
	    run("Close All");
	
	    // Open the file using Bio-Formats
	    run("Bio-Formats Windowless Importer", "open=[" + fullFilePath + "]");
	    orig_name = getTitle();
	    print("Original file name: " + orig_name);
	
	    num_slices = nSlices;
		//run("Size...", "width=512 height=512 depth="+num_slices+" constrain average interpolation=Bilinear");
		orig_id = getImageID();
		run("Duplicate...", "duplicate channels=1");
		rename("ch1");
		ch1_name = getTitle();
		ch1_id = getImageID();
		selectImage(orig_id);
		run("Duplicate...", "duplicate channels=2");
		rename("ch2");
		ch2_name = getTitle();
		ch2_id = getImageID();
		imageCalculator("Divide create 32-bit stack", ch2_name,ch1_name);
		
		rename("ratio");
		ch3_name = getTitle();
		setSlice(num_slices/2);
		run("Enhance Contrast", "saturated=0.35");
		run("16-bit");
		run("Fire");
		ch3_id = getImageID();
		//run("Size...", "width=512 height=512 depth="+num_slices+" constrain average interpolation=Bilinear");
		run("Merge Channels...", "c1="+ch1_name+" c2="+ch2_name+" c3="+ch3_name+" create");
		//new_name=orig_name+"ratio";

		//print(saving_path);
		//replace(saving_path, "*\*" , "*/*");
		//print(saving_path);
	    // Save the processed file
	    saveAs("Tiff", outputFilePath);
	    print("Saved processed file: " + outputFilePath);
			 
	
}



run("Close All");
run("Collect Garbage");



//setBatchMode(false);
print("done");


