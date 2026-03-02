// Define the input directory and the output directory
inputDir = getDirectory("Choose a Directory");
outputDir = inputDir + "merged_files/";
setBatchMode(true);{
// Create the output directory if it doesn't exist
if (!File.exists(outputDir)) {
    File.makeDirectory(outputDir);
}

// Get a list of all the files in the input directory
list = getFileList(inputDir);

// Loop through each file in the input directory
for (i = 0; i < list.length; i++) {
    filename = list[i];

    // Check if the file is an NRRD file
    if (endsWith(filename, ".nrrd")) {

        // Get the source and channel from the filename
        parts = split(filename, "_");
        source = parts[2];
        channel = parts[3];

        // Check if the file is the first channel for the source
        if (channel == "01") {

            // Open the three channels for the source as stacks
            ch1 = inputDir +parts[0]+"_"+parts[1]+"_" + source + "_01_"+parts[4]+"_"+parts[5];
            print(ch1);
            ch2 = inputDir +parts[0]+"_"+parts[1]+"_" + source + "_02_"+parts[4]+"_"+parts[5];
            print(ch2);
            ch3 = inputDir +parts[0]+"_"+parts[1]+"_" + source + "_03_"+parts[4]+"_"+parts[5];
            print(ch3);
            //run("Image Sequence...", "open=[" + ch1 + "] sort use");
            //run("Image Sequence...", "open=[" + ch2 + "] sort use");
            //run("Image Sequence...", "open=[" + ch3 + "] sort use");
			open(ch1);
			title1 = getTitle();
			open(ch2);
            title2 = getTitle();
			open(ch3);
            title3 = getTitle();
            // Merge the hyperstacks into a multichannel hyperstack
     //       run("Merge Channels...", "c1=[Stack 1] c2=[Stack 2] c3=[Stack 3] create");
run("Merge Channels...", "c2="+title2+" c6="+title1+" c7="+title3+" create keep");
saveAs("tiff", outputDir+title1+"_merged");
            // Save the multichannel hyperstack as a TIFF file in the output directory
           // outputFilename = "merged_" + source + ".tif";
          //  saveAs("Tiff", outputDir + outputFilename);

            // Close all the images
            run("Close All");
            run("Collect Garbage");
        }
    }
}
}
print("done");