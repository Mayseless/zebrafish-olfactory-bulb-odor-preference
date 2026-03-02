



// function description

path = getDirectory("Choose a directory");
	filename = getFileList(path);
	resDir = path + "nrrd_files" + File.separator;
	File.makeDirectory(resDir);
	run("Close All");
setBatchMode(true);{
for (i=0; i<filename.length; i++) {
	if(endsWith(filename[i], ".czi")|| (endsWith(filename[i], ".tif"))|| (endsWith(filename[i], ".vsi") )) {
		open(path+filename[i]);
		getDimensions(width, height, channels, slices, frames);
		orig_name = File.nameWithoutExtension;
		name = replace(orig_name, "_", "");
		name = replace(name, " ", ""); // Remove spaces from the filename
		j=i+1;
		origIm=getImageID();
		print("processing file "+filename[i]);
		rename("image"+i);
		run("Split Channels");
		for (k=1;k<=channels;k++) {
			image_name="C"+k+"-image"+i;
			selectImage(image_name);
			run("Nrrd ... ", "nrrd=["+resDir + name +"_"+j+"_0"+k+".nrrd]");

			close();
			
		}
	}run("Close All");
	
}
run("Close All");
run("Collect Garbage");

}
setBatchMode(false);
print("done");
