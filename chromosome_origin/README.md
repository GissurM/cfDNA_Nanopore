# Chromosomal origin of cfDNA

This is a custom analysis method that uses the header information given to each individual cfDNA thread as part of the dorado alignment process. 
To further explain each thread is given a flag by the likes of chr 1 etc. there are also some abnormal headers which i assume dorado likes using because there are a lot of them so they will be placed in a seperate column in the analysis. 
In all likelihood this will be the tallest column. There are only really 3 ways to use this code which i think are actually beneficial the first is to use this to determine biological gender as female samples consistently have a way smaller chr Y column and a chr X column double the size of male samples, 
secondly this code can be used to test the validity of your run if your run was performed on human samples you should expect to see the highest column for chr1 and then a consistend descend with ascending chromosome numbers (chrX and chrY are not counted in this)
and then finally if you used T2T as your reference you will see a small amount of mitochondrial DNA and since this code gives info on the ratio of chr cfDNA origin you will be able to determine a potentially useful % of cfDNA from mitochondria in the human who you took the sample from. 
Otherwise this code is pretty simple to use and requires a couple of R libraries to be installed the code can be run locally on any computer with decent specs and is optimized slightly for lower end hardware.
