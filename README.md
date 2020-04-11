# ubiquitous-chainsaw
repo for some audio processing code. 

currently ubiquitous chainsaw extracts important sounds from songs and is able to map the extracted sounds from one song to the template of important sounds found in another. 

If I'm able to, it would be cool to toy around with generating new song compositions from some corpus, either using a large soundbank of samples or attempting to synthesize everything from sinewaves and lots of effects. 
But that's pretty far off. 

libraries used in this.

stk, a useful library which provides tools for reading and filtering audio files.

https://ccrma.stanford.edu/software/stk/index.html


eigen, a useful library providing matrices

https://bitbucket.org/eigen/eigen/


kiss_fft, a nice library for computing fourier transforms, found at

https://github.com/mborgerding/kissfft

some code from Iowa Hills Software, found for "free" at

http://www.iowahills.com/

more specifically, the parts for creation of FIR filters

http://www.iowahills.com/A7ExampleCodePage.html

Original code from kiss_fft and Iowa Hills Software are respectively contained in "kissfft" and "Iowa_Hills_Source_Code_Kit" folders
Made very little adjustments to makefiles and such to get the code to compile and link
