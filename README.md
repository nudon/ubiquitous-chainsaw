# ubiquitous-chainsaw
Small repo for some audio processing code. 
Github suggested the name and it sounded cool.

Currently playing around with sound library basics
next goal is to be able to parse audio files and extract groups of sounds, such as a a series of notes
Upon success of that goal, come up with some metric for defining and comparing groups of sounds
Upon success of that goal, attempt to reconstruct songs soley using "sound-fonts" extracted from other songs

If all that works, it would be cool to toy around with generating new song compositions from some corpus, either using a large soundbank of samples or attempting to synthesize everything from sinewaves and lots of effects. 
But that's pretty far off. 

libraries used in this.
kiss_fft, a nice library for computing fourier transforms, found at
https://github.com/mborgerding/kissfft

Second, some code from Iowa Hills Software, found for "free" at
http://www.iowahills.com/ 
more specifically, the parts for creation of FIR filters 
http://www.iowahills.com/A7ExampleCodePage.html

Original code is for both libraries are respectively contained in "kissfft" and "Iowa_Hills_Source_Code_Kit" folders
Made very little adjustments to makefiles and such to get the code to compile and link
