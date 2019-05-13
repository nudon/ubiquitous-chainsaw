CC = g++
MATH_LINK = -lm

EIGEN_CONFIG = `pkg-config --cflags eigen3`
FT_CONFIG = `pkg-config --cflags --libs freetype2`

STK_LINK = -lstk
STK_PATH = -I/usr/include/stk/
STK_COMPILE_FLAGS = $(STK_PATH) -D__LITTLE_ENDIAN__

KISSFFT_DIR = ./kissfft/tools/
KISSFFT_PATH = "-I"$(KISSFFT_DIR)
KISSFFT_TARGET = kiss_fftr_objects
KISSFFT_OBJECTS = $(KISSFFT_DIR)kiss_fftr.o $(KISSFFT_DIR)kiss_fft.o 
KISSFFT_LIB = -L$(KISSFFT_DIR)  -lfastconvr_float

IHILLS_DIR = ./Iowa_Hills_Source_Code_Kit/
IHILLS_PATH = "-I"$(IHILLS_DIR)
IHILLS_FILTER_TARGET = filter_objects
IHILLS_OBJECTS = $(IHILLS_DIR)FIRFilterCode.o $(IHILLS_DIR)FFTCode.o

PNGWRITER_PLEASE =  -lpng -L./png/ -lPNGwriter

MY_OBJECTS = driver.o Chunk.o ChunkCompare.o ChunkMatch.o ChunkStats.o Song.o SongEmbedder.o ChunkFilter.o  util.o eigen_to_image.o

TOT_PATH = $(IHILLS_PATH) $(KISSFFT_PATH) $(STK_PATH)
TOT_LINK = $(MATH_LINK) $(STK_LINK) $(PNGWRITER_PLEASE)

COMP_FLAGS = -Wall -g -O2 $(EIGEN_CONFIG) 
TARGET = test
ALLOBJ = $(MY_OBJECTS) $(KISSFFT_OBJECTS) $(IHILLS_OBJECTS)

all: $(ALLOBJ) dumby_eigen_to_image
	$(CC) $(COMP_FLAGS) $(ALLOBJ) $(TOT_LINK) -o $(TARGET)

#alternative build which outputs picture to build, required a wonky configuration on my end to work
#also seems to only output picture to root project dir
picture: $(ALLOBJ) working_eigen_to_image
	$(CC) $(COMP_FLAGS) $(ALLOBJ) $(TOT_LINK) -o $(TARGET)


test.o: test.cpp
	$(CC) $(COMP_FLAGS) -c test.cpp $(TOT_PATH)

dumby_eigen_to_image: eigen_to_image.cpp
	$(CC) $(COMP_FLAGS) -c $< $(TOT_PATH)

working_eigen_to_image: eigen_to_image.cpp
	echo "expect picture output to clog the project dir"
	$(CC) $(COMP_FLAGS)  -D PICTURE -c $< $(TOT_PATH)

$(KISSFFT_OBJECTS):
	cd $(KISSFFT_DIR) && make $(KISSFFT_TARGET)

$(IHILLS_OBJECTS):
	cd $(IHILLS_DIR) && make $(IHILLS_FILTER_TARGET)

%.o : %.cpp
	$(CC) $(COMP_FLAGS) -c $< $(TOT_PATH)
clean:
	cd $(KISSFFT_DIR) && make clean
	cd $(IHILLS_DIR) && make clean                  	
	rm -f *.o $(TARGET)


