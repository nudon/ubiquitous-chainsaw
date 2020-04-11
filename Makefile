CC = g++
MATH_LINK = -lm

EIGEN_CONFIG = `pkg-config --cflags eigen3`
FT_CONFIG = `pkg-config --cflags --libs freetype2`

SRC_DIR = ./src/

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

MY_OBJECTS = driver.o Chunk.o ChunkCompare.o ChunkMatch.o ChunkStats.o Song.o SongEmbedder.o ChunkFilter.o  util.o 
EXT_OBJECTS = $(KISSFFT_OBJECTS) $(IHILLS_OBJECTS)
ALLOBJ = $(MY_OBJECTS) $(EXT_OBJECTS)

TOT_PATH = $(IHILLS_PATH) $(KISSFFT_PATH) $(STK_PATH)
TOT_LINK = $(MATH_LINK) $(STK_LINK) $(PNGWRITER_PLEASE)

COMP_FLAGS = -Wall -g -O1 $(EIGEN_CONFIG) 
TARGET = test


all: $(addprefix $(SRC_DIR), $(MY_OBJECTS) dumby_eigen_to_image.o) $(EXT_OBJECTS)
#	cd $(SRC_DIR) && $(CC) $(COMP_FLAGS) $(ALLOBJ) $(TOT_LINK) -o ../$(TARGET)
	$(CC) $(COMP_FLAGS) $^ $(TOT_LINK) -o $(TARGET)

#alternative build which outputs picture to build, required a wonky configuration on my end to work
#also seems to only output picture to root project dir
picture: $(addprefix $(SRC_DIR), $(MY_OBJECTS) working_eigen_to_image.o ) $(EXT_OBJECTS)
	$(CC) $(COMP_FLAGS) $^ $(TOT_LINK) -o $(TARGET)

$(SRC_DIR)dumby_eigen_to_image.o: $(SRC_DIR)eigen_to_image.cpp
	$(CC) $(COMP_FLAGS) -c -o $@ $< $(TOT_PATH)

$(SRC_DIR)working_eigen_to_image.o: $(SRC_DIR)eigen_to_image.cpp
	echo "expect picture output to clog the project dir"
	$(CC) $(COMP_FLAGS)  -D PICTURE -c -o $@ $< $(TOT_PATH)

$(KISSFFT_OBJECTS):
	cd $(KISSFFT_DIR) && make $(KISSFFT_TARGET)

$(IHILLS_OBJECTS):
	cd $(IHILLS_DIR) && make $(IHILLS_FILTER_TARGET)

$(SRCDIR)%.o : $(SRCDIR)%.cpp
	$(CC) $(COMP_FLAGS) $(TOT_PATH) -c -o $@ $< 
clean:
	cd $(KISSFFT_DIR) && make clean
	cd $(IHILLS_DIR) && make clean                  	
	cd $(SRC_DIR) && rm -f *.o $(TARGET)


