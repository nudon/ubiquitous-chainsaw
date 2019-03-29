CC = g++
MATH_LINK = -lm

EIGEN_CONFIG = `pkg-config --cflags eigen3`

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


TOT_PATH = $(IHILLS_PATH) $(KISSFFT_PATH) $(STK_PATH)
TOT_LINK = $(MATH_LINK) $(STK_LINK) 

COMP_FLAGS = -Wall -g $(EIGEN_CONFIG)
TARGET = test
ALLOBJ = test.o $(KISSFFT_OBJECTS) $(IHILLS_OBJECTS)

all: $(ALLOBJ)
	$(CC) $(COMP_FLAGS) $(ALLOBJ) $(TOT_LINK) -o $(TARGET)

test.o: test.cpp
	$(CC) $(COMP_FLAGS) -c test.cpp $(TOT_PATH)

$(KISSFFT_OBJECTS):
	cd $(KISSFFT_DIR) && make $(KISSFFT_TARGET)

$(IHILLS_OBJECTS):
	cd $(IHILLS_DIR) && make $(IHILLS_FILTER_TARGET)

%.o : %.cpp
	$(CC) -c $<
clean:
	cd $(KISSFFT_DIR) && make clean
	cd $(IHILLS_DIR) && make clean                  	
	rm -f *.o $(TARGET)


