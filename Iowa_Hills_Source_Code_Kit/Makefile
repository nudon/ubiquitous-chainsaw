CC = g++
MATH_LINK = -lm
TOT_LINK = $(MATH_LINK) $(STK_LINK) 

COMP_FLAGS = -Wall -g

all: $(ALLOBJ)
	$(CC) $(COMP_FLAGS) $(ALLOBJ) $(TOT_LINK) -o $(TARGET)

filter_objects: FIRFilterCode.cpp FFTCode.cpp
	$(CC) $(COMP_FLAGS) -c $^

clean:
	rm -f *.o $(TARGET)


