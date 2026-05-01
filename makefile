# Компилятор и флаги
CXX = g++
CXXFLAGS =  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
LDFLAGS = -lm

TARGET = a.out
SOURCES = main.cpp Discrepancy.cpp Rotation.cpp Amethods.cpp Multiplication.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = header.h

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

# Правила компиляции
main.o: main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c main.cpp

Discrepancy.o: Discrepancy.cpp 
	$(CXX) $(CXXFLAGS) -c Discrepancy.cpp

Rotation.o: Rotation.cpp 
	$(CXX) $(CXXFLAGS) -c Rotation.cpp

Amethods.o: Amethods.cpp 
	$(CXX) $(CXXFLAGS) -c Amethods.cpp
Multiplication.o: Multiplication.cpp
	$(CXX) $(CXXFLAGS) -c Multiplication.cpp


# Утилиты
clean:
	rm -f $(OBJECTS) $(TARGET)

rebuild: clean all

run: $(TARGET)
	./$(TARGET) 100 50 5 1

debug: CXXFLAGS += -g -DDEBUG
debug: rebuild
.PHONY: Multiplication.cppall clean rebuild run debug
