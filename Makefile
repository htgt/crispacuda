CC        :=gcc
CXXFLAGS  :=-std=c++11 -O3 -Wall
NVCC      :=nvcc
NVCCFLAGS :=-std=c++11 -O3
BUILD     :=build
SRCDIR    :=src
TARGET    :=crispacuda
CU_SRC    :=$(wildcard src/*.cu)
CPP_SRC   :=$(wildcard src/*.cpp)
CU_OBJS   :=$(patsubst $(SRCDIR)/%,$(BUILD)/%,$(CU_SRC:.cu=.o))
CPP_OBJS  :=$(patsubst $(SRCDIR)/%,$(BUILD)/%,$(CPP_SRC:.cpp=.o))
OBJECTS   :=$(CU_OBJS) $(CPP_OBJS)

all: $(TARGET)

$(BUILD)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

$(BUILD)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJECTS)
	@echo "Linking crispacuda $^"
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) $(OBJECTS)

clean:
	rm $(TARGET)
	rm -rvf $(BUILD)

