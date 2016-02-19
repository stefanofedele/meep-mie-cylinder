CC = g++
ROOT_DIR = ~/research/mymeep
BUILD = build
CFLAGS = -std=c++11 -I $(ROOT_DIR) -c


all: foo

foo: $(BUILD)/main.o $(BUILD)/materials.o $(BUILD)/functions.o $(BUILD)/parameters.o
	$(CC) -std=c++11 `pkg-config --cflags meep` $(BUILD)/main.o $(BUILD)/materials.o $(BUILD)/functions.o $(BUILD)/parameters.o -o foo `pkg-config --libs meep`

$(BUILD)/main.o: main.cc  $(ROOT_DIR)/materials.cc $(ROOT_DIR)/functions.cc $(ROOT_DIR)/parameters.cc
	$(CC) $(CFLAGS) main.cc
	mv main.o $(BUILD)

$(BUILD)/materials.o: $(ROOT_DIR)/materials.cc
	$(CC) $(CFLAGS) $(ROOT_DIR)/materials.cc
	mv materials.o $(BUILD)

$(BUILD)/functions.o: $(ROOT_DIR)/functions.cc
	$(CC) $(CFLAGS) $(ROOT_DIR)/functions.cc
	mv functions.o $(BUILD)

$(BUILD)/parameters.o: $(ROOT_DIR)/parameters.cc $(ROOT_DIR)/functions.cc
	$(CC) $(CFLAGS) $(ROOT_DIR)/parameters.cc
	mv parameters.o $(BUILD)

clean:
	rm *.o
