MAKE = make
TARGET = life
SOURCE = Main.cpp

default:
mpicxx -std=c++14 -o $(TARGET) $(SOURCE