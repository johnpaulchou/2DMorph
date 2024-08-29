# 2DMorph
See https://www.sciencedirect.com/science/article/pii/S0168900214011814 for details

This needs to be built against ROOT 6.30 with a patch to the RooMomentMorphFuncND copy constructor, which needs the member _isPdfMode to be initialized

# build instructions
```mkdir build
cd build
cmake ..
make
./2DMorph
```
