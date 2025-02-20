Chemistry Tests
===============

Tests in this folder are used to support the porting of CAM-Chem functionality to CAM-SIMA.

To build and run the test suite, from this folder run:

```
mkdir build
cd build
cmake ..
make
make test
```

To run the tests in a Docker container, from the `CAM/test/chemistry` folder run:

```
docker build -t cam-test ../../ -f Dockerfile
docker run -it cam-test bash
make test
```