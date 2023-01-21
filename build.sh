#!/bin/bash
SITE_PACKAGE_PATH=$(pip show pip|grep "Location"|awk -F ': ' '{print $2}')
if [ -z $SITE_PACKAGE_PATH ]; then
    echo "pip not found"
    exit 1
fi
pip show scikit-build>/dev/null
if [ $? -ne 0 ]; then
    pip install scikit-build
fi
echo $SITE_PACKAGE_PATH
CMAKE_MODULE_PATH=$SITE_PACKAGE_PATH/skbuild/resources/cmake
CMAKE_INSTALL_PREFIX=`pwd`
if [ ! -d $CMAKE_MODULE_PATH ]; then
    echo "cmake module $CMAKE_MODULE_PATH not existed"
    exit 1
fi
CMAKE_CXX_COMPILER=$(which x86_64-conda_cos6-linux-gnu-g++)
# CMAKE_CXX_COMPILER=`which g++`
rm -rf build && mkdir build && pushd build
cmake \
    -DCMAKE_MODULE_PATH=$CMAKE_MODULE_PATH \
    -DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} \
    ..
make -j8
if [ $? -ne 0 ]; then
    echo "make failed"
    exit 1
fi
rm diamondpy/*so
make install
if [ $? -ne 0 ]; then
    echo "make install failed"
    exit 1
fi
popd
python3 test.py