echo "Setting up g++..."
sudo apt-get install build-essential #install gcc

sudo apt-get install zlib1g-dev #install zlib
sudo apt-get install libeigecdn3-dev #install eigen

echo "Setting up OpenGL..."
sudo apt-get install mesa-common-dev
sudo apt-get install libglu1-mesa
sudo apt-get install libxmu-dev libxi-dev
sudo apt-get install freeglut3 freeglut3-dev

echo "Installing Cmake..."
sudo apt-get install cmake

echo "Installing Boost..."
sudo apt-get install libboost-all-dev #install Boost

echo "Setting up GMP..."
sudo apt-get install libgmp3-dev #install GMP

echo "Setting up MPFR..."
sudo apt-get install libmpfr-dev #install MPFR

echo "Setting up QT4..."
sudo apt-get install qt4-default #install QT4
sudo apt-get install qtdeclarative4-dev

echo "Installing CGAL"
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.6.3/CGAL-4.6.3.tar.xz
tar -xf CGAL-4.6.3.tar.xz
cd CGAL-4.6.3 # go to CGAL directory
cmake . # configure CGAL
make # build the CGAL libraries
sudo make install # install

echo "Compiling"
cd ..
cmake .
make
