printf "Downloading libsnarks"
git clone https://github.com/scipr-lab/libsnark.git
cd ./libsnark
git checkout staging

cd ./depends

rm -rf ./ate-pairing
rm -rf ./gtest
rm -rf ./libff
rm -rf ./libfqfft
rm -rf ./libsnark-supercop
rm -rf ./xbyak

printf "\nDownloading ate-pairing"
git clone https://github.com/herumi/ate-pairing.git
cd ./ate-pairing
git checkout e698901
cd ..

printf "\nDownloading gtest"
git clone https://github.com/google/googletest.git
mv googletest gtest
cd ./gtest
git checkout 3a4cf1a
cd ..

printf "\nDownloading libff"
git clone https://github.com/scipr-lab/libff.git
cd ./libff
git checkout 674e437
cd ..

printf "\nDownloading libsnark-supercop"
git clone https://github.com/mbbarbosa/libsnark-supercop.git
cd ./libsnark-supercop
git checkout b04a0ea
cd ..

printf "\nDownloading xbyak"
git clone https://github.com/herumi/xbyak.git
cd ./xbyak
git checkout f0a8f7f
cd ..

printf "\nDownloading libfqfft"
git clone https://github.com/scipr-lab/libfqfft.git
cd ./libfqfft
git checkout 7d460ca
printf "\nInstalling libfqfft"
cd ./depends
rm -rf ./ate-pairing
rm -rf ./gtest
rm -rf ./libff
rm -rf ./xbyak
cp -a ../../ate-pairing/ .
cp -a ../../gtest/ .
cp -a ../../libff/ .
cp -a ../../xbyak/ .
cd ..
mkdir build && cd build && cmake ..
make -j16 && make check 
sudo make install
cd .. && cd .. && cd ..

printf "\nInstalling libsnark"
mkdir build && cd build && cmake ..
make -j16 && make check && make doc
sudo make install
cd .. && cd ..

printf "\nBuilding hash function project"
mkdir -p build
mkdir -p bin
mkdir -p lib
make -j16
