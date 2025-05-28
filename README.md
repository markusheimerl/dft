# dft
A minimal dft implementation

```
sudo apt update && sudo apt install time libopenblas-dev liblapacke-dev
git clone http://github.com/sunqm/libcint.git
cd libcint/
mkdir build; cd build
make -DCMAKE_INSTALL_PREFIX=/usr/ ..
make install
cd ..
git clone https://github.com/markusheimerl/dft.git
cd dft/
make run
```